//! Integration tests for the gravity_tree_collision module group of rebound_rs.
//! Part of rebound_rs, GPL-3.0-or-later.
#![allow(non_snake_case)]
// Clippy waivers. A test/example is its own crate and does not inherit
// the crate root's waivers, so they are repeated here. Same justification:
// this code mirrors the C source's idioms, and applying clippy's
// suggestions would obscure the correspondence that makes the port
// reviewable. See rebound_rust.md section 17.
#![allow(clippy::manual_clamp)] // mirrors the C's explicit min/max tests
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::identity_op)]
#![allow(clippy::erasing_op)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::manual_range_contains)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::manual_swap)]
#![allow(clippy::manual_div_ceil)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::misrefactored_assign_op)]
#![allow(clippy::neg_multiply)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::needless_late_init)]
#![allow(clippy::while_let_loop)]
#![allow(clippy::unnecessary_cast)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::seek_from_current)]
#![allow(clippy::drop_non_drop)]
#![allow(clippy::redundant_closure)]
#![allow(clippy::approx_constant)]
#![allow(clippy::useless_vec)]
#![allow(clippy::type_complexity)]
use rebound_rs::*;

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

/// Build a particle without touching the (deliberately C-shaped) struct
/// literal at every call site.
fn mkp(m: f64, rad: f64, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> reb_particle {
    let mut q = reb_particle::default();
    q.m = m;
    q.r = rad;
    q.x = x;
    q.y = y;
    q.z = z;
    q.vx = vx;
    q.vy = vy;
    q.vz = vz;
    q
}

/// A simulation with a cubic root box of `root_size * n_root` per side.
fn box_sim(root_size: f64, nx: usize, ny: usize, nz: usize) -> reb_simulation {
    let mut r = reb_simulation_create();
    r.G = 1.0;
    r.softening = 0.0;
    r.root_size = root_size;
    r.N_root_x = nx;
    r.N_root_y = ny;
    r.N_root_z = nz;
    r.N_ghost_x = 0;
    r.N_ghost_y = 0;
    r.N_ghost_z = 0;
    r.boundary = REB_BOUNDARY::NONE;
    r.save_messages = 1; // keep stderr clean; messages are inspectable instead
    r
}

/// The exact monopole kick tree.c applies for one accepted cell:
///   `_r  = sqrt(r2 + softening^2); prefact = -G/(_r*_r*_r)*m; a += prefact*d`
/// reproduced operation-for-operation so the comparison can be bit-exact.
fn monopole(
    G: f64,
    softening: f64,
    px: f64,
    py: f64,
    pz: f64,
    m: f64,
    cx: f64,
    cy: f64,
    cz: f64,
) -> (f64, f64, f64) {
    let softening2 = softening * softening;
    let dx = px - cx;
    let dy = py - cy;
    let dz = pz - cz;
    let r2 = dx * dx + dy * dy + dz * dz;
    let _r = (r2 + softening2).sqrt();
    let prefact = -G / (_r * _r * _r) * m;
    (prefact * dx, prefact * dy, prefact * dz)
}

/// Accumulate monopole kicks in the given order, starting from 0.0 exactly as
/// `reb_gravity_tree_calculate_acceleration` does.
fn walk_sum(
    G: f64,
    softening: f64,
    p: (f64, f64, f64),
    sources: &[(f64, f64, f64, f64)],
) -> (f64, f64, f64) {
    let mut ax = 0.0f64;
    let mut ay = 0.0f64;
    let mut az = 0.0f64;
    for &(m, cx, cy, cz) in sources {
        let (tx, ty, tz) = monopole(G, softening, p.0, p.1, p.2, m, cx, cy, cz);
        ax += tx;
        ay += ty;
        az += tz;
    }
    (ax, ay, az)
}

fn accels(r: &reb_simulation) -> Vec<(f64, f64, f64)> {
    (0..r.N)
        .map(|i| (r.particles[i].ax, r.particles[i].ay, r.particles[i].az))
        .collect()
}

fn max_abs(a: &[(f64, f64, f64)]) -> f64 {
    let mut m: f64 = 0.0;
    for &(x, y, z) in a {
        m = m.max(x.abs()).max(y.abs()).max(z.abs());
    }
    m
}

/// Number of leaves in the sub-tree rooted at arena index `c`.
fn count_leaves(cells: &[reb_treecell], c: usize) -> usize {
    let cell = cells[c];
    if cell.pt >= 0 {
        return 1;
    }
    let mut n = 0;
    for o in 0..8 {
        if cell.oct[o] != REB_TREECELL_NONE {
            n += count_leaves(cells, cell.oct[o]);
        }
    }
    n
}

/// Check the structural invariants of one sub-tree:
///  * every non-leaf stores `pt == -(number of particles below it)`
///  * every child has half the parent's width and is centred a quarter of
///    the parent's width away along each axis, with the sign given by the
///    octant bits
///  * every leaf's particle lies inside its own cell
fn check_subtree(cells: &[reb_treecell], particles: &[reb_particle], c: usize) {
    let cell = cells[c];
    if cell.pt >= 0 {
        let p = particles[cell.pt as usize];
        let h = cell.w / 2.0;
        assert!(
            (p.x - cell.x).abs() <= h && (p.y - cell.y).abs() <= h && (p.z - cell.z).abs() <= h,
            "leaf cell {} (centre {:?} w {}) does not contain its particle {} at ({}, {}, {})",
            c,
            (cell.x, cell.y, cell.z),
            cell.w,
            cell.pt,
            p.x,
            p.y,
            p.z
        );
        return;
    }
    let n = count_leaves(cells, c);
    assert_eq!(
        cell.pt,
        -(n as i32),
        "non-leaf cell {} stores pt={} but holds {} particles",
        c,
        cell.pt,
        n
    );
    for o in 0..8 {
        let d = cell.oct[o];
        if d == REB_TREECELL_NONE {
            continue;
        }
        let child = cells[d];
        assert_eq!(
            child.w.to_bits(),
            (cell.w / 2.0).to_bits(),
            "child width of cell {} octant {} is {} not {}",
            c,
            o,
            child.w,
            cell.w / 2.0
        );
        let sx = if (o >> 0) % 2 == 0 { 1.0 } else { -1.0 };
        let sy = if (o >> 1) % 2 == 0 { 1.0 } else { -1.0 };
        let sz = if (o >> 2) % 2 == 0 { 1.0 } else { -1.0 };
        assert_eq!(
            child.x.to_bits(),
            (cell.x + child.w / 2.0 * sx).to_bits(),
            "child centre x of cell {} octant {}",
            c,
            o
        );
        assert_eq!(
            child.y.to_bits(),
            (cell.y + child.w / 2.0 * sy).to_bits(),
            "child centre y of cell {} octant {}",
            c,
            o
        );
        assert_eq!(
            child.z.to_bits(),
            (cell.z + child.w / 2.0 * sz).to_bits(),
            "child centre z of cell {} octant {}",
            c,
            o
        );
        check_subtree(cells, particles, d);
    }
}

/// Replay `reb_collision_search`'s post-search Fisher-Yates-style shuffle,
/// which draws exactly `n` values from the simulation RNG.
fn replay_shuffle(seed: &mut u32, order: &mut Vec<(usize, usize)>) {
    let n = order.len();
    for i in 0..n {
        let new = (rand_r(seed) as usize) % n;
        order.swap(i, new);
    }
}

fn found_pairs(r: &reb_simulation) -> Vec<(usize, usize)> {
    (0..r.N_collisions)
        .map(|i| (r.collisions[i].p1, r.collisions[i].p2))
        .collect()
}

fn normalized_sorted(pairs: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut v: Vec<(usize, usize)> = pairs
        .iter()
        .map(|&(a, b)| if a <= b { (a, b) } else { (b, a) })
        .collect();
    v.sort();
    v
}

// ---------------------------------------------------------------------------
// octree construction
// ---------------------------------------------------------------------------

/// A lone particle produces one root cell whose geometry is fixed by
/// `root_size` alone: centre at the root-box centre, width == root_size,
/// `pt` == the particle index.
#[test]
fn tree_single_particle_root_cell_is_the_root_box() {
    let mut r = box_sim(10.0, 1, 1, 1);
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 1.25, -2.5, 3.75, 0., 0., 0.));
    tree::reb_tree_construct(&mut r);

    assert_eq!(r.tree_cells.len(), 1, "one particle must create exactly one cell");
    assert_eq!(r.tree_root.len(), 1, "1x1x1 root grid must have one root slot");
    assert_eq!(r.tree_root[0], 0, "the only cell must be the root of box 0");
    let c = r.tree_cells[0];
    assert_eq!(c.pt, 0, "root cell must be a leaf holding particle 0");
    assert_eq!(c.w.to_bits(), 10.0f64.to_bits(), "root cell width must equal root_size");
    assert_eq!(c.x.to_bits(), 0.0f64.to_bits(), "root cell centre x");
    assert_eq!(c.y.to_bits(), 0.0f64.to_bits(), "root cell centre y");
    assert_eq!(c.z.to_bits(), 0.0f64.to_bits(), "root cell centre z");
    for o in 0..8 {
        assert_eq!(c.oct[o], REB_TREECELL_NONE, "leaf must have no daughter in octant {}", o);
    }
}

/// Two particles in different octants split the root exactly once. The
/// arena order is fixed by tree.c: the *resident* particle's new cell is
/// allocated before the incoming particle's.
#[test]
fn tree_two_particles_split_root_in_resident_first_order() {
    let mut r = box_sim(10.0, 1, 1, 1);
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 1.0, 1.0, 1.0, 0., 0., 0.)); // octant 0
    reb_simulation_add(&mut r, mkp(1.0, 0.0, -1.0, -1.0, -1.0, 0., 0., 0.)); // octant 7
    tree::reb_tree_construct(&mut r);

    assert_eq!(r.tree_cells.len(), 3, "root + two daughters");
    assert_eq!(r.tree_root[0], 0, "root must stay at arena index 0");
    assert_eq!(r.tree_cells[0].pt, -2, "root pt must encode -(2 particles)");
    // Daughter for the resident particle 0 (octant 0) comes first.
    assert_eq!(r.tree_cells[1].pt, 0, "arena cell 1 must hold the resident particle 0");
    assert_eq!(r.tree_cells[2].pt, 1, "arena cell 2 must hold the incoming particle 1");
    assert_eq!(r.tree_cells[0].oct[0], 1, "root octant 0 -> cell 1");
    assert_eq!(r.tree_cells[0].oct[7], 2, "root octant 7 -> cell 2");
    for (idx, sign) in [(1usize, 1.0f64), (2usize, -1.0f64)] {
        let c = r.tree_cells[idx];
        assert_eq!(c.w.to_bits(), 5.0f64.to_bits(), "daughter width must be root_size/2");
        assert_eq!(c.x.to_bits(), (2.5f64 * sign).to_bits(), "daughter {} centre x", idx);
        assert_eq!(c.y.to_bits(), (2.5f64 * sign).to_bits(), "daughter {} centre y", idx);
        assert_eq!(c.z.to_bits(), (2.5f64 * sign).to_bits(), "daughter {} centre z", idx);
    }
}

/// Three particles on the coordinate axes all fall in root octant 0 and
/// force two extra levels of subdivision. Every cell index, width, centre
/// and `pt` value below is derived by hand from tree.c's recursion.
#[test]
fn tree_deep_recursion_arena_layout_is_exact() {
    let mut r = box_sim(8.0, 1, 1, 1);
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 1.0, 0.0, 0.0, 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 0.0, 1.0, 0.0, 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 0.0, 0.0, 1.0, 0., 0., 0.));
    tree::reb_tree_construct(&mut r);

    assert_eq!(r.tree_cells.len(), 6, "root + 2 internal + 3 leaves");
    // pt of the three internal cells is -(particles below)
    assert_eq!(r.tree_cells[0].pt, -3, "root holds 3 particles");
    assert_eq!(r.tree_cells[1].pt, -3, "octant-0 cell holds 3 particles");
    assert_eq!(r.tree_cells[2].pt, -3, "octant-7-of-cell-1 holds 3 particles");
    assert_eq!(r.tree_cells[3].pt, 0, "leaf 3 holds particle 0");
    assert_eq!(r.tree_cells[4].pt, 1, "leaf 4 holds particle 1");
    assert_eq!(r.tree_cells[5].pt, 2, "leaf 5 holds particle 2");

    let expect = [
        (0usize, 8.0f64, (0.0f64, 0.0f64, 0.0f64)),
        (1, 4.0, (2.0, 2.0, 2.0)),
        (2, 2.0, (1.0, 1.0, 1.0)),
        (3, 1.0, (1.5, 0.5, 0.5)),
        (4, 1.0, (0.5, 1.5, 0.5)),
        (5, 1.0, (0.5, 0.5, 1.5)),
    ];
    for (idx, w, (cx, cy, cz)) in expect {
        let c = r.tree_cells[idx];
        assert_eq!(c.w.to_bits(), w.to_bits(), "cell {} width", idx);
        assert_eq!(c.x.to_bits(), cx.to_bits(), "cell {} centre x", idx);
        assert_eq!(c.y.to_bits(), cy.to_bits(), "cell {} centre y", idx);
        assert_eq!(c.z.to_bits(), cz.to_bits(), "cell {} centre z", idx);
    }
    // daughter links of the deepest internal cell, in octant order
    assert_eq!(r.tree_cells[2].oct[3], 5, "octant 3 of cell 2 -> particle 2");
    assert_eq!(r.tree_cells[2].oct[5], 4, "octant 5 of cell 2 -> particle 1");
    assert_eq!(r.tree_cells[2].oct[6], 3, "octant 6 of cell 2 -> particle 0");
}

/// Structural invariants on a larger, pseudo-random cloud: leaf count,
/// `pt` bookkeeping, cell geometry and leaf containment.
#[test]
fn tree_structural_invariants_hold_for_a_random_cloud() {
    let mut r = box_sim(8.0, 1, 1, 1);
    let mut seed: u32 = 9001;
    let n = 40usize;
    for _ in 0..n {
        let x = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let y = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let z = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let m = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) + 0.5;
        reb_simulation_add(&mut r, mkp(m, 0.0, x, y, z, 0., 0., 0.));
    }
    tree::reb_tree_construct(&mut r);
    assert!(r.messages.is_empty(), "tree construction reported: {:?}", r.messages);

    let root = r.tree_root[0];
    assert!(root != REB_TREECELL_NONE, "root box must have a tree");
    assert_eq!(count_leaves(&r.tree_cells, root), n, "one leaf per particle");
    check_subtree(&r.tree_cells, &r.particles, root);

    // every arena cell must be a live cell reachable from the root, so the
    // number of leaf cells in the arena equals N and their pt values are a
    // permutation of 0..N
    let mut pts: Vec<i32> = r.tree_cells.iter().filter(|c| c.pt >= 0).map(|c| c.pt).collect();
    pts.sort();
    assert_eq!(pts, (0..n as i32).collect::<Vec<i32>>(), "leaf pt values must be 0..N");
}

/// `reb_get_rootbox_for_particle` and the root-cell placement agree: with
/// a 2x2x2 root grid, particle b sits at the centre of root box b.
#[test]
fn tree_multiple_root_boxes_index_and_centre_mapping() {
    let mut r = box_sim(1.0, 2, 2, 2);
    for b in 0..8usize {
        let i = (b % 2) as f64;
        let j = ((b / 2) % 2) as f64;
        let k = (b / 4) as f64;
        let p = mkp(1.0, 0.0, -0.5 + i, -0.5 + j, -0.5 + k, 0., 0., 0.);
        assert_eq!(
            reb_get_rootbox_for_particle(&r, p) as usize,
            b,
            "root box index for particle at ({}, {}, {})",
            p.x,
            p.y,
            p.z
        );
        reb_simulation_add(&mut r, p);
    }
    tree::reb_tree_construct(&mut r);

    assert_eq!(r.tree_root.len(), 8, "2x2x2 root grid");
    assert_eq!(r.tree_cells.len(), 8, "one leaf cell per root box");
    for b in 0..8usize {
        assert_eq!(r.tree_root[b], b, "root box {} must own arena cell {}", b, b);
        let c = r.tree_cells[b];
        assert_eq!(c.pt as usize, b, "root box {} must hold particle {}", b, b);
        assert_eq!(c.w.to_bits(), 1.0f64.to_bits(), "root cell {} width", b);
        assert_eq!(c.x.to_bits(), r.particles[b].x.to_bits(), "root cell {} centre x", b);
        assert_eq!(c.y.to_bits(), r.particles[b].y.to_bits(), "root cell {} centre y", b);
        assert_eq!(c.z.to_bits(), r.particles[b].z.to_bits(), "root cell {} centre z", b);
    }
}

/// The hierarchical mass/centre-of-mass pass. For the axis configuration
/// the intermediate values are exact small integers, so this pins the
/// summation down completely.
#[test]
fn tree_gravity_data_center_of_mass_is_exact_for_unit_masses() {
    let mut r = box_sim(16.0, 1, 1, 1);
    reb_simulation_add(&mut r, mkp(1.0, 0.0, -6.0, -6.0, -6.0, 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 6.0, 6.0, 6.0, 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 2.0, 2.0, 2.0, 0., 0., 0.));
    tree::reb_tree_construct(&mut r);
    tree::reb_tree_calculate_gravity_data(&mut r);

    // arena: 0 root, 1 leaf(p0), 2 internal, 3 leaf(p1), 4 leaf(p2)
    assert_eq!(r.tree_cells.len(), 5, "root + leaf(p0) + internal + leaf(p1) + leaf(p2)");
    let inner = r.tree_cells[2];
    assert_eq!(inner.m.to_bits(), 2.0f64.to_bits(), "inner cell total mass");
    assert_eq!(inner.mx.to_bits(), 4.0f64.to_bits(), "inner cell com x = (6+2)/2");
    assert_eq!(inner.my.to_bits(), 4.0f64.to_bits(), "inner cell com y = (6+2)/2");
    assert_eq!(inner.mz.to_bits(), 4.0f64.to_bits(), "inner cell com z = (6+2)/2");

    let root = r.tree_cells[0];
    assert_eq!(root.m.to_bits(), 3.0f64.to_bits(), "root total mass");
    let expect = 2.0f64 / 3.0;
    assert_eq!(root.mx.to_bits(), expect.to_bits(), "root com x = (6+2-6)/3");
    assert_eq!(root.my.to_bits(), expect.to_bits(), "root com y = (6+2-6)/3");
    assert_eq!(root.mz.to_bits(), expect.to_bits(), "root com z = (6+2-6)/3");

    // a leaf carries the particle's own position, never the cell centre
    assert_eq!(r.tree_cells[1].mx.to_bits(), (-6.0f64).to_bits(), "leaf com is the particle");
}

/// For a random cloud the tree's root monopole must reproduce the total
/// mass exactly and the brute-force centre of mass to round-off.
#[test]
fn tree_root_monopole_matches_brute_force_center_of_mass() {
    let mut r = box_sim(8.0, 1, 1, 1);
    let mut seed: u32 = 24601;
    let n = 25usize;
    for _ in 0..n {
        let x = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let y = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let z = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let m = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) + 0.5;
        reb_simulation_add(&mut r, mkp(m, 0.0, x, y, z, 0., 0., 0.));
    }
    tree::reb_tree_construct(&mut r);
    tree::reb_tree_calculate_gravity_data(&mut r);
    let root = r.tree_cells[r.tree_root[0]];

    let mut M = 0.0f64;
    let mut cx = 0.0f64;
    let mut cy = 0.0f64;
    let mut cz = 0.0f64;
    for i in 0..n {
        let p = r.particles[i];
        M += p.m;
        cx += p.m * p.x;
        cy += p.m * p.y;
        cz += p.m * p.z;
    }
    cx /= M;
    cy /= M;
    cz /= M;

    assert!(
        (root.m - M).abs() <= 1e-13 * M,
        "tree root mass {} vs brute force {}",
        root.m,
        M
    );
    assert!(
        (root.mx - cx).abs() < 1e-12 && (root.my - cy).abs() < 1e-12 && (root.mz - cz).abs() < 1e-12,
        "tree root com ({}, {}, {}) vs brute force ({}, {}, {})",
        root.mx,
        root.my,
        root.mz,
        cx,
        cy,
        cz
    );
}

/// `reb_tree_delete` really releases the arena, and the gravity module
/// leaves no tree behind (so the next step rebuilds from scratch).
#[test]
fn tree_is_deleted_after_a_gravity_evaluation() {
    let mut r = box_sim(8.0, 2, 1, 1);
    reb_simulation_add(&mut r, mkp(1.0, 0.0, -3.0, 0.5, 0.25, 0., 0., 0.));
    reb_simulation_add(&mut r, mkp(2.0, 0.0, 3.0, -0.5, 0.75, 0., 0., 0.));
    r.gravity = REB_GRAVITY::TREE;
    r.opening_angle2 = 0.0;
    reb_gravity_tree_calculate_acceleration(&mut r);

    assert!(r.tree_cells.is_empty(), "cell arena must be empty after reb_tree_delete");
    assert_eq!(r.tree_root.len(), 2, "root slot vector is kept, not freed");
    for i in 0..r.tree_root.len() {
        assert_eq!(r.tree_root[i], REB_TREECELL_NONE, "root slot {} must be cleared", i);
    }
    // and the accelerations are still the two-body ones
    let a = accels(&r);
    assert!(a[0].0 > 0.0, "particle 0 must be pulled toward +x, got ax={}", a[0].0);
    assert!(a[1].0 < 0.0, "particle 1 must be pulled toward -x, got ax={}", a[1].0);
}

// ---------------------------------------------------------------------------
// tree walk / opening angle
// ---------------------------------------------------------------------------

/// Five particles, one per root octant, so the walk for particle 1 visits
/// leaves strictly in octant order 0, 2, 4, 7 (octant 1 is itself and is
/// skipped). Reproducing that accumulation order gives bit-identical
/// accelerations; reversing it does not.
#[test]
fn tree_walk_visits_daughters_in_octant_order() {
    let mut r = box_sim(8.0, 1, 1, 1);
    let bodies = [
        (1.0f64, 1.0f64, 1.5f64, 2.0f64),   // p0 -> octant 0
        (0.7, -1.5, 1.0, 1.5),              // p1 -> octant 1 (the test particle)
        (1.3, 2.0, -1.0, 1.0),              // p2 -> octant 2
        (0.3, 1.5, 2.0, -1.5),              // p3 -> octant 4
        (2.1, -2.0, -1.5, -1.0),            // p4 -> octant 7
    ];
    for &(m, x, y, z) in bodies.iter() {
        reb_simulation_add(&mut r, mkp(m, 0.0, x, y, z, 0., 0., 0.));
    }
    tree::reb_tree_construct(&mut r);
    assert_eq!(r.tree_cells.len(), 6, "root + five leaves, one per occupied octant");
    for k in 1..6 {
        assert_eq!(r.tree_cells[k].pt as usize, k - 1, "arena cell {} holds particle {}", k, k - 1);
    }
    assert_eq!(r.tree_cells[0].pt, -5, "root pt encodes -(5 particles)");
    for (o, cell) in [(0usize, 1usize), (1, 2), (2, 3), (4, 4), (7, 5)] {
        assert_eq!(r.tree_cells[0].oct[o], cell, "root octant {} -> arena cell {}", o, cell);
    }
    tree::reb_tree_delete(&mut r);

    r.opening_angle2 = 0.0; // never accept a cell: descend to every leaf
    reb_gravity_tree_calculate_acceleration(&mut r);

    let p1 = (bodies[1].1, bodies[1].2, bodies[1].3);
    let in_order = [
        (bodies[0].0, bodies[0].1, bodies[0].2, bodies[0].3),
        (bodies[2].0, bodies[2].1, bodies[2].2, bodies[2].3),
        (bodies[3].0, bodies[3].1, bodies[3].2, bodies[3].3),
        (bodies[4].0, bodies[4].1, bodies[4].2, bodies[4].3),
    ];
    let (ax, ay, az) = walk_sum(r.G, r.softening, p1, &in_order);
    assert_eq!(r.particles[1].ax.to_bits(), ax.to_bits(), "ax of particle 1 in octant order");
    assert_eq!(r.particles[1].ay.to_bits(), ay.to_bits(), "ay of particle 1 in octant order");
    assert_eq!(r.particles[1].az.to_bits(), az.to_bits(), "az of particle 1 in octant order");

    // the assertion above is genuinely order sensitive
    let mut reversed = in_order;
    reversed.reverse();
    let (rx, ry, rz) = walk_sum(r.G, r.softening, p1, &reversed);
    assert!(
        rx.to_bits() != ax.to_bits() || ry.to_bits() != ay.to_bits() || rz.to_bits() != az.to_bits(),
        "the four-term sum must not be order invariant, otherwise the octant-order check is vacuous"
    );
}

/// The opening-angle criterion `w*w > opening_angle2 * r2` selects a level
/// of the tree. Three particles give a two-level tree whose accepted
/// monopoles are exact numbers, so all three regimes can be pinned down.
#[test]
fn opening_angle_criterion_selects_the_tree_level() {
    let mut r = box_sim(16.0, 1, 1, 1);
    reb_simulation_add(&mut r, mkp(1.0, 0.0, -6.0, -6.0, -6.0, 0., 0., 0.)); // p0, octant 7
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 6.0, 6.0, 6.0, 0., 0., 0.)); // p1, octant 0
    reb_simulation_add(&mut r, mkp(1.0, 0.0, 2.0, 2.0, 2.0, 0., 0., 0.)); // p2, octant 0
    r.gravity = REB_GRAVITY::TREE;

    // read the root monopole the walk will see
    tree::reb_tree_construct(&mut r);
    tree::reb_tree_calculate_gravity_data(&mut r);
    let root = r.tree_cells[r.tree_root[0]];
    let (W, M, cx, cy, cz) = (root.w, root.m, root.mx, root.my, root.mz);
    tree::reb_tree_delete(&mut r);

    let p0 = (-6.0f64, -6.0f64, -6.0f64);
    let dx = p0.0 - cx;
    let dy = p0.1 - cy;
    let dz = p0.2 - cz;
    let r2_root = dx * dx + dy * dy + dz * dz;
    let crit = W * W / r2_root; // opening_angle2 at which the root flips

    // (a) well above the critical angle: the root itself is accepted, and
    // the test particle is (famously) pulled by its own monopole share.
    r.opening_angle2 = crit * 2.0;
    reb_gravity_tree_calculate_acceleration(&mut r);
    let accepted = monopole(r.G, r.softening, p0.0, p0.1, p0.2, M, cx, cy, cz);
    assert_eq!(r.particles[0].ax.to_bits(), accepted.0.to_bits(), "root-monopole ax");
    assert_eq!(r.particles[0].ay.to_bits(), accepted.1.to_bits(), "root-monopole ay");
    assert_eq!(r.particles[0].az.to_bits(), accepted.2.to_bits(), "root-monopole az");

    // (b) just below it: the root opens, its octant-0 daughter (mass 2 at
    // (4,4,4)) is still accepted, and p0's own leaf is skipped.
    r.opening_angle2 = crit * 0.5;
    reb_gravity_tree_calculate_acceleration(&mut r);
    let one_level = monopole(r.G, r.softening, p0.0, p0.1, p0.2, 2.0, 4.0, 4.0, 4.0);
    assert_eq!(r.particles[0].ax.to_bits(), one_level.0.to_bits(), "one-level-down ax");
    assert_eq!(r.particles[0].ay.to_bits(), one_level.1.to_bits(), "one-level-down ay");
    assert_eq!(r.particles[0].az.to_bits(), one_level.2.to_bits(), "one-level-down az");

    // (c) opening_angle2 == 0: descend to every leaf, in octant order
    // (octant 0 of the daughter is p1, octant 7 is p2).
    r.opening_angle2 = 0.0;
    reb_gravity_tree_calculate_acceleration(&mut r);
    let leaves = walk_sum(
        r.G,
        r.softening,
        p0,
        &[(1.0, 6.0, 6.0, 6.0), (1.0, 2.0, 2.0, 2.0)],
    );
    assert_eq!(r.particles[0].ax.to_bits(), leaves.0.to_bits(), "exact leaf-sum ax");
    assert_eq!(r.particles[0].ay.to_bits(), leaves.1.to_bits(), "exact leaf-sum ay");
    assert_eq!(r.particles[0].az.to_bits(), leaves.2.to_bits(), "exact leaf-sum az");

    // the three regimes really are different
    assert!(
        (accepted.0 - one_level.0).abs() > 0.1 * accepted.0.abs(),
        "accepting the root must differ from opening it: {} vs {}",
        accepted.0,
        one_level.0
    );
    assert!(
        (one_level.0 - leaves.0).abs() > 1e-4 * leaves.0.abs(),
        "one level down must differ from the exact leaf sum: {} vs {}",
        one_level.0,
        leaves.0
    );
}

/// Larger opening angles cost accuracy monotonically: relative to the
/// exact direct sum, the tree error must not shrink as theta^2 grows.
#[test]
fn opening_angle_error_grows_with_the_angle() {
    let mut r = box_sim(8.0, 1, 1, 1);
    let mut seed: u32 = 777;
    let n = 30usize;
    for _ in 0..n {
        let x = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let y = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let z = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let m = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) + 0.5;
        reb_simulation_add(&mut r, mkp(m, 0.0, x, y, z, 0., 0., 0.));
    }
    r.gravity = REB_GRAVITY::BASIC;
    reb_gravity_basic_calculate_acceleration(&mut r);
    let exact = accels(&r);
    let scale = max_abs(&exact);
    assert!(scale > 0.0, "the reference direct sum must be non-trivial");

    let mut last = -1.0f64;
    for &oa2 in [0.0f64, 0.01, 0.09, 0.25, 0.64].iter() {
        r.opening_angle2 = oa2;
        reb_gravity_tree_calculate_acceleration(&mut r);
        let tree_a = accels(&r);
        let mut err = 0.0f64;
        for i in 0..n {
            err = err
                .max((tree_a[i].0 - exact[i].0).abs())
                .max((tree_a[i].1 - exact[i].1).abs())
                .max((tree_a[i].2 - exact[i].2).abs());
        }
        let rel = err / scale;
        if oa2 == 0.0 {
            assert!(
                rel < 1e-13,
                "opening_angle2=0 must reproduce the direct sum, relative error {}",
                rel
            );
        }
        assert!(
            rel >= last - 1e-15,
            "tree error must not decrease with opening_angle2 (oa2={} rel={} previous={})",
            oa2,
            rel,
            last
        );
        last = rel;
    }
    assert!(last > 1e-6, "the largest opening angle must show a real error, got {}", last);
}

/// With `opening_angle2 == 0` the tree walk degenerates into the exact
/// pairwise sum, including a non-zero softening length.
#[test]
fn tree_gravity_with_softening_matches_direct_sum() {
    let mut r = box_sim(8.0, 1, 1, 1);
    let bodies = [
        (1.0f64, 1.0f64, 1.5f64, 2.0f64),
        (0.7, -1.5, 1.0, 1.5),
        (1.3, 2.0, -1.0, 1.0),
        (0.3, 1.5, 2.0, -1.5),
        (2.1, -2.0, -1.5, -1.0),
    ];
    for &(m, x, y, z) in bodies.iter() {
        reb_simulation_add(&mut r, mkp(m, 0.0, x, y, z, 0., 0., 0.));
    }
    r.softening = 0.35;
    r.opening_angle2 = 0.0;

    reb_gravity_tree_calculate_acceleration(&mut r);
    let tree_a = accels(&r);
    reb_gravity_basic_calculate_acceleration(&mut r);
    let basic_a = accels(&r);

    let scale = max_abs(&basic_a);
    assert!(scale > 0.0, "reference accelerations must be non-trivial");
    for i in 0..r.N {
        assert!(
            (tree_a[i].0 - basic_a[i].0).abs() <= 1e-14 * scale
                && (tree_a[i].1 - basic_a[i].1).abs() <= 1e-14 * scale
                && (tree_a[i].2 - basic_a[i].2).abs() <= 1e-14 * scale,
            "softened tree acceleration of particle {} {:?} vs direct {:?}",
            i,
            tree_a[i],
            basic_a[i]
        );
    }
}

/// Rebuilding the tree from identical state gives bit-identical forces.
#[test]
fn tree_gravity_is_bit_reproducible() {
    let mut r = box_sim(8.0, 1, 1, 1);
    let mut seed: u32 = 5150;
    for _ in 0..12 {
        let x = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let y = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        let z = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 6.0 - 3.0;
        reb_simulation_add(&mut r, mkp(1.0, 0.0, x, y, z, 0., 0., 0.));
    }
    r.opening_angle2 = 0.4;
    reb_gravity_tree_calculate_acceleration(&mut r);
    let first = accels(&r);
    reb_gravity_tree_calculate_acceleration(&mut r);
    let second = accels(&r);
    for i in 0..r.N {
        assert_eq!(first[i].0.to_bits(), second[i].0.to_bits(), "ax of particle {} must be reproducible", i);
        assert_eq!(first[i].1.to_bits(), second[i].1.to_bits(), "ay of particle {} must be reproducible", i);
        assert_eq!(first[i].2.to_bits(), second[i].2.to_bits(), "az of particle {} must be reproducible", i);
    }
}

// ---------------------------------------------------------------------------
// ghost boxes
// ---------------------------------------------------------------------------

/// Periodic ghost boxes are pure integer translations with no velocity.
#[test]
fn periodic_ghostbox_is_an_integer_translation() {
    let mut r = box_sim(2.0, 3, 1, 5);
    r.boundary = REB_BOUNDARY::PERIODIC;
    for i in -1..=1i32 {
        for j in -1..=1i32 {
            for k in -1..=1i32 {
                let gb = reb_boundary_get_ghostbox(&r, i, j, k);
                assert_eq!(gb.x.to_bits(), (6.0 * i as f64).to_bits(), "ghostbox x for i={}", i);
                assert_eq!(gb.y.to_bits(), (2.0 * j as f64).to_bits(), "ghostbox y for j={}", j);
                assert_eq!(gb.z.to_bits(), (10.0 * k as f64).to_bits(), "ghostbox z for k={}", k);
                assert_eq!(gb.vx, 0.0, "periodic ghostbox must not move in x");
                assert_eq!(gb.vy, 0.0, "periodic ghostbox must not move in y");
                assert_eq!(gb.vz, 0.0, "periodic ghostbox must not move in z");
            }
        }
    }
}

/// Shearing-sheet ghost boxes at t=0: the azimuthal shift vanishes, but
/// the radial neighbours already carry the Keplerian shear velocity
/// -3/2 * i * OMEGA * Lx.
#[test]
fn shear_ghostbox_at_t_zero_has_no_shift_but_full_shear_velocity() {
    let mut r = box_sim(1.0, 1, 1, 1);
    r.boundary = REB_BOUNDARY::SHEAR;
    r.OMEGA = 1.0;
    r.t = 0.0;
    for i in -1..=1i32 {
        for j in -1..=1i32 {
            let gb = reb_boundary_get_ghostbox(&r, i, j, 0);
            assert_eq!(gb.x.to_bits(), (i as f64).to_bits(), "shear ghostbox x for i={}", i);
            assert_eq!(gb.y, j as f64, "shear ghostbox y must be unshifted at t=0, i={}", i);
            assert_eq!(gb.z, 0.0, "shear ghostbox z");
            assert_eq!(gb.vx, 0.0, "shear ghostbox must not drift radially");
            assert_eq!(
                gb.vy,
                -1.5 * (i as f64),
                "shear ghostbox vy must be -3/2 i OMEGA Lx, i={}",
                i
            );
            assert_eq!(gb.vz, 0.0, "shear ghostbox vz");
        }
    }
}

/// At t = 1/4 with OMEGA = Lx = Ly = 1 every intermediate value is a dyadic
/// rational, so the azimuthal offset is exactly +-3/8 and antisymmetric.
#[test]
fn shear_ghostbox_offset_at_quarter_time_is_exact() {
    let mut r = box_sim(1.0, 1, 1, 1);
    r.boundary = REB_BOUNDARY::SHEAR;
    r.OMEGA = 1.0;
    r.t = 0.25;

    let gp = reb_boundary_get_ghostbox(&r, 1, 0, 0);
    assert_eq!(gp.vy.to_bits(), (-1.5f64).to_bits(), "outer ghostbox shear velocity");
    // shift = -fmod(vy*t - Ly/2, Ly) - Ly/2 = 0.875 - 0.5 = 0.375, y = -shift
    assert_eq!(gp.y.to_bits(), (-0.375f64).to_bits(), "outer ghostbox azimuthal offset");

    let gm = reb_boundary_get_ghostbox(&r, -1, 0, 0);
    assert_eq!(gm.vy.to_bits(), 1.5f64.to_bits(), "inner ghostbox shear velocity");
    assert_eq!(gm.y.to_bits(), 0.375f64.to_bits(), "inner ghostbox azimuthal offset");

    let g0 = reb_boundary_get_ghostbox(&r, 0, 1, 0);
    assert_eq!(g0.y.to_bits(), 1.0f64.to_bits(), "co-radial ghostbox is never sheared");
}

/// Whatever the time, the shear offset only ever differs from the exact
/// drift `-vy*t` by a whole number of box lengths. That is what makes the
/// ghost boxes tile the shearing sheet.
#[test]
fn shear_ghostbox_offset_is_the_drift_modulo_the_box() {
    let mut r = box_sim(1.0, 1, 1, 1);
    r.boundary = REB_BOUNDARY::SHEAR;
    r.OMEGA = 0.9;
    let Ly = 1.0f64;
    for step in 0..40 {
        r.t = 0.137 * (step as f64);
        for i in -1..=1i32 {
            let gb = reb_boundary_get_ghostbox(&r, i, 0, 0);
            assert_eq!(
                gb.vy,
                -1.5 * (i as f64) * r.OMEGA * 1.0,
                "shear velocity at t={} for i={}",
                r.t,
                i
            );
            let shift = -gb.y; // gb.y = Ly*0 - shift
            let q = (shift + gb.vy * r.t) / Ly;
            assert!(
                (q - q.round()).abs() < 1e-10,
                "shear offset {} at t={} for i={} is not the drift {} modulo Ly (q={})",
                shift,
                r.t,
                i,
                -gb.vy * r.t,
                q
            );
            assert!(
                shift.abs() <= 1.5 * Ly,
                "shear offset {} at t={} must stay of order the box size",
                shift,
                r.t
            );
        }
    }
}

/// Tree gravity and direct gravity must agree even when the force is
/// summed over 27 sheared ghost boxes at a time when the offset is
/// non-zero: both modules call `reb_boundary_get_ghostbox`.
#[test]
fn tree_gravity_over_shear_ghost_boxes_matches_direct_sum() {
    let mut r = box_sim(1.0, 1, 1, 1);
    r.boundary = REB_BOUNDARY::SHEAR;
    r.OMEGA = 1.0;
    r.t = 0.25;
    r.N_ghost_x = 1;
    r.N_ghost_y = 1;
    r.N_ghost_z = 1;
    r.opening_angle2 = 0.0;
    let bodies = [
        (1.0f64, 0.10f64, 0.21f64, 0.05f64),
        (0.6, -0.23, 0.11, -0.17),
        (1.4, 0.31, -0.14, 0.22),
        (0.9, -0.07, -0.28, 0.13),
        (0.4, 0.19, 0.02, -0.31),
    ];
    for &(m, x, y, z) in bodies.iter() {
        reb_simulation_add(&mut r, mkp(m, 0.0, x, y, z, 0., 0., 0.));
    }
    assert_eq!(r.N, bodies.len(), "all particles must have been inside the box");

    reb_gravity_tree_calculate_acceleration(&mut r);
    let tree_a = accels(&r);
    reb_gravity_basic_calculate_acceleration(&mut r);
    let basic_a = accels(&r);

    let scale = max_abs(&basic_a);
    assert!(scale > 0.0, "sheared ghost box forces must be non-trivial");
    for i in 0..r.N {
        assert!(
            (tree_a[i].0 - basic_a[i].0).abs() <= 1e-12 * scale
                && (tree_a[i].1 - basic_a[i].1).abs() <= 1e-12 * scale
                && (tree_a[i].2 - basic_a[i].2).abs() <= 1e-12 * scale,
            "sheared-ghost tree acceleration of particle {} {:?} vs direct {:?} (scale {})",
            i,
            tree_a[i],
            basic_a[i],
            scale
        );
    }
}

// ---------------------------------------------------------------------------
// collision search
// ---------------------------------------------------------------------------

/// Three mutually overlapping, mutually approaching particles. The DIRECT
/// module scans projectiles i then targets j != i, so it records every
/// pair twice, in the order (0,1) (0,2) (1,0) (1,2) (2,0) (2,1); the
/// post-search shuffle then draws exactly N_collisions values from the
/// simulation RNG.
#[test]
fn direct_collision_search_order_and_shuffle_are_exact() {
    let mut r = box_sim(8.0, 1, 1, 1);
    r.collision = REB_COLLISION::DIRECT;
    r.rand_seed = 42;
    reb_simulation_add(&mut r, mkp(1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0));
    reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0));
    reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0));

    reb_collision_search(&mut r);

    assert_eq!(r.N_collisions, 6, "each of the three pairs must be found in both directions");
    let mut expected = vec![(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)];
    let mut seed: u32 = 42;
    replay_shuffle(&mut seed, &mut expected);
    assert_eq!(found_pairs(&r), expected, "shuffled DIRECT collision list");
    assert_eq!(r.rand_seed, seed, "the search must draw exactly N_collisions random numbers");
    assert_eq!(
        r.status, REB_STATUS_COLLISION,
        "the default resolver halts the simulation"
    );
}

/// The TREE module scans the same pairs but reaches them through the
/// octree, so its raw order follows the daughter octants (3, 5, 6 of the
/// deepest cell = particles 2, 1, 0) rather than the particle index.
#[test]
fn tree_collision_search_follows_the_octree_walk_order() {
    let mut r = box_sim(8.0, 1, 1, 1);
    r.collision = REB_COLLISION::TREE;
    r.rand_seed = 42;
    reb_simulation_add(&mut r, mkp(1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0));
    reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0));
    reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0));

    reb_collision_search(&mut r);

    assert_eq!(r.N_collisions, 6, "tree search must find the same six directed pairs");
    let mut expected = vec![(0, 2), (0, 1), (1, 2), (1, 0), (2, 1), (2, 0)];
    let mut seed: u32 = 42;
    replay_shuffle(&mut seed, &mut expected);
    assert_eq!(found_pairs(&r), expected, "shuffled TREE collision list");
    assert_eq!(r.rand_seed, seed, "the search must draw exactly N_collisions random numbers");
    assert!(r.tree_cells.is_empty(), "tree collision search must delete its tree");

    // and the walk order really is different from the index order
    let direct_order = vec![(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)];
    let mut direct_shuffled = direct_order;
    let mut s2: u32 = 42;
    replay_shuffle(&mut s2, &mut direct_shuffled);
    assert!(
        found_pairs(&r) != direct_shuffled,
        "the tree walk order must not coincide with the direct index order"
    );
}

/// On a random cloud with equal radii the tree neighbour search and the
/// O(N^2) direct search must find exactly the same set of colliding pairs,
/// which is also derivable by brute force.
#[test]
fn tree_and_direct_collision_searches_find_the_same_pairs() {
    let n = 14usize;
    let radius = 0.45f64;
    let build = || {
        let mut r = box_sim(8.0, 1, 1, 1);
        let mut seed: u32 = 20240607;
        for _ in 0..n {
            let x = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 3.0 - 1.5;
            let y = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 3.0 - 1.5;
            let z = (rand_r(&mut seed) as f64) / (REB_RAND_MAX as f64) * 3.0 - 1.5;
            // every particle falls toward the origin, so every overlapping
            // pair is automatically an approaching pair
            reb_simulation_add(&mut r, mkp(1.0, radius, x, y, z, -0.1 * x, -0.1 * y, -0.1 * z));
        }
        r
    };

    // independent brute-force expectation
    let probe = build();
    let mut expected: Vec<(usize, usize)> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let a = probe.particles[i];
            let b = probe.particles[j];
            let d2 = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
            if d2 <= (2.0 * radius) * (2.0 * radius) {
                expected.push((i, j));
            }
        }
    }
    assert!(!expected.is_empty(), "the test configuration must contain overlaps");
    assert!(
        expected.len() < n * (n - 1) / 2,
        "the test configuration must not be a single blob"
    );
    let mut expected_both: Vec<(usize, usize)> = Vec::new();
    for &(a, b) in expected.iter() {
        expected_both.push((a, b));
        expected_both.push((a, b));
    }
    expected_both.sort();

    let mut rd = build();
    rd.collision = REB_COLLISION::DIRECT;
    reb_collision_search(&mut rd);
    let mut rt = build();
    rt.collision = REB_COLLISION::TREE;
    reb_collision_search(&mut rt);

    assert_eq!(
        normalized_sorted(&found_pairs(&rd)),
        expected_both,
        "DIRECT search must find every overlapping approaching pair twice"
    );
    assert_eq!(
        normalized_sorted(&found_pairs(&rt)),
        expected_both,
        "TREE search must find exactly the same pairs as the brute-force scan"
    );
    assert_eq!(rd.N_collisions, rt.N_collisions, "tree and direct collision counts");
}

/// DIRECT requires the pair to be approaching; LINE (with no elapsed
/// timestep) only requires overlap, and it scans j > i so it reports each
/// pair once instead of twice.
#[test]
fn line_search_counts_pairs_once_and_ignores_the_approach_test() {
    let separating = |c: REB_COLLISION| {
        let mut r = box_sim(8.0, 1, 1, 1);
        r.collision = c;
        r.dt_last_done = 0.0;
        r.rand_seed = 7;
        // overlapping (distance 1, radii 1+1) but flying apart
        reb_simulation_add(&mut r, mkp(1.0, 1.0, -0.5, 0.0, 0.0, -1.0, 0.0, 0.0));
        reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0));
        reb_collision_search(&mut r);
        r.N_collisions
    };
    assert_eq!(separating(REB_COLLISION::DIRECT), 0, "DIRECT must reject a separating pair");
    assert_eq!(separating(REB_COLLISION::LINE), 1, "LINE must report the overlap once");

    let approaching = |c: REB_COLLISION| {
        let mut r = box_sim(8.0, 1, 1, 1);
        r.collision = c;
        r.dt_last_done = 0.0;
        r.rand_seed = 7;
        reb_simulation_add(&mut r, mkp(1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0));
        reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0));
        reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0));
        reb_collision_search(&mut r);
        r.N_collisions
    };
    assert_eq!(approaching(REB_COLLISION::DIRECT), 6, "DIRECT reports 3 pairs twice");
    assert_eq!(approaching(REB_COLLISION::LINE), 3, "LINE reports 3 pairs once");
}

/// No collisions means no random numbers are drawn and the status is left
/// alone; REB_COLLISION::NONE short-circuits entirely.
#[test]
fn collision_search_edge_cases_leave_the_rng_and_status_untouched() {
    // module off
    let mut r = box_sim(8.0, 1, 1, 1);
    r.collision = REB_COLLISION::NONE;
    r.rand_seed = 314159;
    reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.2, 0.0, 0.0, -1.0, 0.0, 0.0));
    reb_simulation_add(&mut r, mkp(1.0, 1.0, -0.2, 0.0, 0.0, 1.0, 0.0, 0.0));
    reb_collision_search(&mut r);
    assert_eq!(r.N_collisions, 0, "REB_COLLISION::NONE must find nothing");
    assert_eq!(r.rand_seed, 314159, "REB_COLLISION::NONE must not touch the RNG");
    assert_eq!(r.status, 0, "status must be untouched");

    // N == 0 and N == 1, direct and tree
    for c in [REB_COLLISION::DIRECT, REB_COLLISION::TREE] {
        let mut e = box_sim(8.0, 1, 1, 1);
        e.collision = c;
        e.rand_seed = 271828;
        reb_collision_search(&mut e);
        assert_eq!(e.N_collisions, 0, "empty simulation, module {:?}", c);
        assert_eq!(e.rand_seed, 271828, "empty simulation must not draw randoms, module {:?}", c);

        let mut o = box_sim(8.0, 1, 1, 1);
        o.collision = c;
        o.rand_seed = 271828;
        reb_simulation_add(&mut o, mkp(1.0, 1.0, 0.3, -0.2, 0.1, 0., 0., 0.));
        reb_collision_search(&mut o);
        assert_eq!(o.N_collisions, 0, "a particle must not collide with itself, module {:?}", c);
        assert_eq!(o.rand_seed, 271828, "single particle must not draw randoms, module {:?}", c);
        assert!(o.tree_cells.is_empty(), "tree must be released again, module {:?}", c);
    }
}

/// Shear ghost boxes let particles on opposite radial edges collide, and
/// the ghost box stored with the collision carries the shear velocity that
/// the resolver needs.
#[test]
fn shear_ghost_boxes_enable_collisions_across_the_radial_boundary() {
    let build = |n_ghost_x: i32| {
        let mut r = box_sim(1.0, 1, 1, 1);
        r.boundary = REB_BOUNDARY::SHEAR;
        r.OMEGA = 1.0;
        r.t = 0.0;
        r.collision = REB_COLLISION::DIRECT;
        r.rand_seed = 11;
        r.N_ghost_x = n_ghost_x;
        r.N_ghost_y = 1;
        r.N_ghost_z = 0;
        reb_simulation_add(&mut r, mkp(1.0, 0.1, 0.45, 0.0, 0.0, -0.5, 0.0, 0.0));
        reb_simulation_add(&mut r, mkp(1.0, 0.1, -0.45, 0.0, 0.0, -0.6, 0.0, 0.0));
        reb_collision_search(&mut r);
        r
    };

    let none = build(0);
    assert_eq!(
        none.N_collisions, 0,
        "without radial ghost boxes the two particles are 0.9 apart and cannot collide"
    );

    let ghosted = build(1);
    assert_eq!(
        ghosted.N_collisions, 2,
        "with one radial ghost ring the pair collides through both boundaries"
    );
    for i in 0..2 {
        let c = ghosted.collisions[i];
        let (want_x, want_vy) = if c.p1 == 0 { (-1.0f64, 1.5f64) } else { (1.0f64, -1.5f64) };
        assert_eq!(c.gb.x.to_bits(), want_x.to_bits(), "ghost box x stored with collision {:?}", (c.p1, c.p2));
        assert_eq!(c.gb.y.to_bits(), 0.0f64.to_bits(), "ghost box y at t=0");
        assert_eq!(
            c.gb.vy.to_bits(),
            want_vy.to_bits(),
            "ghost box shear velocity stored with collision {:?}",
            (c.p1, c.p2)
        );
    }
    let mut seen = normalized_sorted(&found_pairs(&ghosted));
    seen.dedup();
    assert_eq!(seen, vec![(0usize, 1usize)], "only the 0-1 pair may collide");
}

// ---------------------------------------------------------------------------
// collision resolution
// ---------------------------------------------------------------------------

/// Head-on elastic collision of two equal-mass spheres: the velocities are
/// exchanged exactly, and the transverse impulse cancels bit for bit.
#[test]
fn hardsphere_head_on_equal_masses_exchanges_velocities_exactly() {
    let mut r = box_sim(8.0, 1, 1, 1);
    r.collision = REB_COLLISION::DIRECT;
    r.collision_resolve = Some(reb_collision_resolve_hardsphere);
    r.rand_seed = 3;
    reb_simulation_add(&mut r, mkp(1.0, 1.5, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0));
    reb_simulation_add(&mut r, mkp(1.0, 1.5, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0));

    reb_collision_search(&mut r);
    assert_eq!(r.N_collisions, 2, "the pair is reported in both directions");

    assert_eq!(r.particles[0].vx.to_bits(), (-1.0f64).to_bits(), "particle 0 vx after bounce");
    assert_eq!(r.particles[1].vx.to_bits(), 1.0f64.to_bits(), "particle 1 vx after bounce");
    // equal masses => the transverse impulses are equal and opposite bit for bit
    assert_eq!(
        (r.particles[0].vy + r.particles[1].vy).to_bits(),
        0.0f64.to_bits(),
        "transverse momentum must cancel exactly"
    );
    let ke: f64 = (0..2)
        .map(|i| {
            let p = r.particles[i];
            0.5 * p.m * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
        })
        .sum();
    assert!(
        (ke - 1.0).abs() < 1e-25,
        "a perfectly elastic bounce must conserve kinetic energy, got {}",
        ke
    );
}

/// A coefficient of restitution of 1/2 rescales the relative normal
/// velocity by exactly -1/2 and the kinetic energy by exactly 1/4.
#[test]
fn hardsphere_coefficient_of_restitution_scales_the_rebound() {
    fn eps_half(_r: &reb_simulation, _v: f64) -> f64 {
        0.5
    }
    let mut r = box_sim(8.0, 1, 1, 1);
    r.collision = REB_COLLISION::DIRECT;
    r.collision_resolve = Some(reb_collision_resolve_hardsphere);
    r.coefficient_of_restitution = Some(eps_half);
    r.rand_seed = 3;
    reb_simulation_add(&mut r, mkp(1.0, 1.5, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0));
    reb_simulation_add(&mut r, mkp(1.0, 1.5, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0));

    reb_collision_search(&mut r);

    assert_eq!(r.particles[0].vx.to_bits(), (-0.5f64).to_bits(), "particle 0 vx with eps=1/2");
    assert_eq!(r.particles[1].vx.to_bits(), 0.5f64.to_bits(), "particle 1 vx with eps=1/2");
    let dv_after = r.particles[0].vx - r.particles[1].vx;
    assert_eq!(
        dv_after.to_bits(),
        (-1.0f64).to_bits(),
        "relative normal velocity must be reversed and halved (was +2)"
    );
    let ke: f64 = (0..2)
        .map(|i| {
            let p = r.particles[i];
            0.5 * p.m * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
        })
        .sum();
    assert!(
        (ke - 0.25).abs() < 1e-25,
        "kinetic energy must be scaled by eps^2 = 1/4, got {}",
        ke
    );
}

/// Two overlapping particles at rest are pushed apart at
/// `minimum_collision_velocity` scaled by the overlap depth. With
/// radii 2 and separation 3 the depth factor is exactly 1/2 and the
/// resulting separation speed is exactly `minimum_collision_velocity`.
#[test]
fn hardsphere_minimum_collision_velocity_separates_resting_overlaps() {
    let mut r = box_sim(16.0, 1, 1, 1);
    r.collision = REB_COLLISION::DIRECT;
    r.collision_resolve = Some(reb_collision_resolve_hardsphere);
    r.minimum_collision_velocity = 1.0;
    r.rand_seed = 5;
    reb_simulation_add(&mut r, mkp(1.0, 2.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0));
    reb_simulation_add(&mut r, mkp(1.0, 2.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0));

    reb_collision_search(&mut r);

    // mindv = minr*mcv * (1 - (r-maxr)/minr) = 2*1*(1 - (3-2)/2) = 1
    assert_eq!(r.particles[0].vx.to_bits(), (-0.5f64).to_bits(), "left particle pushed to -x");
    assert_eq!(r.particles[1].vx.to_bits(), 0.5f64.to_bits(), "right particle pushed to +x");
    let separation_speed = r.particles[1].vx - r.particles[0].vx;
    assert_eq!(
        separation_speed.to_bits(),
        1.0f64.to_bits(),
        "separation speed must equal minimum_collision_velocity"
    );
}

/// An oblique, unequal-mass bounce conserves linear momentum and (with
/// eps = 1) kinetic energy.
#[test]
fn hardsphere_oblique_bounce_conserves_momentum_and_energy() {
    let mut r = box_sim(16.0, 1, 1, 1);
    r.collision = REB_COLLISION::DIRECT;
    r.collision_resolve = Some(reb_collision_resolve_hardsphere);
    r.rand_seed = 9;
    reb_simulation_add(&mut r, mkp(1.0, 1.0, -0.7, -0.4, 0.3, 1.3, 0.5, -0.2));
    reb_simulation_add(&mut r, mkp(2.5, 1.0, 0.6, 0.5, -0.4, -0.7, -0.9, 0.6));

    let mom = |r: &reb_simulation| {
        let mut m = (0.0f64, 0.0f64, 0.0f64);
        for i in 0..r.N {
            let p = r.particles[i];
            m.0 += p.m * p.vx;
            m.1 += p.m * p.vy;
            m.2 += p.m * p.vz;
        }
        m
    };
    let ke = |r: &reb_simulation| {
        let mut e = 0.0f64;
        for i in 0..r.N {
            let p = r.particles[i];
            e += 0.5 * p.m * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
        }
        e
    };
    let p0 = mom(&r);
    let e0 = ke(&r);

    reb_collision_search(&mut r);
    assert_eq!(r.N_collisions, 2, "the oblique pair must be detected");

    let p1 = mom(&r);
    let e1 = ke(&r);
    let scale = p0.0.abs().max(p0.1.abs()).max(p0.2.abs()).max(1.0);
    assert!(
        (p1.0 - p0.0).abs() < 1e-14 * scale
            && (p1.1 - p0.1).abs() < 1e-14 * scale
            && (p1.2 - p0.2).abs() < 1e-14 * scale,
        "momentum {:?} -> {:?} across an oblique hard-sphere bounce",
        p0,
        p1
    );
    assert!(
        (e1 - e0).abs() < 1e-13 * e0,
        "kinetic energy {} -> {} across a perfectly elastic bounce",
        e0,
        e1
    );
    // the bounce actually did something
    assert!(
        (r.particles[0].vx - 1.3).abs() > 1e-6,
        "the resolver must have changed the velocities"
    );
}

/// `reb_collision_resolve_merge` conserves mass, momentum and volume, and
/// removes the higher of the two indices. All the numbers here are chosen
/// so that the merge is exact in binary floating point.
#[test]
fn merge_conserves_mass_momentum_and_volume_exactly() {
    let mut r = box_sim(16.0, 1, 1, 1);
    reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.5, 0.0, 0.0, 4.0, 0.0, 0.0));
    reb_simulation_add(&mut r, mkp(3.0, 2.0, -0.25, 0.0, 0.0, 8.0, 0.0, 0.0));

    let c = reb_collision {
        p1: 0,
        p2: 1,
        gb: reb_vec6d::default(),
        ri: 0,
    };
    let outcome = reb_collision_resolve_merge(&mut r, c);
    assert_eq!(
        outcome, REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2,
        "merging (0,1) must ask for the higher index p2 to go"
    );

    let m = r.particles[0];
    assert_eq!(m.m.to_bits(), 4.0f64.to_bits(), "merged mass must be 1+3");
    // v = (4*1 + 8*3)/4 = 7, so the total momentum 28 is preserved exactly
    assert_eq!(m.vx.to_bits(), 7.0f64.to_bits(), "merged vx");
    assert_eq!((m.m * m.vx).to_bits(), 28.0f64.to_bits(), "merged momentum must be 4*1+8*3");
    // x = (0.5*1 + (-0.25)*3)/4 = -0.0625
    assert_eq!(m.x.to_bits(), (-0.0625f64).to_bits(), "merged position is the centre of mass");
    let vol = m.r * m.r * m.r;
    assert!(
        (vol - 9.0).abs() < 1e-14,
        "merged radius must conserve volume 1^3 + 2^3 = 9, got r^3 = {}",
        vol
    );

    // reversing the pair merges into the same particle but flags p1
    let mut r2 = box_sim(16.0, 1, 1, 1);
    reb_simulation_add(&mut r2, mkp(1.0, 1.0, 0.5, 0.0, 0.0, 4.0, 0.0, 0.0));
    reb_simulation_add(&mut r2, mkp(3.0, 2.0, -0.25, 0.0, 0.0, 8.0, 0.0, 0.0));
    let outcome2 = reb_collision_resolve_merge(
        &mut r2,
        reb_collision {
            p1: 1,
            p2: 0,
            gb: reb_vec6d::default(),
            ri: 0,
        },
    );
    assert_eq!(
        outcome2, REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P1,
        "merging (1,0) must ask for the higher index p1 to go"
    );
    assert_eq!(
        r2.particles[0].vx.to_bits(),
        r.particles[0].vx.to_bits(),
        "the merge product must not depend on the order of the pair"
    );
    assert_eq!(
        r2.particles[0].x.to_bits(),
        r.particles[0].x.to_bits(),
        "the merge product position must not depend on the order of the pair"
    );
}

/// With `track_energy_offset` set, merging an isolated pair leaves the
/// total simulation energy unchanged: the offset absorbs the kinetic and
/// potential energy lost to the merger.
#[test]
fn merge_energy_offset_preserves_the_total_energy_of_an_isolated_pair() {
    let mut r = box_sim(16.0, 1, 1, 1);
    r.track_energy_offset = 1;
    reb_simulation_add(&mut r, mkp(1.0, 0.5, -0.6, 0.1, 0.0, 0.4, -0.3, 0.2));
    reb_simulation_add(&mut r, mkp(2.0, 0.5, 0.3, -0.2, 0.1, -0.5, 0.7, -0.1));

    let e_before = reb_simulation_energy(&r);
    let outcome = reb_collision_resolve_merge(
        &mut r,
        reb_collision {
            p1: 0,
            p2: 1,
            gb: reb_vec6d::default(),
            ri: 0,
        },
    );
    assert_eq!(outcome, REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2, "p2 must be flagged");
    assert_eq!(reb_simulation_remove_particle(&mut r, 1), 0, "removing p2 must succeed");
    assert_eq!(r.N, 1, "one merged particle must remain");

    let e_after = reb_simulation_energy(&r);
    assert!(
        (e_after - e_before).abs() < 1e-12 * e_before.abs().max(1.0),
        "energy_offset must keep the total energy: {} -> {} (offset {})",
        e_before,
        e_after,
        r.energy_offset
    );
    assert!(
        r.energy_offset != 0.0,
        "an inelastic merger must have recorded a non-zero energy offset"
    );
}

/// Two independent colliding pairs merged in the same search. Each pair is
/// found twice, so the resolver loop has to renumber and invalidate the
/// pending collisions after every removal. The outcome must not depend on
/// the random order in which the four entries are resolved.
#[test]
fn merge_renumbers_pending_collisions_after_each_removal() {
    for seed in [1u32, 2, 3, 17, 12345, 987654321] {
        let mut r = box_sim(64.0, 1, 1, 1);
        r.collision = REB_COLLISION::DIRECT;
        r.collision_resolve = Some(reb_collision_resolve_merge);
        r.rand_seed = seed;
        reb_simulation_add(&mut r, mkp(1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0));
        reb_simulation_add(&mut r, mkp(1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0));
        reb_simulation_add(&mut r, mkp(2.0, 1.0, 10.0, 0.0, 0.0, 1.0, 0.0, 0.0));
        reb_simulation_add(&mut r, mkp(2.0, 1.0, 11.0, 0.0, 0.0, -1.0, 0.0, 0.0));

        reb_collision_search(&mut r);

        assert_eq!(r.N_collisions, 4, "two pairs, each reported twice (seed {})", seed);
        assert_eq!(r.N, 2, "both pairs must have merged (seed {})", seed);
        assert_eq!(r.particles.len(), 2, "particle array must shrink with N (seed {})", seed);

        assert_eq!(r.particles[0].m.to_bits(), 2.0f64.to_bits(), "merged mass of pair 0-1 (seed {})", seed);
        assert_eq!(r.particles[0].x.to_bits(), 0.5f64.to_bits(), "merged position of pair 0-1 (seed {})", seed);
        assert_eq!(r.particles[0].vx.to_bits(), 0.0f64.to_bits(), "merged velocity of pair 0-1 (seed {})", seed);
        assert_eq!(r.particles[1].m.to_bits(), 4.0f64.to_bits(), "merged mass of pair 2-3 (seed {})", seed);
        assert_eq!(r.particles[1].x.to_bits(), 10.5f64.to_bits(), "merged position of pair 2-3 (seed {})", seed);
        assert_eq!(r.particles[1].vx.to_bits(), 0.0f64.to_bits(), "merged velocity of pair 2-3 (seed {})", seed);

        let px: f64 = (0..r.N).map(|i| r.particles[i].m * r.particles[i].vx).sum();
        assert_eq!(px.to_bits(), 0.0f64.to_bits(), "total momentum must stay zero (seed {})", seed);
    }
}

/// The default resolver halts the integration and leaves the particles
/// completely untouched.
#[test]
fn default_resolver_halts_without_touching_the_particles() {
    let mut r = box_sim(8.0, 1, 1, 1);
    r.collision = REB_COLLISION::DIRECT;
    r.rand_seed = 13;
    reb_simulation_add(&mut r, mkp(1.0, 1.0, -0.4, 0.0, 0.0, 0.3, 0.1, -0.2));
    reb_simulation_add(&mut r, mkp(2.0, 1.0, 0.4, 0.0, 0.0, -0.3, -0.1, 0.2));
    let before: Vec<reb_particle> = r.particles.clone();
    assert_eq!(r.status, 0, "a fresh simulation starts with status 0");

    reb_collision_search(&mut r);

    assert_eq!(r.N_collisions, 2, "the overlapping approaching pair must be found");
    assert_eq!(r.status, REB_STATUS_COLLISION, "reb_collision_resolve_halt must set the status");
    assert_eq!(r.N, 2, "halting must not remove particles");
    for i in 0..2 {
        assert_eq!(r.particles[i].vx.to_bits(), before[i].vx.to_bits(), "particle {} vx untouched", i);
        assert_eq!(r.particles[i].vy.to_bits(), before[i].vy.to_bits(), "particle {} vy untouched", i);
        assert_eq!(r.particles[i].vz.to_bits(), before[i].vz.to_bits(), "particle {} vz untouched", i);
        assert_eq!(r.particles[i].x.to_bits(), before[i].x.to_bits(), "particle {} x untouched", i);
    }
}
