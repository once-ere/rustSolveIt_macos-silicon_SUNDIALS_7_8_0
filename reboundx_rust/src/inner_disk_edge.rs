//! inner_disk_edge.rs — translation of REBOUNDx inner_disk_edge.c
//! Inner disk edge implemention at a chosen location, while planets are undergoing migration.
//!
//! Part of reboundx_rs, a GPL-3.0-or-later translation of REBOUNDx 5.1.0
//! (c) Dan Tamayo, Hanno Rein and contributors. See crate root.
//!
//! @author Kaltrina Kajtazi <1kaltrinakajtazi@gmail.com>
//!
//! # $Orbit Modifications$
//!
//! ======================= ================================================================================================================
//! Authors                 Kajtazi, Kaltrina and D. Petit, C. Antoine
//! Implementation Paper    `Kajtazi et al 2022 <https://ui.adsabs.harvard.edu/abs/2022arXiv221106181K/abstract>`_.
//! Based on                `Pichierri et al 2018 <https://ui.adsabs.harvard.edu/abs/2018CeMDA.130...54P/abstract>`_.
//! C example               :ref:`c_example_inner_disk_edge`
//! Python example          `InnerDiskEdge.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/InnerDiskEdge.ipynb>`_.
//! ======================= ================================================================================================================
//!
//! This applies an inner disk edge that functions as a planet trap. Within its width the planet's
//! migration is reversed by an opposite and roughly equal magnitude torque. Thus, stopping further
//! migration and trapping the planet within the width of the trap.
//! The functions here provide a way to modify the tau_a timescale in modify_orbits_forces,
//! modify_orbit_direct, and type_I_migration.
//! Note that the present prescription is very useful for simple simulations when an inner trap is
//! needed during the migration but it shouldn't be considered as a realistic model of the inner
//! edge of a disk.
//!
//! **Effect Parameters**
//!
//! ============================ =========== ===================================================================================
//! Field (C type)               Required    Description
//! ============================ =========== ===================================================================================
//! ide_position (double)        Yes         The position of the inner disk edge in code units
//! ide_width (double)           Yes         The disk edge width (planet will stop within ide_width of ide_position)
//! ============================ =========== ===================================================================================

use std::f64::consts::PI;

/// A planet trap that is active only at the inner disk edge, to reverse the planetary migration
/// and prevent migration onto the star.
///
/// C: `const double rebx_calculate_planet_trap(const double r, const double dedge, const double hedge)`
pub fn rebx_calculate_planet_trap(r: f64, dedge: f64, hedge: f64) -> f64 {
    let tau_a_red: f64;

    if r > dedge * (1.0 + hedge) {
        tau_a_red = 1.0;
    } else if dedge * (1.0 - hedge) < r {
        tau_a_red =
            5.5 * (((dedge * (1.0 + hedge) - r) * 2.0 * PI) / (4.0 * hedge * dedge)).cos() - 4.5;
    } else {
        tau_a_red = -10.0;
    }

    tau_a_red
}
