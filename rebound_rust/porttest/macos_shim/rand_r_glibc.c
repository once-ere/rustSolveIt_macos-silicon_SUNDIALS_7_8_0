/* rand_r_glibc.c — the glibc rand_r, for the macOS C reference build.
 *
 * Why this file exists. REBOUND wants its random initial conditions to be
 * identical on every platform, so on Windows — whose C library has no
 * rand_r at all — rebound.c vendors glibc's rand_r directly, guarded by
 * `#ifdef _WIN32` (see src/rebound.c, "Source: codebrowser.dev/glibc/...").
 * On macOS that guard is false, so the C build silently calls Apple's
 * libc rand_r instead — a DIFFERENT generator that produces a different
 * random stream from the same seed. The physics is fine, but the random
 * initial conditions no longer match the ones every recorded reference
 * run used, and a bit-for-bit comparison dies at particle 1.
 *
 * This file is the same vendored glibc algorithm, byte for byte the same
 * arithmetic as the `#ifdef _WIN32` block in rebound.c, compiled as its
 * own object file. Linking it into a harness makes the linker resolve
 * rebound.o's rand_r calls here instead of in Apple's libc — the same
 * stream Windows and glibc-Linux get, without modifying one line of the
 * upstream source. It is the macOS twin of the MSVC shim in
 * reboundx_rust/porttest/msvc_shim/ (which exists for the opposite
 * reason: MSVC lacks a C99 feature that clang has).
 *
 * The Rust port implements this same algorithm (src/tools.rs), which is
 * why the Rust side needs no shim on any platform.
 *
 * Algorithm source, as recorded upstream:
 *   https://codebrowser.dev/glibc/glibc/stdlib/rand_r.c.html
 * GPL-3.0-or-later, same as REBOUND and this repository.
 */

int rand_r(unsigned int *seed) {
    unsigned int next = *seed;
    int result;

    next *= 1103515245;
    next += 12345;
    result = (unsigned int) (next / 65536) % 2048;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;

    *seed = next;

    return result;
}
