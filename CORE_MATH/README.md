# The CORE-MATH project

CORE-MATH Mission: provide on-the-shelf open-source mathematical
functions with correct rounding that can be integrated into current
mathematical libraries (GNU libc, Intel Math Library, AMD Libm,
Newlib, OpenLibm, Musl, Apple Libm, llvm-libc, CUDA libm, ROCm).

Homepage: https://core-math.gitlabpages.inria.fr/

In this project, we have as log in Binary64 of CORE-MATH project. 

## Quick guide

### Worst case checks

These checks are available for double-precision functions.

    ./check.sh [rounding_modes] log

where:
- `[rounding_modes]` can be a selection of `--rndn` (round to
  nearest), `--rndz` (toward zero), `--rndu` (upwards), `--rndd`
  (downwards). The default is all four.

This command is sensitive to the following environment variables:
- `CC`
- `CFLAGS`
- OpenMP variables such as `OMP_NUM_THREADS`

Note: on Debian, you need the libomp-dev package to use clang.

### Performance measurements

Performance measurement scripts rely on the Linux perf framework. You
might need to run the command:

    echo -1 > /proc/sys/kernel/perf_event_paranoid

as root before running the following commands.

To evaluate the reciprocal throughput of some function, run:

    ./perf.sh log

It outputs two numbers: the performance of core-math, and the one the
libc, in cycles/call. You can also evaluate the performance of some
other libm (given by a static .a archive, typically llvmlibc), by
setting the `LIBM` environment variable to the absolute path of the .a
file. In the case, `./perf.sh` will output three numbers: the
performances of core-math, standard libm, given libm.
