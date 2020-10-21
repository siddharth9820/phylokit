# phylokit
C++ library for high performance phylogenetics. Used in (e.g.) ASTRID, FastRFS, SIESTA.

To build, you will need bazel (https://www.bazel.build/):

    bazel build phylokit:all
    bazel test test:all

It's easiest to use phylokit from within a bazel project;
copy the contents of phylokit's WORKSPACE into your WORKSPACE file.

If you want to use phylokit from a cmake project, you could try https://github.com/google/bazel-to-cmake.