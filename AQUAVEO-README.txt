This optimizer uses a slightly modified version of CGAL and a heavily modified example.

Modifications made by Aquaveo were made to Mesh_global_optimizer_2.h. Insertions to std::cerr were replaced with insertions to std::cout, and the string written at the end of each iteration was altered so it wouldn't erase the previous output. The optimizer was also modified to look for a cancelOptimization variable so it can be aborted from another thread.

The file lloyd_optimize.cpp is loosely based on the draw_triangulation_2 example.

To build:
1. In terminal, run:
- mkdir build
- cd build
- conan install .. -pr ../profile.txt
2. Run CMake, checking the WITH_examples option
3. Open build/CGAL.sln in Visual Studio
4. Build the lloyd_optimize project
5. Go to https://github.com/CGAL/cgal/releases and find the release that matches this version of CGAL. Currently, we're using CGAL-5.4
6. In the assets section for the release, download the GMP and MPFR libraries
7. Copy <zip-file>/auxiliary/gmp/lib/libgmp-10.dll to <readme-folder>/build/examples/lloyd_optimize/<Debug or Release>/

To package:
Package the executable and DLL from step 7 in building.