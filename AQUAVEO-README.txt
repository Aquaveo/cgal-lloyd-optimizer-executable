This optimizer uses a slightly modified version of CGAL and a heavily modified example.

Modifications made by Aquaveo were made to Mesh_global_optimizer_2.h. Insertions to std::cerr were replaced with insertions to std::cout, and the string written at the end of each iteration was altered so it wouldn't erase the previous output. The optimizer was also modified to look for a cancelOptimization variable so it can be aborted from another thread.

The file lloyd_optimize.cpp is loosely based on the draw_triangulation_2 example.