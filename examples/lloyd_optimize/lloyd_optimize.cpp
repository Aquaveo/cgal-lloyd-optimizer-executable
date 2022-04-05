#define CGAL_MESH_2_OPTIMIZER_VERBOSE

constexpr int mode_none = 0;
constexpr int mode_lock = 1;
constexpr int mode_delete = 2;

int cgal_constraint_mode = mode_lock;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <CGAL/lloyd_optimize_mesh_2.h>

#include <CGAL/boost/graph/IO/VTK.h>

#include <CGAL/draw_triangulation_2.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;


void loadMesh(std::ifstream& a_in, CDT& a_cdt)
{
	std::vector<CDT::Point_2> points;
	size_t numPoints;
	a_in >> numPoints;
	for (size_t i = 0; i < numPoints; i++)
	{
		double x, y;
		a_in >> x >> y;
		CDT::Point_2 point(x, y);
		a_cdt.insert(point);
		points.push_back(point);
	}

	size_t numBoundaries;
	a_in >> numBoundaries;
	for (size_t i = 0; i < numBoundaries; i++)
	{
		size_t a, b;
		a_in >> a >> b;
		a_cdt.insert_constraint(points[a], points[b]);
	}
}

void saveMesh(const std::string& a_file, CDT& a_cdt)
{
	std::ofstream out(a_file);

	out << a_cdt.number_of_vertices() << '\n';

	for (auto& it = a_cdt.finite_vertices_begin(); it != a_cdt.finite_vertices_end(); ++it)
	{
		out << it->point().x() << ' ' << it->point().y() << '\n';
	}
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " filename\n";
		return 1;
	}

	std::string inFile(argv[1]), outFile(argv[1]);
	inFile += ".in";
	outFile += ".out";

	int iterations = 0;
	double convergenceRatio = 0;
	double freezeBound = 0;
	int timeLimit = 0;
	CDT cdt;

	std::ifstream in(inFile);

	in >> iterations >> timeLimit >> convergenceRatio >> freezeBound;

	loadMesh(in, cdt);

	CGAL::draw(cdt);

	CGAL::lloyd_optimize_mesh_2(cdt,
		CGAL::parameters::time_limit = timeLimit,
		CGAL::parameters::max_iteration_number = iterations,
		CGAL::parameters::convergence = convergenceRatio,
		CGAL::parameters::freeze_bound = freezeBound,
		CGAL::parameters::mark = true);

	CGAL::draw(cdt);

	saveMesh(outFile, cdt);

	return EXIT_SUCCESS;
}
