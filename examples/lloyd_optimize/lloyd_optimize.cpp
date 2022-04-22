#define CGAL_MESH_2_OPTIMIZER_VERBOSE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>

#include <CGAL/lloyd_optimize_mesh_2.h>

#include <fstream>
#include <iostream>
#include <string>
#include <thread>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

bool stopOptimizing = false;

void watchForCancel()
{
	std::cin.get();
	stopOptimizing = true;
}

void loadMesh(std::ifstream& a_in, CDT& a_cdt)
{
	std::vector<CDT::Point_2> points;
	size_t numPoints = 0;
	a_in >> numPoints;
	for (size_t i = 0; i < numPoints; i++)
	{
		double x, y;
		a_in >> x >> y;
		CDT::Point_2 point(x, y);
		a_cdt.insert(point);
		points.push_back(point);
	}

	size_t numBoundaries = 0;
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
	std::map<CDT::Point_2, size_t> map;
	size_t id = 0;
	std::ofstream out(a_file);

	for (auto& vertex : a_cdt.finite_vertex_handles())
	{
		map[vertex->point()] = id++;
		out << vertex->point().x() << ' ' << vertex->point().y() << '\n';
	}

	out << "-\n";

	for (auto& face : a_cdt.finite_face_handles())
	{
		size_t a = map[face->vertex(0)->point()];
		size_t b = map[face->vertex(1)->point()];
		size_t c = map[face->vertex(2)->point()];
		out << a << ' ' << b << ' ' << c << '\n';
	}
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " filename\n";
		return 1;
	}

	std::string directory = argv[1];
	if (directory[directory.length() - 1] != '/' && directory[directory.length() - 1] != '\\')
	{
		directory += '\\';
	}

	std::string inFile = directory + "to-cgal", outFile = directory + "to-xms";

	int iterations = 0;
	double convergenceRatio = 0;
	double freezeBound = 0;
	int timeLimit = 0;
	CDT cdt;

	std::ifstream in(inFile);
	if (!in.is_open())
	{
		std::cerr << "Unable to open file: " << inFile << '\n';
		return 1;
	}

	in >> iterations >> timeLimit >> convergenceRatio >> freezeBound;
	
	loadMesh(in, cdt);

	std::thread watcher(watchForCancel);

	CGAL::lloyd_optimize_mesh_2(cdt,
		CGAL::parameters::time_limit = timeLimit,
		CGAL::parameters::max_iteration_number = iterations,
		CGAL::parameters::convergence = convergenceRatio,
		CGAL::parameters::freeze_bound = freezeBound,
		CGAL::parameters::mark = true);

	saveMesh(outFile, cdt);

	watcher.detach();

	return EXIT_SUCCESS;
}
