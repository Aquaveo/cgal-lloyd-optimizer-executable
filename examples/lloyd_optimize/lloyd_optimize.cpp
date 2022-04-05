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


template <typename T>
T max()
{
	return std::numeric_limits<T>::max();
}

template <typename T>
bool isMax(T value)
{
	return value == std::numeric_limits<T>::max();
}

/// <summary>
/// Try parsing an argument. If the argument was present but invalid, throws std::invalid_argument.
/// </summary>
/// <typeparam name="T">Type of argument.</typeparam>
/// <param name="a_args">Vector of arguments to try and parse from.</param>
/// <param name="a_arg">Name of argument.</param>
/// <param name="a_result">The parsed value.</param>
/// <param name="a_description">Description of the argument for error messages.</param>
/// <returns>Whether the argument was parsed.</returns>
template <typename T>
bool tryParseArgument(
	std::vector<std::string>& a_args,
	const std::string& a_arg,
	T& a_result,
	const std::string& a_description)
{
	if (a_args.front() == a_arg)
	{
		if (a_args.size() > 1)
		{
			std::stringstream stream(a_args[1]);
			stream >> a_result;

			a_args.erase(a_args.begin());
			a_args.erase(a_args.begin());

			return true;
		}
		else
		{
			std::cout << a_arg << " argument requires " << a_description << "\n";
			throw std::invalid_argument("Invalid argument");
		}
	}
	return false;
}

/// <summary>
/// Print usage message.
/// </summary>
void printUsage()
{
	std::cout <<
		"Usage:\n"
		"--mesh <path>          : Path to file containing mesh. Required.\n"
		"--out <path>           : Path to output file. Required.\n"
		"--iterations <integer> : Number of iterations to stop after. Defaults to no limit.\n"
		"--timeout <integer>    : Number of seconds to stop after. Defaults to never.\n"
		"--ratio <double>       : Fraction of shortest edge incident to vertex that the vertex must move by.\n"
		"                         If any vertex moves by less than this fraction of its shortest edge, the\n"
		"                         optimization is considered to have converged. Defaults to no minimum.\n"
		"                         Requires 0 <= ratio <=1.\n"
		"--freeze <double>      : Fraction of shortest edge incident to vertex that the vertex must move by.\n"
		"                         If any vertex would move by less than this fraction of its shortest edge,\n"
		"                         then the move will be ignored to save optimization time. Defaults to no\n"
		"                         minimum movement. Requires 0 <= freeze <= 1.\n";
}

/// <summary>
/// Parse arguments out of the command line.
/// </summary>
/// <param name="argc">Number of arguments.</param>
/// <param name="argv">Argument values.</param>
/// <param name="a_mesh">Path to mesh.</param>
/// <param name="a_out">Path to output file.</param>
/// <param name="a_iterations">Number of iterations.</param>
/// <param name="a_convergenceRatio">Minimum fraction movement for non-convergence.</param>
/// <param name="a_freezeBound">Minimum fraction movement to move.</param>
/// <param name="a_timeLimit">Number of seconds to optimize for.</param>
/// <returns>Whether arguments were parsed successfully.</returns>
bool parseArguments(
	int argc,
	char* argv[],
	std::string& a_mesh,
	std::string& a_out,
	int& a_iterations,
	double& a_convergenceRatio,
	double& a_freezeBound,
	int& a_timeLimit)
{
	std::vector<std::string> args(argv, argv+argc);
	args.erase(args.begin()); // Discard program name

	bool parsed = true;

	try
	{
		while (!args.empty() && parsed)
		{
			parsed = false;
			parsed = parsed || tryParseArgument(args, "--mesh", a_mesh, "path to mesh");
			parsed = parsed || tryParseArgument(args, "--out", a_out, "path to output");
			parsed = parsed || tryParseArgument(args, "--iterations", a_iterations, "number of iterations");
			parsed = parsed || tryParseArgument(args, "--timeout", a_timeLimit, "number of seconds to timeout");
			parsed = parsed || tryParseArgument(args, "--ratio", a_convergenceRatio, "convergence ratio");
			parsed = parsed || tryParseArgument(args, "--freeze", a_freezeBound, "freeze bound");
		}

		if (!parsed)
		{
			std::cout << "Unknown argument: " << args[0] << "\n";
			return false;
		}

		if (a_mesh == "")
		{
			std::cout << "No mesh specified.\n";
			return false;
		}
	}
	catch (std::invalid_argument&)
	{
		return false;
	}
	
	return true;
}


bool loadMesh(const std::string& a_file, CDT& a_cdt, std::vector<CDT::Point_2>& a_seeds)
{
	std::ifstream in(a_file);

	if (!in)
	{
		std::cerr << "Unable to open input file.\n";
		return false;
	}

	size_t numPoints = 0;
	in >> numPoints;
	std::vector<Vertex_handle> points;
	points.reserve(numPoints);
	for (size_t i = 0; i < numPoints; i++)
	{
		double x = max<double>(),
			   y = max<double>();
		in >> x >> y;

		if (isMax(x) || isMax(y))
		{
			std::cerr << "Point " << i << " was bad.\n";
			return false;
		}

		std::cout << "Inserting point (" << x << ", " << y << ")\n";
		points.push_back(a_cdt.insert(Point(x, y)));
	}

	size_t numBoundaries = 0;
	in >> numBoundaries;
	for (size_t i = 0; i < numBoundaries; i++)
	{
		size_t a = max<size_t>(),
			   b = max<size_t>();
		in >> a >> b;

		if (isMax(a) || isMax(b))
		{
			std::cerr << "Boundary " << i << " was bad.\n";
			return false;
		}

		std::cout << "Inserting boundary (" << points[a]->point().x() << ", " << points[a]->point().y() << "), (" << points[b]->point().x() << ", " << points[b]->point().y() << ")\n";
		a_cdt.insert_constraint(points[a], points[b]);
	}

	size_t numSeeds = 0;
	in >> numSeeds;
	a_seeds.clear();
	a_seeds.reserve(numSeeds);
	for (size_t i = 0; i < numSeeds; i++)
	{
		size_t s = max<size_t>();
		in >> s;

		if (isMax(s))
		{
			std::cerr << "Seed " << i << " was bad.\n";
			return false;
		}

		std::cout << "Seed point: (" << points[s]->point().x() << ", " << points[s]->point().y() << ")\n";
		a_seeds.push_back(points[s]->point());
	}

	return true;
}

bool saveMesh(const std::string& a_file, CDT& a_cdt)
{
	std::ofstream out(a_file);

	if (!out)
	{
		std::cerr << "Unable to open output file.\n";
		return false;
	}

	std::map<Vertex_handle, size_t> vMap;

	size_t numPoints = a_cdt.number_of_vertices(),
		   numCells = a_cdt.number_of_faces(),
		   numDimensions = a_cdt.dimension();

	out << numPoints << ' ' << numCells << ' ' << numDimensions << '\n';

	for (auto& it = a_cdt.finite_vertices_begin(); it != a_cdt.finite_vertices_end(); ++it)
	{
		vMap.emplace(it, vMap.size());
		out << it->point().x() << ' ' << it->point().y() << '\n';
	}
	out << '\n';

	for (auto& it = a_cdt.finite_faces_begin(); it != a_cdt.finite_faces_end(); ++it)
	{
		for (int i = 0; i < 3; i++)
		{
			auto v = it->vertex(i);
			size_t index = vMap[v];
			out << index << (i == 2 ? '\n' : ' ');
		}
	}

	return true;
}

int main(int argc, char* argv[])
{
	std::string mesh;
	std::string result;
	int iterations = 0;
	double convergenceRatio = 0;
	double freezeBound = 0;
	int timeLimit = 0;

	bool parsed = parseArguments(argc, argv, mesh, result, iterations, convergenceRatio, freezeBound, timeLimit);

	if (!parsed)
	{
		return 1;
	}

	CDT cdt;
	std::vector<CDT::Point_2> seeds;
	loadMesh(mesh, cdt, seeds);

	CGAL::draw(cdt);

	CGAL::lloyd_optimize_mesh_2(cdt,
		CGAL::parameters::time_limit = timeLimit,
		CGAL::parameters::max_iteration_number = iterations,
		CGAL::parameters::convergence = convergenceRatio,
		CGAL::parameters::freeze_bound = freezeBound,
		CGAL::parameters::mark = true);
	
	CGAL::draw(cdt);

	saveMesh(result, cdt);

	return EXIT_SUCCESS;
}
