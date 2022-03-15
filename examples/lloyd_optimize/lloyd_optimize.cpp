#define CGAL_MESH_2_OPTIMIZER_VERBOSE
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
	std::ifstream in(mesh);
	in >> cdt;

	// CGAL::draw(cdt);

	CGAL::lloyd_optimize_mesh_2(cdt,
		CGAL::parameters::time_limit = timeLimit,
		CGAL::parameters::max_iteration_number = iterations,
		CGAL::parameters::convergence = convergenceRatio,
		CGAL::parameters::freeze_bound = freezeBound);

	// CGAL::draw(cdt);

	std::ofstream out(result);
	out << cdt;

	return EXIT_SUCCESS;
}
