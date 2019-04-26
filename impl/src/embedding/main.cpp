// Produces a spherical embedding of a mesh.
#include "main.h"

#include <CGAL/Surface_mesh_deformation.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/property_map.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/Inverse_index.h>
#include <CGAL/enum.h>  // ORIGIN
#include <unordered_set>

#include "../common/typedefs.h"
#include "../common/constants.h"

#include "approx_bounding_sphere.h"
#include "debug_geometry.h"
#include "geometry_utils.h"
#include "io.h"
#include "random_tetrahedron.h"
#include "relaxation.h"

using namespace std;

#define EPSILON_MULTIPLIER_STEP 0.5
#define ORIENTATION_CHECK_FREQUENCY 2000
#define EXPORT_FREQUENCY 2000
#define SWAP_FREQUENCY 100000


// Example arguments:
// -v -s 3.14825 -0.004761 2.77544 0.375797 0.5 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Deer\fixed_deer.obj"
// -v -s 1.091 0.0 0.0 0.0 0.5 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Banana\Banana.obj"

int main_embedding(int argc, char **argv) {
	// Setup.
	std::string input_path, output_path;
	bool verbose = false;
	double RELAXATION_EPSILON = 0.0;
	Sphere bounding_sphere_to_use{ 0.0, Point{} };
	parse_args(argc, argv, input_path, output_path, verbose, RELAXATION_EPSILON, bounding_sphere_to_use);
	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;
	// End of setup.


	// Read in the mesh into points, faces, and a Polyhedron mesh.
	std::vector<Point> points{};
	std::vector< std::vector<std::size_t> > faces{};
	Polyhedron mesh;
	if (!input_mesh(points, faces, mesh, vout, std::string{ input_path })) {
		cerr << "Error reading mesh.\n";
		ask_before_termination(DISABLE_USER_INTERACTION);
		return -1;
	}

	// debug
	size_t debug_points_size = points.size();

	// debug
	debug_print_faces_and_points(points, faces, vout);

	// debug
	std::vector<Point> points_copy_for_debug{};
	for (Point p : points)
		points_copy_for_debug.push_back(Point{ p });

	// Make sure we're working on a triangle mesh. Modifies points, faces, and mesh.
	// Does nothing if the mesh is already a triangle mesh.
	triangulate_mesh_if_necessary(points, faces, mesh, vout, std::string{ input_path });

	// debug
	for (int i = 0; i < points.size(); i++) {
		CGAL_assertion(points[i] == points_copy_for_debug[i]);
	}

	// Compute a bounding sphere.
	Sphere sph;
	if (bounding_sphere_to_use.radius != 0) {
		sph = bounding_sphere_to_use;
		vout << "Using custom bounding sphere..." << endl;
	}
	else {
		sph = compute_approx_bounding_sphere(points);
		vout << "Automatically computing bounding sphere..." << endl;
	}
	vout << "bounding_sphere.centre = " << sph.centre << endl;
	vout << "bounding_sphere.radius = " << sph.radius << endl;

	for (Point p : points) {
		CGAL_assertion((p - sph.centre).squared_length() <= sph.radius*sph.radius);
	}

	// Normalize the model so that the bounding sphere is a unit sphere centered at zero.
	normalize_bounding_sphere(points, sph);

	// Normalize all the points so that they lie on the sphere. Perhaps at this point the embedding will be invalid (self-intersecting).
	// Don't...
	normalize_all_points(points);

	// debug
	vout << "\nAfter normalize_bounding_sphere(..):" << endl;
	debug_print_faces_and_points(points, faces, vout);
	debug_print_face_types(points, faces, vout);
	CGAL_assertion(points.size() == debug_points_size);

	unsigned int seed = 0;
	CGAL::cpp11::tuple<Point, Point, Point, Point> tetrahedron_vertices;
	while (true) {
		// Initialize the relaxation process.
		tetrahedron_vertices = choose_random_tetrahedron(points, seed++);
		auto fixed_vertices = choose_fixed_vertices(points, tetrahedron_vertices);

		// Work on a copy in case we need to retry with a new seed.
		std::vector<Point> points_copy{};
		for (Point p : points)
			points_copy.push_back(Point{ p });
		Polyhedron mesh_copy;
		make_mesh_from_points_and_faces(points_copy, faces, mesh_copy);
		CGAL_assertion(points_copy.size() == debug_points_size);

		// Relax until the change between steps is small enough.
		size_t iterations = 0;
		while (relaxation_step(points_copy, faces, mesh_copy, fixed_vertices, vout, std::string{ input_path }, ((iterations++) % (EXPORT_FREQUENCY)) == 0) > RELAXATION_EPSILON) {
			;
		}

		// If the embedding is collapsed, retry. Otherwise, proceed to step 7.
		if (!is_embedding_collapsed(points_copy, faces, mesh_copy, tetrahedron_vertices, vout)) {
			points = points_copy;
			mesh = mesh_copy;
			CGAL_assertion(points.size() == debug_points_size);
			break;
		}
		else {
			vout << "The embedding was collapsed. Retrying with a different seed tetrahedron." << endl;
		}
	}
	// debug
	vout << "\nBefore step 7:" << endl;
	debug_print_faces_and_points(points, faces, vout);
	CGAL_assertion(points.size() == debug_points_size);
	for (Point p : points) { CGAL_assertion(abs((p - CGAL::ORIGIN).squared_length() - 1.0) < 0.001); }

	std::string input_name{ input_path };
	//output_mesh(input_name.substr(0, input_name.length() - 4) + std::string{ "[before_step_7].obj" }, mesh, vout);

	// Step 7.	
	auto diametric_tetrahedron_vertices = diametric(tetrahedron_vertices);
	auto diametric_fixed_vertices = choose_fixed_vertices(points, diametric_tetrahedron_vertices);
	bool fixed_vertices_are_original = false;
	size_t swap_counter = 0;

	double epsilon_multiplier = 1.0;
	size_t iterations = 0;
	size_t next_check = ORIENTATION_CHECK_FREQUENCY;
	while (!is_orientation_consistent_on_sphere(points, faces)) {
		// Relax until the change between steps is small enough.
		while (relaxation_step(points, faces, mesh, diametric_fixed_vertices, vout, std::string{ input_path }, ((iterations++) % (EXPORT_FREQUENCY)) == 0) > RELAXATION_EPSILON * epsilon_multiplier && iterations < next_check) {
			;
		}
		// Change fixed vertices.
		//if (swap_counter >= SWAP_FREQUENCY) {
		//	swap_counter = 0;
		//	fixed_vertices_are_original = !fixed_vertices_are_original;
		//	diametric_fixed_vertices = choose_fixed_vertices(points, fixed_vertices_are_original ? tetrahedron_vertices : diametric_tetrahedron_vertices);
		//	cerr << "Swapped fixed vertices." << endl;
		//}
		// Keep going until a smaller relaxation epsilon.
		if (epsilon_multiplier == 0.0) {
			cerr << "Error: epsilon_multipler dropped to 0.0 and we still haven't succeeded" << endl;
			return -1;
		}
		if (iterations >= next_check) {
			next_check += ORIENTATION_CHECK_FREQUENCY;
			vout << "Orientation inconsistent after " << iterations << " relaxation steps." << endl;
			swap_counter += ORIENTATION_CHECK_FREQUENCY;
		}
		else {
			epsilon_multiplier *= EPSILON_MULTIPLIER_STEP;
			vout << "Decreasing the epsilon to " << RELAXATION_EPSILON * epsilon_multiplier << " as the orientation is still inconsistent." << endl;
		}
	}
	// the embedding was successfully created
	vout << "The embedding was successfully created!" << endl;
	
	// debug
	for (Point p : points) {
		CGAL_assertion(abs((p - CGAL::ORIGIN).squared_length() - 1.0) < 0.001);
	}

	// Output
	output_mesh(output_path, mesh, vout);
	ask_before_termination(DISABLE_USER_INTERACTION);
	return 0;
}
// ------------------------------------------------------|================================================================================
// end of the main function                              |================================================================================
// ------------------------------------------------------|================================================================================
// ------------------------------------------------------|================================================================================
//_______________________________________________________|================================================================================
