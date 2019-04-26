#include "main.h"

#include "../common/typedefs.h"
#include "../common/constants.h"
#include "../embedding/io.h"
#include "../embedding/debug_geometry.h"

#include "debug_utils.h"
#include "impersonator.h"
#include "intersections_deprecated.h"
#include "intersections_new.h"
#include "mesh_merger.h"
#include "test_overlay.h"

#include "overlay.h"

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

using namespace std;

#define TEST false

// Example arguments:
// -t
// -v "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Banana\Banana.obj" ""D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Banana\_Banana[output_tweaked].obj"" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Deer\fixed_deer.obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Deer\_fixed_deer[output_tweaked].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Deer\dearPLUSbanana.obj"
// -v "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (1).obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (1).obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (2).obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (2).obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (1,2).obj"

int main_overlay(int argc, char **argv) {
#if TEST
	run_all_tests();
#else	
	if (strcmp(argv[1], "-t") == 0) {
		run_all_tests();
		return 0;
	}


	bool verbose = false;
	int n = 0;
	if (strcmp(argv[1], "-v") == 0) { verbose = true; n++;  }
	std::string input_path_1 = argv[1 + n];
	std::string input_path_1_embedding = argv[2 + n];
	std::string input_path_2 = argv[3 + n];
	std::string input_path_2_embedding = argv[4 + n];
	std::string output_path  = argv[5+n];
	std::string output_path_merged_embedding = output_path.substr(0, output_path.find_last_of(".")) + "[merged-embedding]" + output_path.substr(output_path.find_last_of("."));
	std::string output_path_impersonation_1 = output_path.substr(0, output_path.find_last_of(".")) + "[impers_1]" + output_path.substr(output_path.find_last_of("."));
	std::string output_path_impersonation_2 = output_path.substr(0, output_path.find_last_of(".")) + "[impers_2]" + output_path.substr(output_path.find_last_of("."));
	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;

	std::vector<Point> points_1, points_2;
	std::vector<std::vector<size_t>> faces_1, faces_2;
	Polyhedron mesh_1, mesh_2;

	std::vector<Point> original_mesh_points_1, original_mesh_points_2;             // passed to 'impersonate()'
	std::vector<std::vector<size_t>> original_mesh_faces_1, original_mesh_faces_2; // unused
	Polyhedron original_mesh_1, original_mesh_2;                                   // unused

	bool succ = input_mesh(points_1, faces_1, mesh_1, vout, input_path_1_embedding);
	succ = succ && input_mesh(points_2, faces_2, mesh_2, vout, input_path_2_embedding);
	succ = succ && input_mesh(original_mesh_points_1, original_mesh_faces_1, original_mesh_1, vout, input_path_1);
	succ = succ && input_mesh(original_mesh_points_2, original_mesh_faces_2, original_mesh_2, vout, input_path_2);
	if (!succ) {
		cerr << "Error reading input meshes." << endl;
		ask_before_termination(DISABLE_USER_INTERACTION);
		return -1;
	}

	std::tuple<Polyhedron, Polyhedron, Polyhedron> overlay_results = overlay(mesh_1, mesh_2, original_mesh_1, original_mesh_2, vout);
	Polyhedron mesh_impersonation_1 = std::get<0>(overlay_results);
	Polyhedron mesh_impersonation_2 = std::get<1>(overlay_results);
	Polyhedron merged_embedding = std::get<2>(overlay_results);

	// Output the merged embedding.
	output_mesh(output_path_merged_embedding, merged_embedding, vout, true);

	// Output impersonated mesh.
	output_mesh(output_path_impersonation_1, mesh_impersonation_1, vout, true);

	// Output impersonated mesh.
	output_mesh(output_path_impersonation_2, mesh_impersonation_2, vout, true);

	ask_before_termination(DISABLE_USER_INTERACTION);
	return 0;
#endif
}