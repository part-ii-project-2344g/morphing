#include "main.h"
#include "../common/typedefs.h"
#include "../common/constants.h"
#include "../embedding/io.h"
#include "parse_features.h"
#include "feature_alignment.h"

using namespace std;

// Example arguments:
// -v "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\tetrahedron_features.feat" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (1).obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (2).obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (1)[aligned].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\Simple Tetrahedron (2)[aligned].obj"

#define TEST false

int main_feature_alignment(int argc, char **argv) {
#if TEST

#else
	bool verbose = false;
	int n = 0;
	if (strcmp(argv[1], "-v") == 0) { verbose = true; n++; }
	std::string input_path_features = argv[1 + n];
	std::string input_path_mesh_1 = argv[2 + n];
	std::string input_path_mesh_2 = argv[3 + n];
	std::string output_path_mesh_1 = argv[4 + n];
	std::string output_path_mesh_2 = argv[5 + n];
	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;

	std::vector<Point> points_1, points_2;
	std::vector<std::vector<size_t>> faces_1, faces_2;
	Polyhedron mesh_1, mesh_2;

	std::vector<Point> original_mesh_points_1, original_mesh_points_2;             // passed to 'impersonate()'
	std::vector<std::vector<size_t>> original_mesh_faces_1, original_mesh_faces_2; // unused
	Polyhedron original_mesh_1, original_mesh_2;                                   // unused

	bool succ = input_mesh(points_1, faces_1, mesh_1, vout, input_path_mesh_1);
	succ = succ && input_mesh(points_2, faces_2, mesh_2, vout, input_path_mesh_2);
	if (!succ) {
		cerr << "Error reading input meshes." << endl;
		ask_before_termination(DISABLE_USER_INTERACTION);
		return -1;
	}

	// If feature path is empty, this means we are doing no alignment.
	if (input_path_features == "") {
		std::cout << "Feature alignment doing a NOOP." << std::endl;
		// Output unchanged meshes.
		output_mesh(output_path_mesh_1, mesh_1, vout, true);
		output_mesh(output_path_mesh_2, mesh_2, vout, true);
		return 0;
	}

	// Parse features.
	std::vector<std::pair<size_t, size_t>> features = read_features(input_path_features);

	// Pass everything to the function.
	auto aligned_meshes = align_features(mesh_1, mesh_2, features, output_path_mesh_1, vout);
	Polyhedron aligned_mesh_1 = aligned_meshes.first;
	Polyhedron aligned_mesh_2 = aligned_meshes.second;

	// Output the resulting meshes.
	output_mesh(output_path_mesh_1, aligned_mesh_1, vout, true);
	output_mesh(output_path_mesh_2, aligned_mesh_2, vout, true);

	ask_before_termination(DISABLE_USER_INTERACTION);
	return 0;
#endif
}
