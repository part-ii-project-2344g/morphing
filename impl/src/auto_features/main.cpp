#include "main.h"
#include "common/constants.h"
#include "simple_features.h"
#include "../embedding/io.h"

using namespace std;

// Example arguments:
// -v poisson mesh1 mesh2 mesh1_lowres mesh2_lowres outpath 0.1 3 0.3 1.0 0.3

int main_simple_features(int argc, char **argv) {
	bool verbose = false;
	int n = 0;
	if (strcmp(argv[1], "-v") == 0) { verbose = true; n++; }
	std::string mode = argv[1 + n];
	std::string input_path_mesh_1 = argv[2 + n];
	std::string input_path_mesh_2 = argv[3 + n];
	std::string input_path_low_res_1 = argv[4 + n];
	std::string input_path_low_res_2 = argv[5 + n];
	std::string output_path_features = argv[6 + n];
	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;


	std::vector<Point> points_1, points_2, points_lr_1, points_lr_2;
	std::vector<std::vector<size_t>> faces_1, faces_2, faces_lr_1, faces_lr_2;
	Polyhedron mesh_1, mesh_2, low_res_1, low_res_2;

	bool succ = input_mesh(points_1, faces_1, mesh_1, vout, input_path_mesh_1);
	succ = succ && input_mesh(points_2, faces_2, mesh_2, vout, input_path_mesh_2);
	succ = succ && input_mesh(points_lr_1, faces_lr_1, low_res_1, vout, input_path_low_res_1);
	succ = succ && input_mesh(points_lr_2, faces_lr_2, low_res_2, vout, input_path_low_res_2);
	if (!succ) {
		cerr << "Error reading input meshes." << endl;
		ask_before_termination(DISABLE_USER_INTERACTION);
		return -1;
	}

	if (mode == "simple")
		simple_features_to_file(mesh_1, mesh_2, output_path_features);
	else if (mode == "spatial")
		spatial_simple_features_to_file(mesh_1, mesh_2, output_path_features);
	else if (mode == "poisson") {
		double PERCENTAGE_OF_VERTS_TO_CHECK = 0.1;
		size_t N_FEATURES = 3;
		double POISSON_RADIUS = 0.3;
		double DIST_WEIGHT = 1.0;
		double AVERAGING_RADIUS = 0.3;

		if (argc > 7 + n)
			PERCENTAGE_OF_VERTS_TO_CHECK = std::stod(argv[7 + n]);
		if (argc > 8 + n)
			N_FEATURES = std::stoi(argv[8 + n]);
		if (argc > 9 + n)
			POISSON_RADIUS = std::stod(argv[9 + n]);
		if (argc > 10 + n)
			DIST_WEIGHT = std::stod(argv[10 + n]);
		if (argc > 11 + n)
			AVERAGING_RADIUS = std::stod(argv[11 + n]);

		poisson_features_to_file(mesh_1, mesh_2, low_res_1, low_res_2, PERCENTAGE_OF_VERTS_TO_CHECK, N_FEATURES, POISSON_RADIUS, DIST_WEIGHT, AVERAGING_RADIUS, output_path_features, vout);
	}
	else {
		std::cerr << "Unknown 'mode' argument: '" << mode << "'. Should be either 'simple' or 'spatial' or 'poisson'." << std::endl;
		return -1;
	}

	ask_before_termination(DISABLE_USER_INTERACTION);
	return 0;
}