#include"main.h"

#include "../common/typedefs.h"
#include "../common/constants.h"
#include "../embedding/io.h"
#include "../overlay/cpp_utils.h"

#include "interpolate.h"

using namespace std;

#define TEST false

// Example arguments:
// -v 20 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Deer\_dearPLUSbanana[impers_1].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Deer\_dearPLUSbanana[impers_2].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Deer\_dearPLUSbanana[morph].obj"
// -v 10 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\_Simple Tetrahedron (1,2)[impers_1].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\_Simple Tetrahedron (1,2)[impers_2].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Simple Cube\_Simple Tetrahedron (1,2)[morph].obj"

int main_interpolation(int argc, char **argv) {
#if TEST

#else	
	bool verbose = false;
	int n = 0;
	if (strcmp(argv[1], "-v") == 0) { verbose = true; n++; }
	size_t frames = std::stoull(argv[1 + n]);
	std::string input_path_1 = argv[2 + n];
	std::string input_path_2 = argv[3 + n];
	std::string output_path = argv[4 + n];
	std::string output_path_impersonation_triangulated = output_path.substr(0, output_path.find_last_of(".")) + "[triang]" + output_path.substr(output_path.find_last_of("."));
	std::string output_path_impersonation_1 = output_path.substr(0, output_path.find_last_of(".")) + "[impers_1]" + output_path.substr(output_path.find_last_of("."));
	std::string output_path_impersonation_2 = output_path.substr(0, output_path.find_last_of(".")) + "[impers_2]" + output_path.substr(output_path.find_last_of("."));
	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;

	// Read in the impersonations.
	std::vector<Point> points_1, points_2;
	std::vector<std::vector<size_t>> faces_1, faces_2;
	Polyhedron mesh_1, mesh_2;
	bool succ = input_mesh(points_1, faces_1, mesh_1, vout, input_path_1);
	succ = succ && input_mesh(points_2, faces_2, mesh_2, vout, input_path_2);
	if (!succ) {
		cerr << "Error reading input meshes." << endl;
		ask_before_termination(DISABLE_USER_INTERACTION);
		return -1;
	}

	for (size_t i = 0; i < frames; i++) {
		double t = ((double)i) / ((double)(frames - 1));
		std::string output_path_i = generate_path(output_path, i);
		interpolate_and_export(points_1, points_2, faces_1, t, output_path_i, vout);
	}
#endif
	ask_before_termination(DISABLE_USER_INTERACTION);
	return 0;
}