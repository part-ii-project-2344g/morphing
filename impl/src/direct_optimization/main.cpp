#include "main.h"
#include "direct_optimization.h"
#include "../common/constants.h"
#include "../embedding/io.h"

#include "../overlay/overlay.h"
#include "../interpolation/interpolate.h"
#include "../measurement/measurement.h"

using namespace std;
int main_direct_optimization(int argc, char **argv){

	bool verbose = false;
	int n = 0;
	if (strcmp(argv[1], "-v") == 0) { verbose = true; n++; }
	std::string input_path_1 = argv[1 + n];
	std::string input_path_1_embedding = argv[2 + n];
	std::string input_path_2 = argv[3 + n];
	std::string input_path_2_embedding = argv[4 + n];
	std::string output_path_base = argv[5 + n];
	std::string output_path_1 = argv[6 + n];
	std::string output_path_2 = argv[7 + n];
	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;

	std::vector<Point> embedding_points_1, embedding_points_2;
	std::vector<std::vector<size_t>> embedding_faces_1, embedding_faces_2;
	Polyhedron embedding_1, embedding_2;
	std::vector<Point> shape_points_1, shape_points_2;
	std::vector<std::vector<size_t>> shape_faces_1, shape_faces_2;
	Polyhedron shape_1, shape_2;

	bool succ = input_mesh(embedding_points_1, embedding_faces_1, embedding_1, vout, input_path_1_embedding);
	succ = succ && input_mesh(embedding_points_2, embedding_faces_2, embedding_2, vout, input_path_2_embedding);
	succ = succ && input_mesh(shape_points_1, shape_faces_1, shape_1, vout, input_path_1);
	succ = succ && input_mesh(shape_points_2, shape_faces_2, shape_2, vout, input_path_2);
	if (!succ) {
		cerr << "Error reading input meshes." << endl;
		ask_before_termination(DISABLE_USER_INTERACTION);
		return -1;
	}

	// This function needs to take two embeddings and the original meshes, 
	// then perform overlay, interpolation, and measurement, and return the measurement result.
	std::function<double(Polyhedron&, Polyhedron&, Polyhedron&, Polyhedron&)> f_to_minimize = [&vout, &output_path_base](Polyhedron& emb_1, Polyhedron& emb_2, Polyhedron& shape_1, Polyhedron& shape_2 ) {
		
		// 1/3 Overlay.
		std::tuple<Polyhedron, Polyhedron, Polyhedron> overlay_results = overlay(emb_1, emb_2, shape_1, shape_2, vout);
		Polyhedron impers_1 = std::get<0>(overlay_results);
		Polyhedron impers_2 = std::get<1>(overlay_results);


		// 2/3 Interpolate.
		std::vector<Point> points_1, points_2;
		std::vector<std::vector<size_t>> faces_1;
		get_points(points_1, impers_1);
		get_points(points_2, impers_2);
		get_faces(faces_1, impers_1);
		size_t frames = 40;
		for (size_t i = 0; i < frames; i++) {
			double t = ((double)i) / ((double)(frames - 1));
			std::string output_path_i = generate_path(output_path_base, i);
			interpolate_and_export(points_1, points_2, faces_1, t, output_path_i, vout);
		}


		// 3/3 Measure.
		double measurement_result = 0.0;
		measure(measurement_result, frames, output_path_base, vout);
		
		return measurement_result;
	};

	std::pair<Polyhedron, Polyhedron> embeddings_optimized = optimize_embedding(embedding_1, embedding_2, shape_1, shape_2, f_to_minimize);

	output_mesh(output_path_1, embeddings_optimized.first,  vout, false);
	output_mesh(output_path_2, embeddings_optimized.second, vout, false);
	ask_before_termination(DISABLE_USER_INTERACTION);
	return 0;
}
