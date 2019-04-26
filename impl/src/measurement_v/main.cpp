#include "main.h"

#include "../common/constants.h"
#include "../embedding/io.h"
#include "../overlay/cpp_utils.h"

#include "measurement_v.h"

using namespace std;

#define TEST false

#define N_SAMPLES 20000

int main_measurement_v(int argc, char **argv) {
#if TEST
	Polyhedron p1, p2;
	std::vector<Point> points1, points2;
	std::vector<std::vector<size_t>> faces1, faces2;
	input_mesh(points1, faces1, p1, CGAL::Verbose_ostream{}, std::string{ "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\SimpleShapes\\simple_cube.obj" });
	input_mesh(points2, faces2, p2, CGAL::Verbose_ostream{}, std::string{ "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\SimpleShapes\\simple_cube_half.obj" });
	double correct_res = 0.10938333344328756;
	double test_res = -1;
	std::vector<size_t> n_samples_to_try{ 10, 100, 1000, 10000, 20000, 100000, 1000000 };
	for (size_t n_samples : n_samples_to_try) {
		test_res = sym_diff_volume(p1, p2, n_samples);
		std::cout << "Relative error [n=" << n_samples << "]: " << 100.0 * std::abs(test_res - correct_res) / correct_res << "%." << std::endl;
	}
	std::cout << "Expected value: " << std::setprecision(17) << correct_res << ", test-run value [n=" << n_samples_to_try[n_samples_to_try.size()-1] << "]: " << std::setprecision(17) << test_res << std::endl;
	/*
	Results:
Running a subprogram: measurement_v...
Relative error [n=10]: 24.047%.
Relative error [n=100]: 2.52696%.
Relative error [n=1000]: 0.384575%.
Relative error [n=10000]: 0.324319%.
Relative error [n=20000]: 0.311661%.
Relative error [n=100000]: 0.00607624%.
Relative error [n=1000000]: 0%.
Expected value: 0.10938333344328756, test-run value [n=1000000]: 0.10938333344328756
------------------------------------------
[MASTER] Running time of the 'measurement_v' subprogram: 1.7170000000000001s.
	*/
	return 0;
#else
	// Parse args.
	bool verbose = false;
	int n = 0;
	if (strcmp(argv[1], "-v") == 0) { verbose = true; n++; }
	size_t frames = std::stoull(argv[1 + n]);
	std::string base_path = argv[2 + n];

	std::string output_path = "";
	if (argc == 4 + n) {
		output_path = argv[3 + n];
	}

	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;

	double morph_quality = -1.0f;
	if (output_path.length() > 0) {
		if (!measure_v(morph_quality, N_SAMPLES, frames, base_path, vout, output_path)) {
			ask_before_termination(DISABLE_USER_INTERACTION);
			return -1;
		}
	}
	else {
		if (!measure_v(morph_quality, N_SAMPLES, frames, base_path, vout)) {
			ask_before_termination(DISABLE_USER_INTERACTION);
			return -1;
		}
	}
	ask_before_termination(DISABLE_USER_INTERACTION);
	return 0;
#endif
}

double main_measurement_v_double(int argc, char **argv) {
	// Parse args.
	bool verbose = false;
	int n = 0;
	if (strcmp(argv[1], "-v") == 0) { verbose = true; n++; }
	size_t frames = std::stoull(argv[1 + n]);
	std::string base_path = argv[2 + n];
	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;

	double morph_quality = -1.0f;
	if (!measure_v(morph_quality, N_SAMPLES, frames, base_path, vout)) {
		ask_before_termination(DISABLE_USER_INTERACTION);
		return -1.0;
	}
	ask_before_termination(DISABLE_USER_INTERACTION);
	return morph_quality;
}