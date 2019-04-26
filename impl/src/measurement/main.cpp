#include"main.h"

#include "../common/typedefs.h"
#include "../common/constants.h"
#include "../embedding/io.h"
#include "../overlay/cpp_utils.h"

#include "measurement.h"

using namespace std;

#define TEST false

// Example arguments:
// -v 20 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\_Models\Deer\_dearPLUSbanana2[morph].obj"

int main_measurement(int argc, char **argv) {
#if TEST

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
		if (!measure(morph_quality, frames, base_path, vout, output_path)) {
			ask_before_termination(DISABLE_USER_INTERACTION);
			return -1;
		}
	}
	else {
		if (!measure(morph_quality, frames, base_path, vout)) {
			ask_before_termination(DISABLE_USER_INTERACTION);
			return -1;
		}
	}
	ask_before_termination(DISABLE_USER_INTERACTION);
#endif
	return 0;
}

double main_measurement_double(int argc, char **argv) {
#if TEST

#else	
	// Parse args.
	bool verbose = false;
	int n = 0;
	if (strcmp(argv[1], "-v") == 0) { verbose = true; n++; }
	size_t frames = std::stoull(argv[1 + n]);
	std::string base_path = argv[2 + n];
	CGAL::Verbose_ostream vout(verbose);
	vout << "\nVerbosity on." << endl;

	double morph_quality = -1.0f;
	if (!measure(morph_quality, frames, base_path, vout)) {
		ask_before_termination(DISABLE_USER_INTERACTION);
		return -1.0;
	}
	ask_before_termination(DISABLE_USER_INTERACTION);
	return morph_quality;
#endif
	return -1.0;
}