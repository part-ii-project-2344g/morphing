#include "features_io.h"
#include "../embedding/io.h"
#include "../overlay/cpp_utils.h"
#include <fstream>

using namespace std;

bool output_features(const std::vector<std::pair<size_t, size_t>>& features, std::string path) {

	// start names of output files with an underscore if they don't already
	path = get_path_with_underscore(path);

	ostream* p_out = &cout;
	ofstream out;
	out.open(path);
	p_out = &out;
	if (!*p_out) {
		cerr << "ERROR: cannot open file for writing: " << path << endl;
		return false;
	}

	*p_out << features.size() << endl;
	for (size_t i = 0; i < features.size(); i++) {
		*p_out << features[i].first << " " << features[i].second << endl;
	}

	out.close();
	std::cout << "Wrote features to file: " << path << endl;
	return true;
}

// TODO: Implement a version that adds (non-machine-readable) feature x,y,z coords at the end of .feat file.
bool output_features(const std::vector<std::pair<size_t, size_t>>& features, std::string path, const Polyhedron& mesh_1, const Polyhedron& mesh_2) {


	// start names of output files with an underscore if they don't already
	path = get_path_with_underscore(path);

	ostream* p_out = &cout;
	ofstream out;
	out.open(path);
	p_out = &out;
	if (!*p_out) {
		cerr << "ERROR: cannot open file for writing: " << path << endl;
		return false;
	}

	*p_out << features.size() << endl;
	for (size_t i = 0; i < features.size(); i++) {
		*p_out << features[i].first << " " << features[i].second << endl;
	}

	// Non-machine-readable
	*p_out << "\n\n";
	for (size_t i = 0; i < features.size(); i++) {
		*p_out << get_vertex(mesh_1, features[i].first)->point() << " : " << get_vertex(mesh_2, features[i].second)->point() << endl;
	}

	out.close();
	std::cout << "Wrote features to file: " << path << endl;
	return true;

}
