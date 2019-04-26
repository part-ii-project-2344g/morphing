#include"parse_features.h"
#include<fstream>
#include"../embedding/io.h"
using namespace std;

std::vector<std::pair<size_t, size_t>> read_features(const std::string& path) {

	istream* p_in = &std::cin;
	ifstream in;
	ifstream in2;
	std::string name = path;
	in.open(path);
	p_in = &in;
	if (!*p_in) {
		in.close();
		std::string path_with_underscore_toggled = get_path_with_underscore_toggled(path);
		name = path_with_underscore_toggled;
		cerr << "ERROR: cannot open file for reading: " << path << ".\nInstead attempting: " << path_with_underscore_toggled << "..." << std::endl;
		in2.open(get_path_with_underscore_toggled(path));
		p_in = &in2;
		if (!*p_in) {
			cerr << "ERROR: cannot open file for reading: " << path << endl;
			exit(1);
		}
	}

	size_t n_features;
	*p_in >> n_features;

	cout << "Reading " << n_features << " features..." << endl;

	std::vector<std::pair<size_t, size_t>> res{};
	for (size_t i = 0; i < n_features; i++) {
		size_t v1, v2;
		*p_in >> v1 >> v2;
		res.push_back(make_pair(v1, v2));
	}

	cout << "Read " << res.size() << " features from " << name << "." << endl;

	in.close();
	return res;
}

void write_features(const std::string& path, const std::vector<std::pair<size_t, size_t>>& features) {
	ofstream out;
	out.open(path);
	ostream* p_out = &out;
	if (!*p_out) {
		cerr << "ERROR: cannot open file for writing: " << path << endl;
		exit(1);
	}

	size_t n_features = features.size();
	*p_out << n_features << endl;

	for (auto pair : features) {
		size_t v1 = pair.first, v2 = pair.second;
		*p_out << v1 << " " << v2 << endl;
	}

	out.close();
}