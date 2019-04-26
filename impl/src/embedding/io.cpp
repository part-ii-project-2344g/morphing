#include "io.h"

#include <iostream>
#include <fstream>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h> 
#include <boost/algorithm/string.hpp>    

#include "geometry_utils.h"
#include "debug_geometry.h"

using namespace std;

std::string get_path_with_underscore(const std::string& path) {
	if (path[path.find_last_of("\\/") + 1] == '_')
		return path;
	else
		return path.substr(0, path.find_last_of("\\/") + 1) + "_" + path.substr(path.find_last_of("\\/") + 1);
}

std::string get_path_without_underscore(const std::string& path) {
	if (path[path.find_last_of("\\/") + 1] == '_')
		return path.substr(0, path.find_last_of("\\/") + 1) + path.substr(path.find_last_of("\\/") + 2);
	else
		return path;
}

std::string get_path_with_underscore_toggled(const std::string& path) {
	std::string res = get_path_with_underscore(path);
	if (res.compare(path) != 0)
		return res;
	return get_path_without_underscore(path);
}

bool input_mesh_obj(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::Verbose_ostream& vout, std::string& input_filename) {
	// Prepare for input
	std::string name = "cin";
	istream* p_in = &cin;
	ifstream in;
	ifstream in2;
	in.open(input_filename);
	p_in = &in;
	name = input_filename;
	if (!*p_in) {
		in.close();
		std::string path_with_underscore_toggled = get_path_with_underscore_toggled(input_filename);
		cerr << "ERROR: cannot open file for reading: " << input_filename << ".\nInstead attempting: " << path_with_underscore_toggled << "..." << std::endl;
		name = path_with_underscore_toggled;
		in2.open(get_path_with_underscore_toggled(input_filename));
		p_in = &in2;
		if (!*p_in) {
			cerr << "ERROR: cannot open file for reading: " << name << endl;
			exit(1);
		}
	}

	// Read in points and faces.
	vout << "\nReading a Polyhedron from " << name.substr(name.find_last_of("\\/") + 1) << " ..." << endl;
	if (!CGAL::read_OBJ(*p_in, points, faces))
	{
		std::cerr << "Error parsing the OBJ file." << std::endl;
		return false;
	}
	vout << "... done." << endl << endl;

	make_mesh_from_points_and_faces(points, faces, mesh);

	in.close();
	in2.close();
	return true;
}

bool input_mesh_off(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::Verbose_ostream& vout, std::string& input_filename) {	// Prepare for input
	std::ifstream stream(input_filename);
	if (!stream) {
		std::cerr << "ERROR: Cannot open file for reading: " << input_filename << endl;
		return false;
	}

	vout << "\nReading a Polyhedron from " << input_filename.substr(input_filename.find_last_of("\\/") + 1) << " ..." << endl;
	Polyhedron P;
	stream >> P;
	if (!stream) {
		std::cerr << "ERROR: " << input_filename << " is not a polyhedron.";
		return false;
	}
	vout << "... done." << endl << endl;

	mesh = P;
	get_points(points, mesh);
	get_faces(faces, mesh);

	stream.close();
	return true;
}

bool input_mesh(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::Verbose_ostream& vout, std::string& input_filename) {
	points.clear();
	faces.clear();
	std::string format = boost::algorithm::to_lower_copy(input_filename.substr(input_filename.length() - 4));
	if (format == std::string{ ".obj" }) {
		if (!input_mesh_obj(points, faces, mesh, vout, input_filename))
			return false;
	}
	else if (format == std::string{ ".off" }) {
		if (!input_mesh_off(points, faces, mesh, vout, input_filename))
			return false;
	}
	else
		return false;

	// Print information about the mesh.
	vout << "Mesh read from " << input_filename << ":" << endl;
	debug_print_face_types(points, faces, vout);

	return true;
}



bool output_mesh_obj(std::string path, Polyhedron& mesh, CGAL::Verbose_ostream& vout, bool printing_on) {
	std::string name = "cout";
	ostream* p_out = &cout;
	ofstream out;
	out.open(path);
	p_out = &out;
	name = path;
	if (!*p_out) {
		cerr << "ERROR: cannot open file for writing: " << name << endl;
		return false;
	}

	if (printing_on) vout << "\nWriting a Polyhedron to " << name.substr(name.find_last_of("\\/") + 1) << " ..." << endl;
	CGAL::print_polyhedron_wavefront(out, mesh);
	if (printing_on) vout << "... done." << endl << endl;

	if (!*p_out) {
		cerr << "ERROR: while writing file " << name << endl;
		return false;
	}

	out.close();
	return true;
}



bool custom_output_mesh_obj(std::string path, Polyhedron& mesh, CGAL::Verbose_ostream& vout, bool printing_on) {
	std::string name = "cout";
	ostream* p_out = &cout;
	ofstream out;
	out.open(path);
	p_out = &out;
	name = path;
	if (!*p_out) {
		cerr << "ERROR: cannot open file for writing: " << name << endl;
		return false;
	}

	if (printing_on) vout << "\nWriting a Polyhedron to " << name.substr(name.find_last_of("\\/") + 1) << " ..." << endl;

	std::vector<Point> points;
	get_points(points, mesh);
	std::vector<std::vector<size_t>> faces;
	get_faces(faces, mesh);

	*p_out << setprecision(17);
	*p_out << "# verts: " << points.size() << endl;
	*p_out << "# faces: " << faces.size() << endl << endl;
	for (Point p : points) {
		*p_out << "v " << fixed << p.x() << " " << fixed << p.y() << " " << fixed << p.z() << endl;
	}
	*p_out << endl << "# --- " << endl << endl;
	for (std::vector<size_t> f : faces) {
		*p_out << "f  ";
		for (size_t vind : f) {
			*p_out << vind+1 << " ";
		}
		*p_out << endl;
	}

	if (printing_on) vout << "... done." << endl << endl;

	if (!*p_out) {
		cerr << "ERROR: while writing file " << name << endl;
		return false;
	}

	out.close();
	return true;
}

bool custom_output_mesh_obj_from_points_and_faces(std::string path, const std::vector<Point>& points, const std::vector<std::vector<size_t>>& faces, CGAL::Verbose_ostream& vout, bool printing_on) {
	
	// start names of output files with an underscore if they don't already
	path = get_path_with_underscore(path);
	
	std::string name = "cout";
	ostream* p_out = &cout;
	ofstream out;
	out.open(path);
	p_out = &out;
	name = path;
	if (!*p_out) {
		cerr << "ERROR: cannot open file for writing: " << name << endl;
		return false;
	}

	if (printing_on) vout << "\nWriting points and faces to " << name.substr(name.find_last_of("\\/") + 1) << " ..." << endl;

	*p_out << setprecision(17);
	*p_out << "# verts: " << points.size() << endl;
	*p_out << "# faces: " << faces.size() << endl << endl;
	for (Point p : points) {
		*p_out << "v " << fixed << p.x() << " " << fixed << p.y() << " " << fixed << p.z() << endl;
	}
	*p_out << endl << "# --- " << endl << endl;
	for (std::vector<size_t> f : faces) {
		*p_out << "f  ";
		for (size_t vind : f) {
			*p_out << vind + 1 << " ";
		}
		*p_out << endl;
	}

	if (printing_on) vout << "... done." << endl << endl;

	if (!*p_out) {
		cerr << "ERROR: while writing file " << name << endl;
		return false;
	}

	out.close();
	return true;
}


bool output_mesh_off(std::string path, Polyhedron& mesh, CGAL::Verbose_ostream& vout, bool printing_on) {
	std::ofstream stream(path);
	if (!stream) {
		std::cerr << "ERROR: Cannot open file for writing: " << path << endl;
		return false;
	}

	if (printing_on) vout << "\nWriting a Polyhedron to: " << path.substr(path.find_last_of("\\/") + 1) << " ..." << endl;
	stream << mesh;
	if (!stream) {
		std::cerr << "ERROR: while writing file " << path << endl;
		return false;
	}
	if (printing_on) vout << "... done." << endl << endl;

	stream.close();
	return true;
}

bool output_mesh(std::string path, Polyhedron& mesh, CGAL::Verbose_ostream& vout, bool printing_on) {
	std::string format = boost::algorithm::to_lower_copy(path.substr(path.length() - 4));
	// start names of output files with an underscore if they don't already
	path = get_path_with_underscore(path);

	if (format == std::string{ ".obj" }) {
		if (!custom_output_mesh_obj(path, mesh, vout, printing_on))
			return false;
	}
	else if (format == std::string{ ".off" }) {
		if (!output_mesh_off(path, mesh, vout, printing_on))
			return false;
	}
	else
		return false;
	
	return true;
}

void parse_args(const int argc, char** argv, std::string& input_path, std::string& output_path, bool& verbose, double& RELAXATION_EPSILON, Sphere& sph) {
	int n = 0;
	bool help = false;
	for (int i = 1; i < argc; i++) { // check commandline options
		if (strcmp("-v", argv[i]) == 0)
			verbose = true;
		else if ((strcmp("-h", argv[i]) == 0) ||
			(strcmp("-help", argv[i]) == 0))
			help = true;
		else if (strcmp("-s", argv[i]) == 0) { // custom bounding sphere
			sph.radius = std::stod(argv[++i]);
			double x = std::stod(argv[++i]);
			double y = std::stod(argv[++i]);
			double z = std::stod(argv[++i]);
			sph.centre = Point{ x,y,z };
		}
		else if (n == 0) {
			RELAXATION_EPSILON = std::stof(argv[i]);
			n++;
		}
		else if (n == 1) {
			input_path = argv[i];
			n++;
		}
		else if (n == 2) {
			output_path = argv[i];
			n++;
		}
		else {
			++n;
			break;
		}
	}
	if (n == 2) {
		// a default output_path if none was specified
		output_path = input_path.substr(0, input_path.find_last_of(".")) + "[embedding]" + input_path.substr(input_path.find_last_of("."));
		cout << "Set default output path: " << output_path << endl;
		n++;
	}
	if ((n != 3) || help || input_path == output_path) {
		if (!help)
			cerr << "Error: in parameter list" << endl;
		cerr << "Usage: " << argv[0] << " [<options>] <relaxation_epsilon> <input_path> [<output_path>]" << endl;
		cerr << "Option:       -v                            (verbose)." << endl;
		cerr << "Option:       -s <radius> <cx> <cy> <cz>    (custom bounding sphere)." << endl;
		exit(!help);
	}
}

void make_mesh_from_points_and_faces(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh) {
	mesh = Polyhedron{};
	if (!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces)) {
		cerr << "Error: orient_polygon_soup failed." << endl;
		CGAL_assertion(false);
	}
	CGAL_assertion(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(faces));
	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, mesh);
}
void make_mesh_from_points_and_faces_skipping_orientation(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh) {
	mesh = Polyhedron{};
	CGAL_assertion(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(faces));
	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, mesh);
}

void ask_before_termination() {
	std::cout << "Press Enter to close the console." << std::endl;
	std::cin.get();
}
void ask_before_termination(bool disable_user_interaction) {
	if (!disable_user_interaction) ask_before_termination();
}

void get_faces(std::vector<std::vector<size_t>>& faces, const Polyhedron& mesh) {
	VInvIndex index(mesh.vertices_begin(), mesh.vertices_end());
	faces.clear();
	for (FCI fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi) {
		HFCC hc = fi->facet_begin();
		HFCC hc_end = hc;
		std::size_t n = circulator_size(hc);
		std::vector<size_t> face{};
		CGAL_assertion(n >= 3);
		do {
			face.push_back(index[VCI(hc->vertex())]);
			++hc;
		} while (hc != hc_end);
		faces.push_back(face);
	}
}
void get_points(std::vector<Point>& points, const Polyhedron& mesh) {
	points.clear();
	for (VCI vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi) {
		points.push_back(vi->point());
	}
}





// Testing
void test_path_underscore_functions() {
	std::string path_with_underscore = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_Banana2.obj";
	std::string path_without_underscore = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\Banana2.obj";
	
	std::string s1 = get_path_with_underscore(path_with_underscore);
	CGAL_assertion(s1.compare(path_with_underscore) == 0);
	std::string s2 = get_path_with_underscore(path_without_underscore);
	CGAL_assertion(s2.compare(path_with_underscore) == 0);
	std::string s3 = get_path_without_underscore(path_without_underscore);
	CGAL_assertion(s3.compare(path_without_underscore) == 0);
	std::string s4 = get_path_without_underscore(path_with_underscore);
	CGAL_assertion(s4.compare(path_without_underscore) == 0);
}