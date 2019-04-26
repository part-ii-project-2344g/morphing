#include "debug_geometry.h"

using namespace std;

void debug_print_face_types(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, CGAL::Verbose_ostream& vout) {
	vout << "points.size() = " << points.size() << endl;
	vout << "faces.size() = " << faces.size() << endl;
	std::map<size_t, size_t> counts{};
	for (std::vector<size_t> face : faces) {
		size_t verts_of_face = face.size();
		if (counts.find(verts_of_face) == counts.end()) counts[verts_of_face] = 0;
		counts[verts_of_face]++;
	}
	for (auto kvp : counts) {
		size_t verts_of_face = kvp.first;
		size_t count = kvp.second;
		std::string ngon_name[7]{ "", "", "",
			"tris:        ",
			"quads:       ",
			"pentagons:   ",
			"hexagons:    " };
		CGAL_assertion(verts_of_face >= 3);
		if (verts_of_face <= 6)
			vout << "Number of " << ngon_name[verts_of_face] << count << "." << endl;
		else
			vout << "Number of " << verts_of_face << "-gonal faces: " << count << "." << endl;
	}
	vout << endl;
}

void debug_print_faces_and_points(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, CGAL::Verbose_ostream& vout) {
	if (points.size() + faces.size() < 40) {
		int i = 0;
		for (Point p : points) {
			vout << "points[" << i++ << "] = " << p << endl;
		} i = 0;
		for (std::vector<size_t> face : faces) {
			vout << "faces[" << i++ << "] = {";
			for (int i = 0; i < face.size(); i++) {
				vout << face[i];
				if (i != face.size() - 1) vout << ", ";
			}
			vout << "}" << endl;
		}
	}
}