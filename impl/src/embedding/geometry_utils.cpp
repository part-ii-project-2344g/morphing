#include "geometry_utils.h"
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <CGAL/basic.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include "debug_geometry.h"
#include "io.h"

using namespace std;

bool is_orientation_consistent_on_sphere(const std::vector<Point>& points, const std::vector<std::vector<size_t>>& faces) {
	bool orientation_so_far = false;
	bool orientation_so_far_valid = false;

	int fail_counter = 0;

	for (std::vector<size_t> face : faces) {
		CGAL_assertion(face.size() == 3);
		Vector a, b, c;
		a = points[face[0]] - CGAL::ORIGIN;
		b = points[face[1]] - CGAL::ORIGIN;
		c = points[face[2]] - CGAL::ORIGIN;
		bool orientation = (CGAL::cross_product(a, b) * c) > 0;
		if (orientation_so_far_valid) {
			if (orientation != orientation_so_far) {
				fail_counter++;
			}
		}
		else {
			orientation_so_far = orientation;
			orientation_so_far_valid = true;
		}
	}

	if (fail_counter > 440) {
		cerr << min(fail_counter, (int)faces.size() - fail_counter) << "/" << faces.size() <<  " faces have a different orientation." << endl;
		return false;
	}
	return true;
}

void triangulate_mesh_if_necessary(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::Verbose_ostream& vout, std::string& input_filename) {
	if (!CGAL::is_triangle_mesh(mesh)) {
		// Triangulate all faces.
		CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
		CGAL_assertion(CGAL::is_triangle_mesh(mesh));

		// Recover faces of the triangulated Polyhedron.
		get_faces(faces, mesh);
		get_points(points, mesh);

		// debug
		vout << "\nAfter triangulation:" << endl;
		debug_print_face_types(points, faces, vout);

		// save triangulated mesh for inspection
		//std::string input_name{ input_filename };
		//output_mesh(input_name.substr(0, input_name.length() - 4) + "[triangulated]" + input_name.substr(input_name.length()-4), mesh, vout);
	}
	else {
		vout << "Mesh was already triangulated.\n";
		vout << faces[0].size() << " == 3 :)" << endl;
	}
}

size_t find_nearest_neighbour_of_point(std::vector<Point>& points, const Point& target) {
	double best = std::numeric_limits<double>::max();
	size_t result_index = 0;
	for (int i = 0; i < points.size(); i++) {
		Point p = points[i];
		double dist = (p - target).squared_length();
		if (dist < best) {
			result_index = i;
			best = dist;
		}
	}
	return result_index;
}

Point diametric(const Point& original) {
	return CGAL::ORIGIN + (CGAL::ORIGIN - original);
}
CGAL::cpp11::tuple<Point, Point, Point, Point> diametric(CGAL::cpp11::tuple<Point, Point, Point, Point> tetrahedron_vertices) {
	return CGAL::cpp11::make_tuple(
		diametric(CGAL::cpp11::get<0>(tetrahedron_vertices)),
		diametric(CGAL::cpp11::get<1>(tetrahedron_vertices)),
		diametric(CGAL::cpp11::get<2>(tetrahedron_vertices)),
		diametric(CGAL::cpp11::get<3>(tetrahedron_vertices))
	);
}

Transformation rotX(double angle) {
	const double cosa{ cos(angle) };
	const double sina{ sin(angle) };
	return Transformation(
		1.0, 0.0, 0.0,
		0.0, cosa, -sina,
		0.0, sina, cosa);
}
Transformation rotY(double angle) {
	const double cosa{ cos(angle) };
	const double sina{ sin(angle) };
	return Transformation(
		cosa, 0.0, sina,
		0.0, 1.0, 0.0,
		-sina, 0.0, cosa);
}
Transformation rotZ(double angle) {
	const double cosa{ cos(angle) };
	const double sina{ sin(angle) };
	return Transformation(
		cosa, -sina, 0.0,
		sina, cosa, 0.0,
		0.0, 0.0, 1.0);
}

Vector normalized(const Vector& v) {
	if (sqrt(v.squared_length()) == 0.0) return v;
	Vector res = v / sqrt(v.squared_length());
	CGAL_assertion(CGAL::abs(res.squared_length() - 1.0) < 0.001);
	return res;
}

Point normalized(const Point& p) {
	Vector v = normalized(p - CGAL::ORIGIN); return CGAL::ORIGIN + v;
}

void normalize_all_points(std::vector<Point>& points) {
	for (int i = 0; i < points.size(); i++) {
		points[i] = normalized(points[i]);
	}

	// check
	for (Point p : points) {
		CGAL_assertion(abs((p - CGAL::ORIGIN).squared_length() - 1.0) < 0.001);
	}
}



/*
bool is_orientation_consistent(Polyhedron& mesh) {
bool orientation_so_far = false;
bool orientation_so_far_valid = false;


for (FCI fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi) {
HFCC hc = fi->facet_begin();
HFCC hc_end = hc;
std::size_t n = circulator_size(hc);
CGAL_assertion(n == 3);
Vector a, b, c;
a = (hc++)->vertex()->point() - CGAL::ORIGIN;
b = (hc++)->vertex()->point() - CGAL::ORIGIN;
c = (hc++)->vertex()->point() - CGAL::ORIGIN;
CGAL_assertion(hc == hc_end);
bool orientation = (CGAL::cross_product(a, b) * c) > 0;
if (orientation_so_far_valid) {
if (orientation != orientation_so_far)
return false;
}
else {
orientation_so_far = orientation;
orientation_so_far_valid = true;
}
}
return true;
}
*/