#include "random_tetrahedron.h"

#include "geometry_utils.h"

using namespace std;

CGAL::cpp11::tuple<Point, Point, Point, Point> choose_random_tetrahedron(std::vector<Point>& points, unsigned int seed) {
	Point a{ sqrt(8.0 / 9.0),   0.0,             -1.0 / 3.0 };
	Point b{ -sqrt(2.0 / 9.0),  sqrt(2.0 / 3.0), -1.0 / 3.0 };
	Point c{ -sqrt(2.0 / 9.0), -sqrt(2.0 / 3.0), -1.0 / 3.0 };
	Point d{ 0.0,               0.0,              1.0 };

	Random rand{ seed };
	Transformation rx = rotX(rand.get_double(0.0, 2 * CGAL_PI));
	Transformation ry = rotY(rand.get_double(0.0, 2 * CGAL_PI));
	Transformation rz = rotZ(rand.get_double(0.0, 2 * CGAL_PI));
	Transformation r = rx * ry * rz;

	Point tetrahedron_a = a.transform(r);
	Point tetrahedron_b = b.transform(r);
	Point tetrahedron_c = c.transform(r);
	Point tetrahedron_d = d.transform(r);

	return CGAL::cpp11::make_tuple(
		tetrahedron_a,
		tetrahedron_b,
		tetrahedron_c,
		tetrahedron_d
	);
}

CGAL::cpp11::tuple<size_t, size_t, size_t, size_t> choose_fixed_vertices(std::vector<Point>& points, CGAL::cpp11::tuple<Point, Point, Point, Point> tetrahedron_vertices) {
	return CGAL::cpp11::make_tuple(
		find_nearest_neighbour_of_point(points, CGAL::cpp11::get<0>(tetrahedron_vertices)),
		find_nearest_neighbour_of_point(points, CGAL::cpp11::get<1>(tetrahedron_vertices)),
		find_nearest_neighbour_of_point(points, CGAL::cpp11::get<2>(tetrahedron_vertices)),
		find_nearest_neighbour_of_point(points, CGAL::cpp11::get<3>(tetrahedron_vertices))
	);
}