#pragma once
#include "../common/typedefs.h"

CGAL::cpp11::tuple<Point, Point, Point, Point> choose_random_tetrahedron(std::vector<Point>& points, unsigned int seed);
CGAL::cpp11::tuple<size_t, size_t, size_t, size_t> choose_fixed_vertices(std::vector<Point>& points, CGAL::cpp11::tuple<Point, Point, Point, Point> tetrahedron_vertices);
