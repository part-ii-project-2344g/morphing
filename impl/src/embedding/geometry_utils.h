#pragma once

#include "../common/typedefs.h"

Transformation rotX(double angle);
Transformation rotY(double angle);
Transformation rotZ(double angle);

Vector normalized(const Vector& v);
Point normalized(const Point& p);

void normalize_all_points(std::vector<Point>& points);
bool is_orientation_consistent_on_sphere(const std::vector<Point>& points, const std::vector<std::vector<size_t>>& faces);

size_t find_nearest_neighbour_of_point(std::vector<Point>& points, const Point& target);
Point diametric(const Point& original);
CGAL::cpp11::tuple<Point, Point, Point, Point> diametric(CGAL::cpp11::tuple<Point, Point, Point, Point> tetrahedron_vertices); 
void triangulate_mesh_if_necessary(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::Verbose_ostream& vout, std::string& input_filename);
