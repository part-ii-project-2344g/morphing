#pragma once
#include "../common/typedefs.h"

double low_res_curvature_at_vert(const VCI& vert, const std::vector<std::pair<Point, double>>& low_res_curvatures, const double radius, const double sqrt_radius, CGAL::Verbose_ostream& vout);
Polyhedron low_res_mesh(const Polyhedron& mesh, size_t n_edges, CGAL::Verbose_ostream& vout);
std::vector<std::pair<Point, double>> get_low_res_curvatures(const Polyhedron& low_res);