#pragma once
#include "../common/typedefs.h"

// Expose to simple_features.cpp
double curvature_at_vert(const Polyhedron& mesh, const VCI& vert);
double curvature_at_vert(const Polyhedron& mesh, const VCI& vert, double averaging_radius); // Probably doesn't make sense.

bool measure(double& out_result, size_t frames, std::string base_path, CGAL::Verbose_ostream& vout);

bool measure(double& out_result, size_t frames, std::string base_path, CGAL::Verbose_ostream& vout, const std::string& output_path);
