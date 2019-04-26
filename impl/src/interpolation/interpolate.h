#pragma once

#include "../common/typedefs.h"

std::string generate_path(const std::string& base_output_path, size_t ind);
void interpolate_and_export(const std::vector<Point>& points_1, const std::vector<Point>& points_2, const std::vector<std::vector<size_t>>& faces, double t, std::string output_path_i, CGAL::Verbose_ostream& vout);