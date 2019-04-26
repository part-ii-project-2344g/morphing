#pragma once
#include "../common/typedefs.h"

void debug_print_face_types(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, CGAL::Verbose_ostream& vout);
void debug_print_faces_and_points(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, CGAL::Verbose_ostream& vout);