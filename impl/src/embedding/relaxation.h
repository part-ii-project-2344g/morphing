#pragma once

#include "../common/typedefs.h"

// Assumes that the ordering in points corresponds to ordering in mesh.vertices.
// Returns the largest movement of any point during this relaxation step.
double relaxation_step(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::cpp11::tuple<size_t, size_t, size_t, size_t> fixed_vertices,
	CGAL::Verbose_ostream& vout, std::string& input_filename, bool export_intermediate_mesh);

bool is_embedding_collapsed(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::cpp11::tuple<Point, Point, Point, Point> tetrahedron_vertices, CGAL::Verbose_ostream& vout);