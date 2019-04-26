#pragma once
#include "../common/typedefs.h"

std::pair<Polyhedron, Polyhedron> align_features(
	Polyhedron& mesh_1_raw,
	Polyhedron mesh_2, // pass by copy
	std::vector<std::pair<size_t, size_t>> features,
	const std::string& debug_out_path,
	CGAL::Verbose_ostream& vout
);
