#pragma once
#include "../common/typedefs.h"

std::pair<Polyhedron, Polyhedron> optimize_embedding(
	Polyhedron& mesh_1, // the embedding to be modified
	Polyhedron& mesh_2,
	Polyhedron& shape_1,
	Polyhedron& shape_2,
	std::function<double(Polyhedron&, Polyhedron&, Polyhedron&, Polyhedron&)> f_to_minimize
);
