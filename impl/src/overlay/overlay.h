#pragma once
#include "../common/typedefs.h"

// Returns impersonation_1, impersonation_2, merged_embedding.
std::tuple<Polyhedron, Polyhedron, Polyhedron> overlay(Polyhedron& mesh_1 /* embedding */, Polyhedron& mesh_2 /* embedding */,
	Polyhedron& shape_1 /* original model */, Polyhedron& shape_2 /* original model */, CGAL::Verbose_ostream& vout);
