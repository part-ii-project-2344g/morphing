#pragma once
#include "../common/typedefs.h"
#include <unordered_map>

std::vector<Point> impersonate(bool orig_mesh_is_mesh_1, // do orig_mesh_points correspond to mesh_1 embedding or mesh_2 embedding?
	size_t V1, size_t V2, size_t V3, size_t n_merged_points,
	Polyhedron& mesh_1, Polyhedron& mesh_2, // mesh_1 and mesh_2 are the spherical embeddings
	std::vector<Point> orig_mesh_points, // points of the original mesh, i.e. the mesh we're recreating, not the embeddings
	const std::unordered_map<size_t /*VCI mesh_1*/, std::vector<size_t /*VCI mesh_2*/>>& verts_of_containing_face_1,
	const std::unordered_map<size_t /*VCI mesh_2*/, std::vector<size_t /*VCI mesh_1*/>>& verts_of_containing_face_2,
	const std::unordered_map<size_t /*HCI mesh_1*/, std::unordered_map<size_t /*HCI mesh_2*/, size_t /*into intersection_points*/>>& intersections_data,
	const std::vector<Point>& intersection_points,
	const std::unordered_map<size_t, size_t>& vert_index_mapping_for_merged_mesh,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2
);