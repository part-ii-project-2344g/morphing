#pragma once

#include "../common/typedefs.h"

#include <unordered_map>

void intersect_asymmetric_deprecated(Polyhedron& mesh_1, Polyhedron& mesh_2,
	std::unordered_map<size_t /*HCI*/, std::vector<size_t /*HCI*/>>& out_intersections_of_edges_1,
	std::unordered_map<size_t /*VCI*/, std::vector<size_t /*VCI*/>>& out_verts_of_containing_face_1,
	std::unordered_map<size_t /*HCI*/, std::unordered_map<size_t /*HCI*/, size_t /*into intersection_points*/>>& out_intersections_data,
	std::vector<Point>& out_intersection_points,
	bool collect_intersections,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2,
	CGAL::Verbose_ostream& vout);

void generate_consistent_intersections_of_edges_2_deprecated(Polyhedron& mesh_1, Polyhedron& mesh_2,
	std::unordered_map<size_t /*HCI*/, std::vector<size_t /*HCI*/>>&       intersections_of_edges_1,
	std::unordered_map<size_t /*HCI*/, std::vector<size_t /*HCI*/>>&       out_intersections_of_edges_2,
	std::unordered_map<size_t /*VCI*/, std::vector<size_t /*VCI*/>>&       out_verts_of_containing_face_2,
	std::unordered_map<size_t /*HCI*/, std::unordered_map<size_t /*HCI*/, size_t /*into intersection_points*/>>& intersections_data,
	std::vector<Point>& intersection_points,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2,
	CGAL::Verbose_ostream& vout);