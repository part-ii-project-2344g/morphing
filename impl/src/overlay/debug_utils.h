#pragma once
#include "../common/typedefs.h"
#include "intersections_new.h" // IntersectionEventStore
#include <unordered_map>

void debug_print_intersections(
	CGAL::Verbose_ostream& vout,
	const std::vector<Point>& intersection_points,
	std::unordered_map<size_t, std::vector<size_t>> intersections_of_edges_1,
	std::unordered_map<size_t, std::vector<size_t>> intersections_of_edges_2,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2,
	size_t V1,
	size_t V2,
	const Polyhedron& mesh_1,
	const Polyhedron& mesh_2,
	std::unordered_map<size_t, std::unordered_map<size_t, size_t>> intersections_data);

void debug_print_hedges(Polyhedron& mesh_1, Polyhedron& mesh_2,
	const VInvIndex& v_index_1, const HInvIndex& h_index_1,
	const VInvIndex& v_index_2, const HInvIndex& h_index_2, size_t V1, CGAL::Verbose_ostream& vout);

void debug_print_intersection_events(size_t debug_counter1, const IntersectionEventStore& intersection_events, CGAL::Verbose_ostream& vout);

bool debug_test_intersection_consistency(Polyhedron& mesh_1, Polyhedron& mesh_2, const IntersectionEventStore& intersection_events,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_1,
	const HInvIndex& h_index_2,
	CGAL::Verbose_ostream& vout); // check whether we magically 'get out of a triangle without crossing its edge or vertex' at any time.

void debug_everything(
	Polyhedron& mesh_1, Polyhedron& mesh_2,
	const std::unordered_map<size_t, size_t>& vert_index_mapping_for_merged_mesh,
	const std::vector<Point>& intersection_points,
	const IntersectionEventStore& intersection_events,
	std::unordered_map<size_t /*hedge*/, std::unordered_map<size_t /*hedge*/, size_t /*into intersection_points*/>>& intersections_data,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	const std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>>& vertex_overlap_map_2,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 hedges*/>>& hedges_through_vertex_1,
	const std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 hedges*/>>& hedges_through_vertex_2,
	const size_t V1, // vertex count of mesh_1
	const size_t V2, // vertex count of mesh_2
	const size_t V3, // count of intersection vertices
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	CGAL::Verbose_ostream& vout
);