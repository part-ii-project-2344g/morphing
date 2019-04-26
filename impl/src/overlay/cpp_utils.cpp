#include "cpp_utils.h"

HCI get_halfedge(const Polyhedron& mesh, size_t index) {
	auto it = mesh.halfedges_begin();
	for (size_t i = 0; i < index; i++)
		it++;
	return it;
}

VCI get_vertex(const Polyhedron& mesh, size_t index) {
	auto it = mesh.vertices_begin();
	for (size_t i = 0; i < index; i++)
		it++;
	return it;
}

size_t/*vert*/ get_intersection(size_t e_1, size_t e_2, size_t V1, size_t V2,
	std::unordered_map<size_t /*HCI*/, std::unordered_map<size_t /*HCI*/, size_t /*into intersection_points*/>>& intersections_data) {
	return intersections_data[e_1][e_2];
}

size_t/*edge*/ get_opposite(Polyhedron& mesh, size_t e, const HInvIndex& h_index) {
	return h_index[get_halfedge(mesh, e)->opposite()];
}

size_t/*edge*/ get_next(Polyhedron& mesh, size_t e, const HInvIndex& h_index) {
	return h_index[get_halfedge(mesh, e)->next()];
}

size_t/*vert*/ get_end_vert(Polyhedron& mesh, size_t e, const VInvIndex& v_index) {
	return v_index[get_halfedge(mesh, e)->vertex()];
}

size_t/*vert*/ get_start_vert(Polyhedron& mesh, size_t e, const VInvIndex& v_index) {
	return v_index[get_halfedge(mesh, e)->opposite()->vertex()];
}

FCI get_facet(const Polyhedron& mesh, size_t index) {
	auto it = mesh.facets_begin();
	for (size_t i = 0; i < index; i++)
		it++;
	return it;
}

std::string pad_number(size_t n, size_t len) {
	std::string res = std::to_string(n);
	if (res.length() < len) {
		size_t zeros = len - res.length();
		for (size_t i = 0; i < zeros; i++) {
			res = "0" + res;
		}
	}
	return res;
}