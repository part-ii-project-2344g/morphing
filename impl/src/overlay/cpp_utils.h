#pragma once

#include "../common/typedefs.h"
#include <deque>
#include <unordered_map>
#include <unordered_set>

template<typename T>
T stack_pop(std::deque<T>& stack) {
	T ret = stack.front();
	stack.pop_front();
	return ret;
}

template<typename T, typename U>
void insert_into_map(std::unordered_map<T, std::vector<U>>& m, const T& key, const U& val) {
	if (m.find(key) == m.end()) {
		m[key] = std::vector<U>{};
	}
	m[key].push_back(val);
}

template<typename T, typename U, typename V>
void insert_into_map(std::unordered_map<T, std::unordered_map<U, V>>& m, const T& key_outer, const U& key_inner, const V& val) {
	if (m.find(key_outer) == m.end()) {
		m[key_outer] = std::unordered_map<U, V>{};
	}
	m[key_outer][key_inner] = val;
}

template<typename T>
bool set_contains(const std::unordered_set<T>& s, const T& key) {
	return s.find(key) != s.end();
}


HCI get_halfedge(const Polyhedron& mesh, size_t index);
VCI get_vertex(const Polyhedron& mesh, size_t index);
Facet get_facet(const Polyhedron& mesh, size_t index);

size_t/*edge*/ get_opposite(Polyhedron& mesh, size_t e, const HInvIndex& h_index);
size_t/*edge*/ get_next(Polyhedron& mesh, size_t e, const HInvIndex& h_index);

size_t/*vert*/ get_end_vert(Polyhedron& mesh, size_t e, const VInvIndex& v_index);
size_t/*vert*/ get_start_vert(Polyhedron& mesh, size_t e, const VInvIndex& v_index);

size_t/*vert*/ get_intersection(size_t e_1, size_t e_2, size_t V1, size_t V2,
	std::unordered_map<size_t /*HCI*/, std::unordered_map<size_t /*HCI*/, size_t /*into intersection_points*/>>& intersections_data);

std::string pad_number(size_t n, size_t len);