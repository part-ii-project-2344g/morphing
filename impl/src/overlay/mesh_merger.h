#pragma once
#include "../common/typedefs.h"
#include "intersections_new.h"
#include<unordered_map>

void generate_vv_overlaps(
	Polyhedron& mesh_1, Polyhedron& mesh_2,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& out_vertex_overlap_map_1,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& out_vertex_overlap_map_2,
	const size_t V1, // vertex count of mesh_1
	const size_t V2, // vertex count of mesh_2
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	CGAL::Verbose_ostream& vout);

class NormalizedVertexList {
public:
	NormalizedVertexList(std::vector<size_t> verts) : face(verts) {}
	std::vector<size_t> get_face() { return face; }
	size_t hash() {
		// cstheory.stackexchange.com/a/3403
		size_t S = face.size();
		size_t miniind = 0;
		size_t mini = face[0];
		for (size_t i = 0; i < S; i++) {
			if (face[i] < mini) {
				mini = face[i];
				miniind = i;
			}
		}
		size_t res = 0;
		size_t mul = 31;
		for (size_t i = 0; i < S; i++) {
			size_t ind = (miniind + i) % S;
			res = (res + face[ind] * mul) % 2147483647;
			mul = (mul * 31) % 2147483647;
		}
		return res;
	}
	bool equals(const NormalizedVertexList& other) {
		size_t S = face.size();
		if (S != other.face.size())
			return false;
		size_t miniind = 0;
		size_t mini = face[0];
		for (size_t i = 0; i < S; i++) {
			if (face[i] < mini) {
				mini = face[i];
				miniind = i;
			}
		}
		size_t miniind_ = 0;
		size_t mini_ = other.face[0];
		for (size_t i = 0; i < S; i++) {
			if (other.face[i] < mini_) {
				mini_ = other.face[i];
				miniind_ = i;
			}
		}

		for (size_t i = 0; i < S; i++) {
			size_t ind = (i + miniind) % S;
			size_t ind_ = (i + miniind_) % S;
			if (face[ind] != other.face[ind_])
				return false;
		}
		return true;
	}
private:
	std::vector<size_t> face;
};

std::unordered_map<size_t, size_t> generate_vert_index_mapping_for_merged_mesh(
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	const std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>>& vertex_overlap_map_2,
	const size_t V1, // original vertex count of mesh_1
	const size_t V2, // original vertex count of mesh_2
	const size_t V3  // count of intersection vertices
);

// New implementation based on IntersectionEvents.
std::vector <std::vector<size_t>> generate_merged_faces_from_events(
	Polyhedron& mesh_1, Polyhedron& mesh_2,
	const std::unordered_map<size_t, size_t>& vert_index_mapping_for_merged_mesh,
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
	CGAL::Verbose_ostream& vout);


std::vector<Point> generate_merged_points(const std::vector<Point>&, const std::unordered_map<size_t, size_t>& vert_index_mapping_for_merged_mesh);

// Exposed only to be tested in 'test_overlay.cpp'.
size_t get_neighbouring_halfedge_by_angle(const Vector& normal, const std::vector<std::pair<Vector, size_t>>& candidate_dirs, const Vector& reference_dir, bool go_the_other_way);

//
//
//
//
//
// Deprecated...
//
//
//
//
//

std::vector <std::vector<size_t>> generate_merged_faces_deprecated(
	Polyhedron& mesh_1, Polyhedron& mesh_2,
	std::unordered_map<size_t /*hedge*/, std::vector<size_t /*hedge*/>>& intersections_of_edges_1,
	std::unordered_map<size_t /*hedge*/, std::vector<size_t /*hedge*/>>& intersections_of_edges_2,
	std::unordered_map<size_t /*hedge*/, std::unordered_map<size_t /*hedge*/, size_t /*into intersection_points*/>>& intersections_data,
	const size_t V1, // vertex count of mesh_1
	const size_t V2, // vertex count of mesh_2
	const size_t V3, // count of intersection vertices
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	CGAL::Verbose_ostream& vout);


// Expose to impersonator.
size_t normalize_vertex_index_based_on_overlaps(size_t v_ind, const std::unordered_map<size_t, size_t> vert_index_mapping_for_merged_mesh);
