#include "impersonator.h"

#include "../embedding/io.h"
#include "cpp_utils.h"
#include "geometry_utils.h"
#include "mesh_merger.h"
#include <unordered_map>
using namespace std;




// Helper method to generate coordinates of intersection points on edges of the orig_mesh.
std::vector<double> generate_intersection_points_weights(bool orig_mesh_is_mesh_1, size_t V1, size_t V2, size_t V3,
	Polyhedron& mesh_1, Polyhedron& mesh_2,
	const std::unordered_map<size_t /*VCI mesh_1*/, std::vector<size_t /*VCI mesh_2*/>>& verts_of_containing_face_1,
	const std::unordered_map<size_t /*VCI mesh_2*/, std::vector<size_t /*VCI mesh_1*/>>& verts_of_containing_face_2,
	const std::unordered_map<size_t /*HCI mesh_1*/, std::unordered_map<size_t /*HCI mesh_2*/, size_t /*into intersection_points*/>>& intersections_data,
	const std::vector<Point>& intersection_points,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2
) {
	// TODO: Verify the assumption that intersections_data contains every index into intersection_points!
	std::vector<double> result{};
	for (size_t i = 0; i < V3; i++) {
		result.push_back(-1.0); // initialize result
	}
	if (orig_mesh_is_mesh_1) {
		for (auto kvp_outer : intersections_data) {
			size_t e_1 = kvp_outer.first;
			size_t v_1s = get_start_vert(mesh_1, e_1, v_index_1);
			size_t v_1e = get_end_vert(mesh_1, e_1, v_index_1);
			Vector v_1s_pos = get_vertex(mesh_1, v_1s)->point() - CGAL::ORIGIN;
			Vector v_1e_pos = get_vertex(mesh_1, v_1e)->point() - CGAL::ORIGIN;
			for (auto kvp_inner : kvp_outer.second) {
				size_t v_inters = kvp_inner.second;
				Vector v_inters_pos = intersection_points[v_inters] - CGAL::ORIGIN;
				double nominator = acos(v_1s_pos * v_inters_pos / sqrt(v_1s_pos.squared_length() * v_inters_pos.squared_length()));
				double denominator = acos(v_1s_pos * v_1e_pos / sqrt(v_1s_pos.squared_length() * v_1e_pos.squared_length()));
				double weight = nominator / denominator;
				result[v_inters] = weight;
			}
		}
	}
	else { // !orig_mesh_is_mesh_1
		for (auto kvp_outer : intersections_data) {
			for (auto kvp_inner : kvp_outer.second) {
				size_t e_2 = kvp_inner.first;
				size_t v_2s = get_start_vert(mesh_2, e_2, v_index_2);
				size_t v_2e = get_end_vert(mesh_2, e_2, v_index_2);
				Vector v_2s_pos = get_vertex(mesh_2, v_2s)->point() - CGAL::ORIGIN;
				Vector v_2e_pos = get_vertex(mesh_2, v_2e)->point() - CGAL::ORIGIN;

				size_t v_inters = kvp_inner.second;
				Vector v_inters_pos = intersection_points[v_inters] - CGAL::ORIGIN;
				double nominator = acos(v_2s_pos * v_inters_pos / sqrt(v_2s_pos.squared_length() * v_inters_pos.squared_length()));
				double denominator = acos(v_2s_pos * v_2e_pos / sqrt(v_2s_pos.squared_length() * v_2e_pos.squared_length()));
				double weight = nominator / denominator;
				result[v_inters] = weight;
			}
		}
	}

	return result;
}


// Helper method to generate coordinates of other_mesh points on faces of the orig_mesh.
std::vector<std::vector<double>> generate_other_mesh_points_weights(bool orig_mesh_is_mesh_1, size_t V1, size_t V2, size_t V3,
	Polyhedron& mesh_1, Polyhedron& mesh_2,
	const std::unordered_map<size_t /*VCI mesh_1*/, std::vector<size_t /*VCI mesh_2*/>>& verts_of_containing_face_1,
	const std::unordered_map<size_t /*VCI mesh_2*/, std::vector<size_t /*VCI mesh_1*/>>& verts_of_containing_face_2,
	const std::unordered_map<size_t /*HCI mesh_1*/, std::unordered_map<size_t /*HCI mesh_2*/, size_t /*into intersection_points*/>>& intersections_data,
	const std::vector<Point>& intersection_points,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2
) {
	std::vector<std::vector<double>> result{};
	if (orig_mesh_is_mesh_1) {
		for (size_t i = 0; i < V2; i++) {
			result.push_back(std::vector<double>{});
		}
		for (auto kvp : verts_of_containing_face_2) {
			size_t v_2 = kvp.first;
			Vector v_2_pos = get_vertex(mesh_2, v_2)->point() - CGAL::ORIGIN;
			std::vector<size_t> v_1s = kvp.second;
			CGAL_assertion(v_1s.size() == 3);
			Vector v1a = get_vertex(mesh_1, v_1s[0])->point() - CGAL::ORIGIN;
			Vector v1b = get_vertex(mesh_1, v_1s[1])->point() - CGAL::ORIGIN;
			Vector v1c = get_vertex(mesh_1, v_1s[2])->point() - CGAL::ORIGIN;
			std::vector<double> coords = barycentric_coords_spherical(v_2_pos, v1a, v1b, v1c);
			result[v_2] = coords;
		}
	}
	else { // !orig_mesh_is_mesh_1
		for (size_t i = 0; i < V1; i++) {
			result.push_back(std::vector<double>{});
		}
		for (auto kvp : verts_of_containing_face_1) {
			size_t v_1 = kvp.first;
			Vector v_1_pos = get_vertex(mesh_1, v_1)->point() - CGAL::ORIGIN;
			std::vector<size_t> v_2s = kvp.second;
			CGAL_assertion(v_2s.size() == 3);
			Vector v2a = get_vertex(mesh_2, v_2s[0])->point() - CGAL::ORIGIN;
			Vector v2b = get_vertex(mesh_2, v_2s[1])->point() - CGAL::ORIGIN;
			Vector v2c = get_vertex(mesh_2, v_2s[2])->point() - CGAL::ORIGIN;
			std::vector<double> coords = barycentric_coords_spherical(v_1_pos, v2a, v2b, v2c);
			result[v_1] = coords;
		}
	}
	return result;
}


// Impersonation stands for shaping the merged mesh to look like one of the sources.
std::vector<Point> impersonate(bool orig_mesh_is_mesh_1, // do orig_mesh_points correspond to mesh_1 embedding or mesh_2 embedding?
	size_t V1, size_t V2, size_t V3, size_t n_merged_points,
	Polyhedron& mesh_1,	Polyhedron& mesh_2, // mesh_1 and mesh_2 are the spherical embeddings
	std::vector<Point> orig_mesh_points, // points of the original mesh, i.e. the mesh we're recreating, not the embeddings
	const std::unordered_map<size_t /*VCI mesh_1*/, std::vector<size_t /*VCI mesh_2*/>>& verts_of_containing_face_1,
	const std::unordered_map<size_t /*VCI mesh_2*/, std::vector<size_t /*VCI mesh_1*/>>& verts_of_containing_face_2,
	const std::unordered_map<size_t /*HCI mesh_1*/, std::unordered_map<size_t /*HCI mesh_2*/, size_t /*into intersection_points*/>>& intersections_data,
	const std::vector<Point>& intersection_points,
	const std::unordered_map<size_t, size_t>& vert_index_mapping_for_merged_mesh,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2
) { 
	CGAL_assertion(orig_mesh_points.size() == (orig_mesh_is_mesh_1 ? mesh_1 : mesh_2).size_of_vertices());

	// intersection points are on edges of the orig_mesh, a single double is enough to describe the position on an edge
	// ordered as 'intersection_points', weights refer to the adjacent vertex of the halfedge of the orig_mesh from the intersection_data map 
	//  (in the entry whose value is the index of this intersection point)
	std::vector<double> intersection_points_weights = generate_intersection_points_weights(orig_mesh_is_mesh_1, V1, V2, V3, mesh_1, mesh_2, verts_of_containing_face_1, verts_of_containing_face_2, intersections_data, intersection_points, v_index_1, v_index_2);
	
	// points of the other mesh are on faces of the orig_mesh, we need three baricentric coordinates to describe the position on a face
	// ordered as points of the other mesh, baricentric weights inside each vector ordered as in verts_of_containing_face
	std::vector<std::vector<double>> other_mesh_points_weights = generate_other_mesh_points_weights(orig_mesh_is_mesh_1, V1, V2, V3, mesh_1, mesh_2, verts_of_containing_face_1, verts_of_containing_face_2, intersections_data, intersection_points, v_index_1, v_index_2);

	std::vector<Point> result{};
	if (orig_mesh_is_mesh_1) {

		// Counters used to display progress percentage.
		// We do not count the vertices copied over as progress.
		size_t debug_verts_done = 0;
		size_t debug_verts_total = n_merged_points - V1;

		CGAL_assertion(orig_mesh_points.size() == V1);

		// Copy the original points over.
		result.insert(result.end(), orig_mesh_points.begin(), orig_mesh_points.end());

		// Populate the vector with zeros.
		for (size_t i = 0; i < n_merged_points - V1; i++) {
			result.push_back(Point{ 0,0,0 });
		}

		// Handle the other_mesh points.
		for (auto kvp : verts_of_containing_face_2) {
			size_t v_2 = kvp.first;
			std::vector<size_t> v_1s = kvp.second;
			std::vector<Vector> position_1s{};
			for (size_t v_1 : v_1s) {
				Point p_1 = orig_mesh_points[v_1];
				position_1s.push_back(p_1 - CGAL::ORIGIN);
			}
			Vector position_inside_face = Vector{ 0,0,0 };
			for (size_t i = 0; i < position_1s.size(); i++) {
				// TODO: Make sure the other_mesh_points_weights are normalized, i.e. <sum over i> of <other_mesh_points_weights[v_2][i]> is <1>.
				position_inside_face += position_1s[i] * other_mesh_points_weights[v_2][i];
			}
			result[normalize_vertex_index_based_on_overlaps(V1 + v_2, vert_index_mapping_for_merged_mesh)] = CGAL::ORIGIN + position_inside_face;

			// Debug.
			debug_verts_done++;
			for (size_t percentage = 5; percentage <= 100; percentage += 5) {
				if (debug_verts_done == (debug_verts_total / 100)*percentage) {
					std::cout << percentage << "% done..." << std::endl;
				}
			}
		}

		// Handle the intersection points.
		for (auto kvp_outer : intersections_data) {
			size_t e_1 = kvp_outer.first;
			size_t v_1s = get_start_vert(mesh_1, e_1, v_index_1);
			size_t v_1e = get_end_vert(mesh_1, e_1, v_index_1);

			Vector pos_1s = orig_mesh_points[v_1s] - CGAL::ORIGIN;
			Vector pos_1e = orig_mesh_points[v_1e] - CGAL::ORIGIN;

			for (auto kvp_inner : kvp_outer.second) {
				size_t v_inters = kvp_inner.second;
				double v_inters_weight = intersection_points_weights.at(v_inters);
				Vector v_inters_pos = v_inters_weight * pos_1e + (1.0 - v_inters_weight) * pos_1s;
				result[normalize_vertex_index_based_on_overlaps(V1 + V2 + v_inters, vert_index_mapping_for_merged_mesh)] = CGAL::ORIGIN + v_inters_pos;

				// Debug.
				debug_verts_done++;
				for (size_t percentage = 5; percentage <= 100; percentage += 5) {
					if (debug_verts_done == (debug_verts_total / 100)*percentage) {
						std::cout << percentage << "% done..." << std::endl;
					}
				}
			}
		}
	}
	else { // !orig_mesh_is_mesh_1

		// Counters used to display progress percentage.
		// We do not count the vertices copied over as progress.
		size_t debug_verts_done = 0;
		size_t debug_verts_total = n_merged_points - V2;

		CGAL_assertion(orig_mesh_points.size() == V2);

		// Populate the vector fully with zeros.
		for (size_t i = 0; i < n_merged_points; i++) {
			result.push_back(Point{ 0,0,0 });
		}

		// Handle the other_mesh points.
		for (auto kvp : verts_of_containing_face_1) {
			size_t v_1 = kvp.first;
			std::vector<size_t> v_2s = kvp.second;
			std::vector<Vector> position_2s{};
			for (size_t v_2 : v_2s) {
				Point p_2 = orig_mesh_points[v_2];
				position_2s.push_back(p_2 - CGAL::ORIGIN);
			}
			Vector position_inside_face = Vector{ 0,0,0 };
			for (size_t i = 0; i < position_2s.size(); i++) {
				// TODO: Make sure the other_mesh_points_weights are normalized, i.e. <sum over i> of <other_mesh_points_weights[v_2][i]> is <1>.
				position_inside_face += position_2s[i] * other_mesh_points_weights[v_1][i];
			}
			result[normalize_vertex_index_based_on_overlaps(0 + v_1, vert_index_mapping_for_merged_mesh)] = CGAL::ORIGIN + position_inside_face;

			// Debug.
			debug_verts_done++;
			for (size_t percentage = 5; percentage <= 100; percentage += 5) {
				if (debug_verts_done == (debug_verts_total / 100)*percentage) {
					std::cout << percentage << "% done..." << std::endl;
				}
			}
		}


		// Copy the original points over.
		for (size_t i = 0; i < V2; i++) {
			result[normalize_vertex_index_based_on_overlaps(V1 + i, vert_index_mapping_for_merged_mesh)] = orig_mesh_points[i];
		}
		//result.insert(result.end(), orig_mesh_points.begin(), orig_mesh_points.end());}

		// Handle the intersection points.
		for (auto kvp_outer : intersections_data) {
			for (auto kvp_inner : kvp_outer.second) {
				size_t e_2 = kvp_inner.first;
				size_t v_2s = get_start_vert(mesh_2, e_2, v_index_2);
				size_t v_2e = get_end_vert(mesh_2, e_2, v_index_2);
				
				Vector pos_2s = orig_mesh_points[v_2s] - CGAL::ORIGIN;
				Vector pos_2e = orig_mesh_points[v_2e] - CGAL::ORIGIN;

				size_t v_inters = kvp_inner.second;
				double v_inters_weight = intersection_points_weights.at(v_inters);
				Vector v_inters_pos = v_inters_weight * pos_2e + (1.0 - v_inters_weight) * pos_2s;
				result[normalize_vertex_index_based_on_overlaps(V1 + V2 + v_inters, vert_index_mapping_for_merged_mesh)] = CGAL::ORIGIN + v_inters_pos;

				// Debug.
				debug_verts_done++;
				for (size_t percentage = 5; percentage <= 100; percentage += 5) {
					if (debug_verts_done == (debug_verts_total / 100)*percentage) {
						std::cout << percentage << "% done..." << std::endl;
					}
				}
			}
		}
	}

	CGAL_assertion(result.size() == n_merged_points);
	return result;
}