#include "overlay.h"

#include "../common/constants.h"
#include "../embedding/io.h"
#include "../embedding/debug_geometry.h"

#include "debug_utils.h"
#include "impersonator.h"
#include "intersections_deprecated.h"
#include "intersections_new.h"
#include "mesh_merger.h"
#include "test_overlay.h"

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>


#define USE_DEPRECATED_INTERSECTION_METHOD false
#define DEBUG_EVERYTHING false

// Returns impersonation_1, impersonation_2, merged_embedding.
std::tuple<Polyhedron, Polyhedron, Polyhedron> overlay(Polyhedron& mesh_1 /* embedding */, Polyhedron& mesh_2 /* embedding */, 
	Polyhedron& shape_1 /* original model */, Polyhedron& shape_2 /* original model */, CGAL::Verbose_ostream& vout) {
	

	std::vector<Point> points_1, points_2;
	std::vector<std::vector<size_t>> faces_1, faces_2;
	std::vector<Point> original_mesh_points_1, original_mesh_points_2;             // passed to 'impersonate()'

	get_points(original_mesh_points_1, shape_1);
	get_points(original_mesh_points_2, shape_2);
	get_points(points_1, mesh_1);
	get_points(points_2, mesh_2);
	get_faces(faces_1, mesh_1);
	get_faces(faces_2, mesh_2);

	// Map a halfedge H to vector of halfedges of the other mesh that intersect H
	std::unordered_map<size_t /*mesh_1 edges*/, std::vector<size_t /*mesh_2 edges*/>> intersections_of_edges_1;
	std::unordered_map<size_t /*mesh_2 edges*/, std::vector<size_t /*mesh_1 edges*/>> intersections_of_edges_2;

	// Map a vertex V to a vector of 3 verts of the other mesh that form a triangle containing V
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>> verts_of_containing_face_1;
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>> verts_of_containing_face_2;

	// Store intersection points of edges of two meshes, together with the generating edges
	std::unordered_map<size_t /*mesh_1 edges*/, std::unordered_map<size_t /*mesh_2 edges*/, size_t /*into intersection_points*/>> intersections_data;

	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>> vertex_overlap_map_1;
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>> vertex_overlap_map_2;

	// Store hedges of one mesh passing through vertex of the other mesh ('passing through' excludes vertex-vertex overlaps).
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 hedges*/>> hedges_through_vertex_1;
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 hedges*/>> hedges_through_vertex_2;

	// Lists generated intersection points
	std::vector<Point> intersection_points;

	// Index vertices and halfedges of mesh_1 and mesh_2
	HInvIndex h_index_1(mesh_1.halfedges_begin(), mesh_1.halfedges_end());
	HInvIndex h_index_2(mesh_2.halfedges_begin(), mesh_2.halfedges_end());
	VInvIndex v_index_1(mesh_1.vertices_begin(), mesh_1.vertices_end());
	VInvIndex v_index_2(mesh_2.vertices_begin(), mesh_2.vertices_end());
	FInvIndex f_index_1(mesh_1.facets_begin(), mesh_1.facets_end());
	FInvIndex f_index_2(mesh_2.facets_begin(), mesh_2.facets_end());

	get_points(points_1, mesh_1);
	get_points(points_2, mesh_2);
	size_t V1 = points_1.size();
	size_t V2 = points_2.size();

	// Collect the vertex-vertex overlaps.
	generate_vv_overlaps(mesh_1, mesh_2, vertex_overlap_map_1, vertex_overlap_map_2, V1, V2, v_index_1, h_index_1, v_index_2, h_index_2, vout);

	// Collect the intersections.
#if USE_DEPRECATED_INTERSECTION_METHOD
	// Alternative approach, works for simple tetrahedra.
	// This approach generates mutually consistent intersection lists, because geometry test are only conducted once, and then the dual intersection lists are generated based on the same results.
	intersect_asymmetric_deprecated(mesh_1, mesh_2, intersections_of_edges_1, verts_of_containing_face_1, intersections_data, intersection_points, true, v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout);
	generate_consistent_intersections_of_edges_2_deprecated(mesh_1, mesh_2, intersections_of_edges_1, intersections_of_edges_2, verts_of_containing_face_2, intersections_data, intersection_points, v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout);
	size_t V3 = intersection_points.size();

	// Even older approach, also works for simple tetrahedra.
//	intersect_asymmetric_deprecated(mesh_2, mesh_1, intersections_of_edges_2, verts_of_containing_face_2, intersections_data, intersection_points, false, v_index_2, h_index_2, f_index_2, v_index_1, h_index_1, f_index_1, vout);
//	intersect_asymmetric_deprecated(mesh_1, mesh_2, intersections_of_edges_1, verts_of_containing_face_1, intersections_data, intersection_points, true, v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout);
#else 
	// A new implementation taking into account vertex-vertex overlaps and vertex-edge intersections.
	IntersectionEventStore intersection_events{};
	collect_intersection_events(mesh_1, mesh_2, intersection_events, intersections_data, intersection_points, vertex_overlap_map_1, vertex_overlap_map_2,
		hedges_through_vertex_1, hedges_through_vertex_2, v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout);
	size_t V3 = intersection_points.size();
	std::unordered_map<size_t, size_t> vert_index_mapping_for_merged_mesh = generate_vert_index_mapping_for_merged_mesh(
		vertex_overlap_map_1, vertex_overlap_map_2, V1, V2, V3);
#if DEBUG_EVERYTHING
	// Debug.
	debug_everything(mesh_1, mesh_2, vert_index_mapping_for_merged_mesh, intersection_points, intersection_events, intersections_data, vertex_overlap_map_1, vertex_overlap_map_2,
		hedges_through_vertex_1, hedges_through_vertex_2, V1, V2, V3, v_index_1, h_index_1, v_index_2, h_index_2, vout);
#endif
	CGAL_assertion(debug_test_intersection_consistency(mesh_1, mesh_2, intersection_events, v_index_1, v_index_2, h_index_1, h_index_2, vout));
#endif

	// Prepare merged_points for later.
	std::vector<Point> merged_points{};
	merged_points.insert(merged_points.end(), points_1.begin(), points_1.end());
	merged_points.insert(merged_points.end(), points_2.begin(), points_2.end());
	merged_points.insert(merged_points.end(), intersection_points.begin(), intersection_points.end());

	// Debug.
//	debug_print_intersections(vout, intersection_points, intersections_of_edges_1, intersections_of_edges_2,
//							  v_index_1, v_index_2, V1, V2, mesh_1, mesh_2, intersections_data);

	// Intersection collected. Now construct the merged mesh.
	// The points of the merged mesh are prepared in 'merged_points'.
	// Now generate the faces of the merged mesh.
#if USE_DEPRECATED_INTERSECTION_METHOD
	std::vector<std::vector<size_t>> merged_faces = generate_merged_faces_deprecated(
		mesh_1, mesh_2, intersections_of_edges_1, intersections_of_edges_2, intersections_data, V1, V2, V3, v_index_1, h_index_1, v_index_2, h_index_2, vout
	);
#else
	// Account for the vertex-vertex overlaps. The merged mesh should have fewer than V1+V2 verts as a result of v-v overlaps and that needs to be implemented carefully!
	merged_points = generate_merged_points(merged_points, vert_index_mapping_for_merged_mesh);

	std::vector<std::vector<size_t>> merged_faces = generate_merged_faces_from_events(
		mesh_1, mesh_2, vert_index_mapping_for_merged_mesh, intersection_events, intersections_data, vertex_overlap_map_1, vertex_overlap_map_2,
		hedges_through_vertex_1, hedges_through_vertex_2, V1, V2, V3, v_index_1, h_index_1, v_index_2, h_index_2, vout
	);

#endif

	// Debug.
	vout << "About to generate merged mesh from:" << std::endl;
	debug_print_face_types(merged_points, merged_faces, vout);
	//	debug_print_faces_and_points(merged_points, merged_faces, vout);

		// Construct the merged embedding:
	Polyhedron merged_mesh;
	//make_mesh_from_points_and_faces_skipping_orientation(merged_points, merged_faces, merged_mesh);
	make_mesh_from_points_and_faces(merged_points, merged_faces, merged_mesh);

	// Debug.
	vout << "Merged mesh generated successfully!" << std::endl;

	// To return.
	Polyhedron result_merged_embedding{ merged_mesh };

	// Triangulate the merged embedding mesh.
	CGAL::Polygon_mesh_processing::triangulate_faces(merged_mesh);
	std::vector<Point> merged_triangulated_points;
	std::vector<std::vector<size_t>> merged_triangulated_faces;
	get_faces(merged_triangulated_faces, merged_mesh);
	get_points(merged_triangulated_points, merged_mesh);

	// Debug.
	vout << "Triangulated the merged mesh successfully!" << std::endl;
	debug_print_face_types(merged_triangulated_points, merged_triangulated_faces, vout);

#if !USE_DEPRECATED_INTERSECTION_METHOD
	generate_verts_of_containing_face(mesh_1, mesh_2, verts_of_containing_face_1, verts_of_containing_face_2, v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2);
#endif

	size_t n_merged_points = merged_triangulated_points.size();
	CGAL_assertion(n_merged_points == merged_points.size());

	// Impersonate mesh_1.
	vout << "Impersonating mesh_1..." << std::endl;
	std::vector<Point> merged_triangulated_points_impersonating_mesh_1 =
		impersonate(true, V1, V2, V3, n_merged_points, mesh_1, mesh_2, original_mesh_points_1, verts_of_containing_face_1, verts_of_containing_face_2,
			intersections_data, intersection_points, vert_index_mapping_for_merged_mesh, v_index_1, v_index_2);
	Polyhedron mesh_impersonation_1;
	make_mesh_from_points_and_faces(merged_triangulated_points_impersonating_mesh_1, merged_triangulated_faces, mesh_impersonation_1);
	vout << "Mesh_1 impersonated successfully!" << std::endl;

	// Impersonate mesh_2.
	vout << "Impersonating mesh_2..." << std::endl;
	std::vector<Point> merged_triangulated_points_impersonating_mesh_2 =
		impersonate(false, V1, V2, V3, n_merged_points, mesh_1, mesh_2, original_mesh_points_2, verts_of_containing_face_1, verts_of_containing_face_2,
			intersections_data, intersection_points, vert_index_mapping_for_merged_mesh, v_index_1, v_index_2);
	Polyhedron mesh_impersonation_2;
	make_mesh_from_points_and_faces(merged_triangulated_points_impersonating_mesh_2, merged_triangulated_faces, mesh_impersonation_2);
	vout << "Mesh_2 impersonated successfully!" << std::endl;

	return std::make_tuple(mesh_impersonation_1, mesh_impersonation_2, result_merged_embedding);
}