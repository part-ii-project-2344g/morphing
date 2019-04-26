#include "debug_utils.h"

#include "cpp_utils.h"
#include "mesh_merger.h"
using namespace std;

std::string hedge_to_str(size_t index, Polyhedron& mesh, const VInvIndex v_index, bool is_mesh_1, size_t V1) {
	HCI h = get_halfedge(mesh, index);
	size_t va = v_index[h->opposite()->vertex()];
	size_t vb = v_index[h->vertex()];
	if (!is_mesh_1) {
		va += V1;
		vb += V1;
	}
	return std::string{ "m" } + (is_mesh_1 ? "1" : "2") + ".h_" + std::to_string(index) + "(v_" + std::to_string(va) + " -> v_" + std::to_string(vb) + ")";
}

void debug_print_intersections(
	CGAL::Verbose_ostream& vout, 
	const std::vector<Point>& intersection_points,
	std::unordered_map<size_t, std::vector<size_t>> intersections_of_edges_1,
	std::unordered_map<size_t, std::vector<size_t>> intersections_of_edges_2,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2,
	size_t V1,
	size_t V2,
	Polyhedron& mesh_1,
	Polyhedron& mesh_2,
	std::unordered_map<size_t, std::unordered_map<size_t, size_t>> intersections_data
) {
	vout << "# Intersection points: " << intersection_points.size() << endl;
	vout << "Intersection points: "; for (Point p : intersection_points) vout << p << "; "; vout << endl;
	vout << "Intersections for edges of mesh 1: \n";
	for (auto kvp : intersections_of_edges_1) {
		vout << hedge_to_str(kvp.first, mesh_1, v_index_1, true, V1) << ":" << endl;
		for (auto ind : kvp.second) {
			vout << hedge_to_str(ind, mesh_2, v_index_2, false, V1) << ",";
		}
		vout << endl;
	}
	vout << endl;
	vout << "Intersections for edges of mesh 2:\n";
	for (auto kvp : intersections_of_edges_2) {
		vout << hedge_to_str(kvp.first, mesh_2, v_index_2, false, V1) << ":" << endl;
		for (auto ind : kvp.second) {
			vout << hedge_to_str(ind, mesh_1, v_index_1, true, V1) << ",";
		}
		vout << endl;
	}
	vout << endl;
	vout << "New vertices:\n";
	for (auto kvp_outer : intersections_data) {
		size_t ind_1 = kvp_outer.first;
		for (auto kvp_inner : kvp_outer.second) {
			size_t ind_2 = kvp_inner.first;
			size_t ind_3 = kvp_inner.second;
			vout << hedge_to_str(ind_1, mesh_1, v_index_1, true, V1) << "  x  ";
			vout << hedge_to_str(ind_2, mesh_2, v_index_2, false, V1) << " --> ";
			vout << V1 + V2 + ind_3;
			vout << endl;
		}
	}
	vout << endl;
}

void debug_print_hedges(Polyhedron& mesh_1, Polyhedron& mesh_2, 
	const VInvIndex& v_index_1, const HInvIndex& h_index_1, 
	const VInvIndex& v_index_2, const HInvIndex& h_index_2, size_t V1, CGAL::Verbose_ostream& vout) {
	// Print hedges.
	for (int i = 0; i < mesh_1.size_of_halfedges(); i++) {
		vout << hedge_to_str(i, mesh_1, v_index_1, true, V1) << endl;
		vout << "m1.h_" << i << ".next() = m1.h_" << get_next(mesh_1, i, h_index_1) << std::endl;
	}
	for (int i = 0; i < mesh_2.size_of_halfedges(); i++) {
		vout << hedge_to_str(i, mesh_2, v_index_2, false, V1) << endl;
		vout << "m2.h_" << i << ".next() = m2.h_" << get_next(mesh_2, i, h_index_2) << std::endl;
	}
}

void debug_print_intersection_events(size_t debug_counter1, const IntersectionEventStore& intersection_events, CGAL::Verbose_ostream& vout) {
	vout << "FINISHED collect_intersection_events(), investigated " << debug_counter1 << " edges." << endl;
	size_t debug_ee_1 = 0, debug_ev_1 = 0, debug_endse_1 = 0, debug_endsv_1 = 0;
	size_t debug_ee_2 = 0, debug_ev_2 = 0, debug_endse_2 = 0, debug_endsv_2 = 0;
	for (auto kvp : intersection_events.intersection_events_1) {
		for (IntersectionEvent ie : kvp.second) {
			switch (ie.tag) {
			case 0: debug_ee_1++; break;
			case 1: debug_ev_1++; break;
			case 2: debug_endse_1++; break;
			case 3: debug_endsv_1++; break;
			default: break;
			}
		}
	}
	for (auto kvp : intersection_events.intersection_events_2) {
		for (IntersectionEvent ie : kvp.second) {
			switch (ie.tag) {
			case 0: debug_ee_2++; break;
			case 1: debug_ev_2++; break;
			case 2: debug_endse_2++; break;
			case 3: debug_endsv_2++; break;
			default: break;
			}
		}
	}
	vout << "\n-- For mesh_1:" << endl;
	vout << "Edge-edge intersections:   " << debug_ee_1 << endl;
	vout << "Edge-vertex intersections: " << debug_ev_1 << endl;
	vout << "Ends-at-edge events:       " << debug_endse_1 << endl;
	vout << "Ends-at-vert events:       " << debug_endsv_1 << endl;
	vout << "\n-- For mesh_2:" << endl;
	vout << "Edge-edge intersections:   " << debug_ee_2 << endl;
	vout << "Edge-vertex intersections: " << debug_ev_2 << endl;
	vout << "Ends-at-edge events:       " << debug_endse_2 << endl;
	vout << "Ends-at-vert events:       " << debug_endsv_2 << endl;
}

bool debug_test_intersection_consistency_h(Polyhedron& mesh_1, Polyhedron& mesh_2, const IntersectionEventStore& intersection_events,
	bool meshes_swapped,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_1,
	const HInvIndex& h_index_2,
	CGAL::Verbose_ostream& vout
) {
	size_t debug_fails = 0;
	size_t debug_successes = 0;
	bool success = true;
	// We want to make sure that we do not 'get out of a triangle' of mesh_2 without noting an intersection with that triangle.
	// An example would be that we intersect T's edge, and then our next intersection is already outside of T (i.e. we somehow got out of T without intersecting it a second time).
	// We can check that by checking that for every EE event, either the predecessor or successor is an event on the same triangle of mesh_2. 
	// If either the predecessor or successor doesn't exist, we don't check anything.
	for (auto kvp : (meshes_swapped ? intersection_events.intersection_events_2 : intersection_events.intersection_events_1)) {
		size_t h_1 = kvp.first;
		std::vector<IntersectionEvent> events = kvp.second;

		for (size_t i = 1; i < events.size()-1; i++) {
			IntersectionEvent ie = events[i];
			if (ie.tag != 0) continue;

			// First suppose that at 'ie' we are entering the triangle. We will see if the next event still is on this triangle.
			bool success_with_next = false;

			// If the next event is an edge intersection, check if that edge is one of the ones we're looking for.
			IntersectionEvent ie_next = events[i + 1];
			if (ie_next.tag == 0) {

				// Identify the halfedges on the inside of the triangle that we are entering.
				// NOTE: This assumes that as we are entering the triangle, the halfedge stored in the e-e event is on the outside of that triangle.
				std::unordered_set<size_t> tri_hedges{};
				HCI h_2 = get_halfedge(mesh_2, ie.other_edge_index())->opposite();
				while (!set_contains<size_t>(tri_hedges, h_index_2[h_2])) {
					tri_hedges.insert(h_index_2[h_2]);
					h_2 = h_2->next();
				}

				// Check if the next event's edge is in that triangle.
				size_t hedge_of_next_event = ie_next.other_edge_index();
				if (set_contains<size_t>(tri_hedges, hedge_of_next_event)) {
					debug_successes++;
					success_with_next = true;
				}
			}
			// If the next event is a vertex intersection, check if that vertex is in the triangle we have identified.
			else if (ie_next.tag == 1) {

				// Identify the vertices on the triangle that we are entering.
				// NOTE: This assumes that as we are entering the triangle, the halfedge stored in the e-e event is on the outside of that triangle.
				std::unordered_set<size_t> tri_verts{};
				HCI h_2 = get_halfedge(mesh_2, ie.other_edge_index())->opposite();
				while (!set_contains<size_t>(tri_verts, v_index_2[h_2->vertex()])) {
					tri_verts.insert(v_index_2[h_2->vertex()]);
					h_2 = h_2->next();
				}

				size_t vert_of_next_event = ie_next.other_vert_index();
				if (set_contains<size_t>(tri_verts, vert_of_next_event)) {
					debug_successes++;
					success_with_next = true;
				}
			}
			else {
				// The next event is an ENDPOINT event, which is also acceptable.
				// Because this means we are never leaving the triangle.
				debug_successes++;
				success_with_next = true;
			}

			if (success_with_next) 
				continue;

			// If we reach this code we've failed with the next event. Let's try with the previous one.
			// Now suppose we are leaving the triangle in 'ie'. We will see if the previous event also is on this triangle.
			bool success_with_prev = false;

			// If the prev event is an edge intersection, check if that edge is one of the ones we're looking for.
			IntersectionEvent ie_prev = events[i - 1];
			if (ie_prev.tag == 0) {

				// Identify the halfedges on the outside of the triangle that we are exiting.
				// NOTE: This assumes that as we are exiting the triangle, the halfedge stored in the e-e event is on the inside of that triangle.
				std::unordered_set<size_t> tri_hedges{};
				HCI h_2 = get_halfedge(mesh_2, ie.other_edge_index());
				while (!set_contains<size_t>(tri_hedges, h_index_2[h_2->opposite()])) {
					tri_hedges.insert(h_index_2[h_2->opposite()]);
					h_2 = h_2->next();
				}

				// Check if the prev event's edge is on the outside of the triangle that we're leaving.
				size_t hedge_of_prev_event = ie_prev.other_edge_index();
				if (set_contains<size_t>(tri_hedges, hedge_of_prev_event)) {
					debug_successes++;
					success_with_prev = true;
				}
			}
			// If the prev event is a vertex intersection, check if that vertex is in the triangle we have identified.
			else if (events[i - 1].tag == 1) {

				// Identify the vertices on the triangle that we are entering.
				// NOTE: This assumes that as we are entering the triangle, the halfedge stored in the e-e event is on the outside of that triangle.
				// That's true due to how the e-e intersections are collected.
				std::unordered_set<size_t> tri_verts{};
				HCI h_2 = get_halfedge(mesh_2, ie.other_edge_index());
				while (!set_contains<size_t>(tri_verts, v_index_2[h_2->vertex()])) {
					tri_verts.insert(v_index_2[h_2->vertex()]);
					h_2 = h_2->next();
				}

				size_t vert_of_prev_event = ie_prev.other_vert_index();
				if (set_contains<size_t>(tri_verts, vert_of_prev_event)) {
					debug_successes++;
					success_with_prev = true;
				}
			}
			else {
				// The next event is an ENDPOINT event, which is also acceptable.
				// Because this means we are never leaving the triangle.
				debug_successes++;
				success_with_prev = true;
			}

			if (success_with_prev) 
				continue;

			// If we reached this code, we have failed the test.
			debug_fails++;

			// Debug.
			vout << "test_intersection_consistency() failed for mesh_" << (meshes_swapped ? "2" : "1") <<"'s hedge #" << h_1 << ". The event has index " << i << " from (0.." << events.size()-1 << ") on the events list for that hedge, event description follows:" << endl;
			vout << ie.to_string((meshes_swapped ? mesh_2 : mesh_1), (meshes_swapped ? mesh_1 : mesh_2), !meshes_swapped, h_1, (meshes_swapped ? h_index_2 : h_index_1), (meshes_swapped ? h_index_1 : h_index_2)) << endl;

			// Exit.
			success = false;
		}
	}
	vout << "test_intersection_consistency() for mesh_" << (meshes_swapped ? "2" : "1") << " registered " << debug_fails << " fails and " << debug_successes << " successes." << endl;
	return success;
}
bool debug_test_intersection_consistency(Polyhedron& mesh_1, Polyhedron& mesh_2, const IntersectionEventStore& intersection_events,
	const VInvIndex& v_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_1,
	const HInvIndex& h_index_2,
	CGAL::Verbose_ostream& vout
) {
	return debug_test_intersection_consistency_h(mesh_1, mesh_2, intersection_events, false, v_index_1, v_index_2, h_index_1, h_index_2, vout)
		&& debug_test_intersection_consistency_h(mesh_2, mesh_1, intersection_events, true,  v_index_2, v_index_1, h_index_2, h_index_1, vout);
}

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
) {
	vout << "debug_everything() running." << endl;
	vout << "\n\n!SECTION #1/14: Vertices of mesh_1:" << endl;
	auto it = mesh_1.vertices_begin();
	for (size_t i = 0; i < mesh_1.size_of_vertices(); i++, it++) {
		vout << "v_" << i << " = " << it->point() << endl;
	}
	vout << "\n\n!SECTION #2/14: Vertices of mesh_2:" << endl;
	it = mesh_2.vertices_begin();
	for (size_t i = 0; i < mesh_2.size_of_vertices(); i++, it++) {
		vout << "v_" << V1 + i << " = m2.v_" << i << " = " << it->point() << endl;
	}
	vout << "\n\n!SECTION #3/14: Intersection vertices:" << endl;
	for (size_t i = 0; i < V3; i++, it++) {
		vout << "v_" << V1 + V2 + i << " = v_(" << V1 + V2 << "+" << i << ") = " << intersection_points[i] << " = " << endl;
	}
	vout << "\n\n!SECTION #4/14: Halfedges of mesh_1:" << endl;
	for (size_t i = 0; i < mesh_1.size_of_halfedges(); i++) {
		vout << hedge_to_str(i, mesh_1, v_index_1, true, V1);
	}
	vout << "\n\n!SECTION #5/14: Halfedges of mesh_2:" << endl;
	for (size_t i = 0; i < mesh_2.size_of_halfedges(); i++) {
		vout << hedge_to_str(i, mesh_2, v_index_2, false, V1);
	}
	vout << "\n\n!SECTION #6/14: Edge-edge intersections:" << endl;
	for (auto kvp_outer : intersections_data) {
		size_t ind_1 = kvp_outer.first;
		for (auto kvp_inner : kvp_outer.second) {
			size_t ind_2 = kvp_inner.first;
			size_t ind_3 = kvp_inner.second;
			vout << hedge_to_str(ind_1, mesh_1, v_index_1, true, V1) << "  x  ";
			vout << hedge_to_str(ind_2, mesh_2, v_index_2, false, V1) << " --> ";
			vout << V1 + V2 + ind_3;
			vout << endl;
		}
	}
	vout << "\n\n!SECTION #7/14: Edge(1)-vertex(2) intersections:" << endl;
	for (auto kvp : intersection_events.intersection_events_1) {
		size_t ind_1 = kvp.first;
		for (IntersectionEvent ie : kvp.second) {
			if (ie.tag == IntersectionEvent::IET_EV) {
				size_t v2 = ie.other_vert_index();
				vout << hedge_to_str(ind_1, mesh_1, v_index_1, true, V1) << " intersects v_" << V1 + v2 << " of mesh 2." << endl;
			}
		}
	}
	vout << "\n\n!SECTION #8/14: Edge(2)-vertex(1) intersections:" << endl;
	for (auto kvp : intersection_events.intersection_events_2) {
		size_t ind_2 = kvp.first;
		for (IntersectionEvent ie : kvp.second) {
			if (ie.tag == IntersectionEvent::IET_EV) {
				size_t v1 = ie.other_vert_index();
				vout << hedge_to_str(ind_2, mesh_2, v_index_2, false, V1) << " intersects v_" << 0 + v1 << " of mesh 1." << endl;
			}
		}
	}
	vout << "\n\n!SECTION #9/14: Vertex(1)-vertex(2) overlaps:" << endl;
	for (auto kvp : vertex_overlap_map_1) {
		size_t v1 = kvp.first;
		for (size_t v2 : kvp.second) {
			vout << "v_" << v1 << " [m1] overlaps with v_" << V1 + v2 << " [m2]" << endl;
		}
	}
	vout << "\n\n!SECTION #10/14: Intersection events for mesh_1:" << endl;
	for (auto kvp : intersection_events.intersection_events_1) {
		size_t ind_1 = kvp.first;
		vout << "Events for " << hedge_to_str(ind_1, mesh_1, v_index_1, true, V1) << endl;
		for (IntersectionEvent ie : kvp.second) {
			switch (ie.tag) {
			case IntersectionEvent::IET_START:
			{
				vout << "  START event";
			}
			break;
			case IntersectionEvent::IET_EE:
			{
				vout << "  EE with " << hedge_to_str(ie.other_edge_index(), mesh_2, v_index_2, false, V1);
			}
			break;
			case IntersectionEvent::IET_EV:
			{
				vout << "  EV with v_" << V1 + ie.other_vert_index();
			}
			break;
			case IntersectionEvent::IET_ENDS_ON_EDGE:
			{
				vout << "  ENDS_ON_EDGE " << hedge_to_str(ie.other_edge_index(), mesh_2, v_index_2, false, V1);
			}
			break;
			case IntersectionEvent::IET_ENDS_ON_VERT:
			{
				vout << "  ENDS_ON_VERT v_" << V1 + ie.other_vert_index() << " of mesh_2.";
			}
			break;
			default:
				break;
			}
			vout << " at " << ie.pos << "." << endl;
		}
	}
	vout << "\n\n!SECTION #11/14: Intersection events for mesh_2:" << endl;
	for (auto kvp : intersection_events.intersection_events_2) {
		size_t ind_2 = kvp.first;
		vout << "Events for " << hedge_to_str(ind_2, mesh_2, v_index_2, false, V1) << endl;
		for (IntersectionEvent ie : kvp.second) {
			switch (ie.tag) {
			case IntersectionEvent::IET_START:
			{
				vout << "  START event";
			}
			break;
			case IntersectionEvent::IET_EE:
			{
				vout << "  EE with " << hedge_to_str(ie.other_edge_index(), mesh_1, v_index_1, true, V1);
			}
			break;
			case IntersectionEvent::IET_EV:
			{
				vout << "  EV with v_" << 0 + ie.other_vert_index();
			}
			break;
			case IntersectionEvent::IET_ENDS_ON_EDGE:
			{
				vout << "  ENDS_ON_EDGE " << hedge_to_str(ie.other_edge_index(), mesh_1, v_index_1, true, V1);
			}
			break;
			case IntersectionEvent::IET_ENDS_ON_VERT:
			{
				vout << "  ENDS_ON_VERT v_" << 0 + ie.other_vert_index() << " of mesh_1.";
			}
			break;
			default:
				break;
			}
			vout << " at " << ie.pos << "." << endl;
		}
	}
	vout << "\n\n!SECTION #12/14: Halfedges(2) of vertices(1):" << endl;
	for (auto kvp : hedges_through_vertex_1) {
		size_t v1 = kvp.first;
		vout << "Halfedges of mesh_2 passing through / ending at v_" << v1 << ":";
		for (size_t h2 : kvp.second) {
			vout << hedge_to_str(h2, mesh_2, v_index_2, false, V1) << ", \n";
		}
		vout << endl;
	}
	vout << "\n\n!SECTION #13/14: Halfedges(1) of vertices(2):" << endl;
	for (auto kvp : hedges_through_vertex_2) {
		size_t v2 = kvp.first;
		vout << "Halfedges of mesh_1 passing through / ending at v_" << v2 << ":";
		for (size_t h1 : kvp.second) {
			vout << hedge_to_str(h1, mesh_1, v_index_1, true, V1) << ", \n";
		}
		vout << endl;
	}
	vout << "\n\n!SECTION #14/14: Vertex index mapping for the merged mesh:" << endl;
	std::map<size_t, size_t> ordered_mapping{};
	for (auto kvp : vert_index_mapping_for_merged_mesh) {
		size_t v_orig = kvp.first;
		size_t v_new = kvp.second;
		ordered_mapping[v_orig] = v_new;
	}
	for (auto kvp : ordered_mapping) {
		size_t v_orig = kvp.first;
		size_t v_new = kvp.second;
		vout << "(v_" << v_orig << " ; v_" << v_new << ") ";
	}
	vout << "\nEnd of debug_everything()." << endl;
}
