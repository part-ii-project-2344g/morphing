#include "intersections_new.h"

#include "../common/constants.h"
#include "cpp_utils.h"
#include "debug_utils.h"
#include "geometry_utils.h"
#include "mesh_merger.h"

using namespace std;

// Quadratic complexity.
void generate_verts_of_containing_face(
	Polyhedron& mesh_1, Polyhedron& mesh_2,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& out_verts_of_containing_face_1,
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>>& out_verts_of_containing_face_2,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2
) {
	out_verts_of_containing_face_1.clear();
	out_verts_of_containing_face_2.clear();

	for (size_t v = 0; v < mesh_1.size_of_vertices(); v++) {
		size_t f_2 = containing_face_arbitrary(get_vertex(mesh_1, v)->point(), mesh_2, f_index_2);
		FCI facet_2 = get_facet(mesh_2, f_2);

		HCI h_it = facet_2->halfedge();
		HCI h_it_end = facet_2->halfedge();
		do {
			insert_into_map<size_t, size_t>(out_verts_of_containing_face_1, v, v_index_2[h_it->vertex()]);
			h_it = h_it->next();
		} while (h_it != h_it_end);
	}

	for (size_t v = 0; v < mesh_2.size_of_vertices(); v++) {
		size_t f_1 = containing_face_arbitrary(get_vertex(mesh_2, v)->point(), mesh_1, f_index_1);
		FCI facet_1 = get_facet(mesh_1, f_1);

		HCI h_it = facet_1->halfedge();
		HCI h_it_end = facet_1->halfedge();
		do {
			insert_into_map<size_t, size_t>(out_verts_of_containing_face_2, v, v_index_1[h_it->vertex()]);
			h_it = h_it->next();
		} while (h_it != h_it_end);
	}
}

bool check_if_we_intersect_a_vertex(Polyhedron& mesh_1, Polyhedron& mesh_2, HCI h_1, VCI vertex_to_check_2,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	std::unordered_set<size_t> hedges_adjacent_to_vertex_to_check_2,
	IntersectionEventStore& out_out_intersection_events,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 hedges*/>>& out_out_hedges_through_vertex_1,
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 hedges*/>>& out_out_hedges_through_vertex_2,
	std::unordered_map<size_t /*mesh_2 vert*/, std::unordered_set<size_t /*mesh_1 edge*/>>& v2_e1_intersections_found,
	bool& out_done_with_h_1,
	bool& out_we_are_intersecting_a_vertex_on_way_out,
	bool& out_we_have_terminated_at_a_vertex_of_mesh_2,
	bool& out_came_in_through_vertex,
	size_t& out_came_in_through_vertex_index,
	bool& out_last_edges_2_valid,
	std::unordered_set<size_t>& out_last_edges_2,
	FCI& out_f_2,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2,
	CGAL::Verbose_ostream& vout
) {
	Point vertex_to_check_2_point = vertex_to_check_2->point();
	size_t vertex_to_check_2_ind = v_index_2[vertex_to_check_2];

	if (!vert_edge_intersection(vertex_to_check_2_point, h_1)) return false;


	out_we_are_intersecting_a_vertex_on_way_out = true;
	// If this runs, we will be done with this iteration, because later on we only do work if (!we_are_intersecting_a_vertex_on_way_out).

	// There are two possibilities. We might be terminating at a vertex of mesh_2, or just passing through it. Check this.
	if (vertex_overlap_map_1.find(get_end_vert(mesh_1, h_index_1[h_1], v_index_1)) != vertex_overlap_map_1.end()) {
		for (size_t vert_2_that_overlaps_with_endpoint_of_h_1 : vertex_overlap_map_1.at(get_end_vert(mesh_1, h_index_1[h_1], v_index_1))) {
			if (vert_2_that_overlaps_with_endpoint_of_h_1 == vertex_to_check_2_ind) {
				// We are terminating at a vertex 'vertex_to_check_2' of mesh_2. 
				// The IE_ENDS_ON_VERT events are added in a separate process, 
				// so we add no events now and are done with h_1.
				out_we_have_terminated_at_a_vertex_of_mesh_2 = true;
				out_done_with_h_1 = true;
				return true;
			}
		}
	}

	if (!out_we_have_terminated_at_a_vertex_of_mesh_2) {
		// We are intersecting a vertex on our way out.
		// 1) Add IntersectionEvent
		if (!set_contains<size_t>(v2_e1_intersections_found[v_index_2[vertex_to_check_2]], h_index_1[h_1])) {
			insert_into_map<size_t, IntersectionEvent>(out_out_intersection_events.intersection_events_1, h_index_1[h_1], IE_EV{ vertex_to_check_2_point, v_index_2[vertex_to_check_2] });
			insert_into_map<size_t, size_t>(out_out_hedges_through_vertex_2, v_index_2[vertex_to_check_2], h_index_1[h_1]);
			v2_e1_intersections_found[v_index_2[vertex_to_check_2]].insert(h_index_1[h_1]);
		}
		// TODO: create a v2e1_found map analogous to the other one
		// 1b) Add IntersectionEvent for the opposite halfedge.
		if (!set_contains<size_t>(v2_e1_intersections_found[v_index_2[vertex_to_check_2]], h_index_1[h_1->opposite()])) {
			insert_into_map<size_t, IntersectionEvent>(out_out_intersection_events.intersection_events_1, h_index_1[h_1->opposite()], IE_EV{ vertex_to_check_2_point, v_index_2[vertex_to_check_2] });
			insert_into_map<size_t, size_t>(out_out_hedges_through_vertex_2, v_index_2[vertex_to_check_2], h_index_1[h_1->opposite()]);
			v2_e1_intersections_found[v_index_2[vertex_to_check_2]].insert(h_index_1[h_1->opposite()]);
		}

		// 2) Figure out which face we're going into next.
		// Do this by checking the containing face an epsilon along h_1 behind the vertex we just intersected
		Vector epsilon_position = (vertex_to_check_2_point - CGAL::ORIGIN) + STEP_SIZE_ALONG_HEDGE * get_edge_direction_normalized(h_1); // TODO: Get the constant right.
		size_t next_face_2 = containing_face_arbitrary(epsilon_position, mesh_2, f_index_2);

		// Advance. 
		out_came_in_through_vertex = true;
		out_came_in_through_vertex_index = v_index_2[vertex_to_check_2];
		out_last_edges_2_valid = true;
		out_last_edges_2.clear();
		for (size_t h_2 : hedges_adjacent_to_vertex_to_check_2) {
			out_last_edges_2.insert(h_2);
			out_last_edges_2.insert(get_opposite(mesh_2, h_2, h_index_2)); // just to be more sure we get it right, any extra edges from other faces don't break anything if present here
		}
		out_f_2 = get_facet(mesh_2, next_face_2);
	}
	return true;
}

// Ensure that if a hedge h starts in a vertex v that overlaps with a vertex w, then it doesn't have any EE events with hedges containing w.
// Should always do nothing because sanitize_events_around_vv_overlaps should deal with everything.
size_t sanitize_events_around_vv_overlaps_2(Polyhedron& mesh_1, Polyhedron& mesh_2,
	IntersectionEventStore& out_intersection_events,
	std::unordered_map<size_t /*mesh_1 edge*/, std::unordered_map<size_t /*mesh_2 edge*/, size_t /*into intersection_points*/>>& out_intersections_data,
	std::vector<Point>& out_intersection_points,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_2,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 hedges*/>>& out_hedges_through_vertex_1,
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 hedges*/>>& out_hedges_through_vertex_2,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2,
	CGAL::Verbose_ostream& vout
) {
	vout << "\nRunning 'sanitize_events_around_vv_overlaps_2'..." << endl;
	size_t debug_counter = 0;

	for (auto kvp : out_intersection_events.intersection_events_1) {
		size_t h1_ind = kvp.first;
		std::vector<IntersectionEvent>* events = &out_intersection_events.intersection_events_1.at(h1_ind);
		size_t h1_start_vert = get_start_vert(mesh_1, h1_ind, v_index_1);

		if (vertex_overlap_map_1.find(h1_start_vert) != vertex_overlap_map_1.end()) {
			for (size_t v2 : vertex_overlap_map_1.at(h1_start_vert)) {
				for (size_t i = events->size() - 1; ; i--) {
					IntersectionEvent event = (*events)[i];
					if (event.tag == IntersectionEvent::IET_EE) {
						size_t other_h_ind = event.other_edge_index();
						size_t other_h_v1 = get_start_vert(mesh_2, other_h_ind, v_index_2);
						size_t other_h_v2 = get_end_vert(mesh_2, other_h_ind, v_index_2);
						if (other_h_v1 == v2 || other_h_v2 == v2) {
							// We need to remove the EE intersections between m2.h_{other_h_ind} and m1.h_{h1_ind}.
							events->erase(events->begin() + i);
							// Erase the dual:
							size_t h1_ind_opp = get_opposite(mesh_1, h1_ind, h_index_1);
							size_t other_h_ind_opp = get_opposite(mesh_2, other_h_ind, h_index_2);
							bool found_it = false;
							std::vector<IntersectionEvent>* other_events = &out_intersection_events.intersection_events_2.at(other_h_ind_opp);
							for (size_t i = 0; i < other_events->size(); i++) {
								if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
									if ((*other_events)[i].other_edge_index() == h1_ind) {
										// Found the corresponding event!
										found_it = true;
										other_events->erase(other_events->begin() + i);
										break;
									}
								}
							}
							CGAL_assertion(found_it);
							debug_counter++;

							// Attempt to erase the other pair.
							found_it = false;
							other_events = &out_intersection_events.intersection_events_2.at(other_h_ind);
							for (size_t i = 0; i < other_events->size(); i++) {
								if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
									if ((*other_events)[i].other_edge_index() == h1_ind_opp) {
										// Found the corresponding event!
										found_it = true;
										other_events->erase(other_events->begin() + i);
										break;
									}
								}
							}
							if (found_it) {
								found_it = false;
								other_events = &out_intersection_events.intersection_events_1.at(h1_ind_opp);
								for (size_t i = 0; i < other_events->size(); i++) {
									if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
										if ((*other_events)[i].other_edge_index() == other_h_ind_opp) {
											// Found the corresponding event!
											found_it = true;
											other_events->erase(other_events->begin() + i);
											break;
										}
									}
								}
								CGAL_assertion(found_it);
								debug_counter++;
							}

						}
					}
					if (i == 0)
						break;
				}
			}
		}
	}

	// do the same thing with meshes swapped
	for (auto kvp : out_intersection_events.intersection_events_2) {
		size_t h2_ind = kvp.first;
		size_t h2_start_vert = get_start_vert(mesh_2, h2_ind, v_index_2);
		std::vector<IntersectionEvent>* events = &out_intersection_events.intersection_events_2.at(h2_ind);

		if (vertex_overlap_map_2.find(h2_start_vert) != vertex_overlap_map_2.end()) {
			for (size_t v1 : vertex_overlap_map_2.at(h2_start_vert)) {
				for (size_t i = events->size() - 1; ; i--) {
					IntersectionEvent event = (*events)[i];
					if (event.tag == IntersectionEvent::IET_EE) {
						size_t other_h_ind = event.other_edge_index();
						size_t other_h_v1 = get_start_vert(mesh_1, other_h_ind, v_index_1);
						size_t other_h_v2 = get_end_vert(mesh_1, other_h_ind, v_index_1);
						if (other_h_v1 == v1 || other_h_v2 == v1) {
							// We need to remove the EE intersections between m2.h_{other_h_ind} and m1.h_{h1_ind}.
							events->erase(events->begin() + i);
							// Erase the dual:
							size_t h2_ind_opp = get_opposite(mesh_2, h2_ind, h_index_2);
							size_t other_h_ind_opp = get_opposite(mesh_1, other_h_ind, h_index_1);
							bool found_it = false;
							std::vector<IntersectionEvent> other_events = out_intersection_events.intersection_events_1.at(other_h_ind_opp);
							for (size_t i = 0; i < other_events.size(); i++) {
								if (other_events[i].tag == IntersectionEvent::IET_EE) {
									if (other_events[i].other_edge_index() == h2_ind) {
										// Found the corresponding event!
										found_it = true;
										other_events.erase(other_events.begin() + i);
										break;
									}
								}
							}
							CGAL_assertion(found_it);
							debug_counter++;

							// Attempt to erase the other pair.
							found_it = false;
							other_events = out_intersection_events.intersection_events_1.at(other_h_ind);
							for (size_t i = 0; i < other_events.size(); i++) {
								if (other_events[i].tag == IntersectionEvent::IET_EE) {
									if (other_events[i].other_edge_index() == h2_ind_opp) {
										// Found the corresponding event!
										found_it = true;
										other_events.erase(other_events.begin() + i);
										break;
									}
								}
							}
							if (found_it) {
								found_it = false;
								other_events = out_intersection_events.intersection_events_2.at(h2_ind_opp);
								for (size_t i = 0; i < other_events.size(); i++) {
									if (other_events[i].tag == IntersectionEvent::IET_EE) {
										if (other_events[i].other_edge_index() == other_h_ind_opp) {
											// Found the corresponding event!
											found_it = true;
											other_events.erase(other_events.begin() + i);
											break;
										}
									}
								}
								CGAL_assertion(found_it);
								debug_counter++;
							}

						}
					}
					if (i == 0)
						break;
				}
			}
		}
	}
	vout << "'sanitize_events_around_vv_overlaps_2' is done and culled " << debug_counter << " E-E event pairs (EEs nearby edge-start verts that have an overlap)!" << endl;
	return debug_counter;
}

// Ensure that if a hedge h has an ENDS_ON_VERTEX v, then it doesn't have any EE events with hedges h' which have v as one of their endpoints.
//                                                               or ENDS_ON_EDGE events with ...
size_t sanitize_events_around_vv_overlaps(Polyhedron& mesh_1, Polyhedron& mesh_2,
	IntersectionEventStore& out_intersection_events,
	std::unordered_map<size_t /*mesh_1 edge*/, std::unordered_map<size_t /*mesh_2 edge*/, size_t /*into intersection_points*/>>& out_intersections_data,
	std::vector<Point>& out_intersection_points,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_2,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 hedges*/>>& out_hedges_through_vertex_1,
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 hedges*/>>& out_hedges_through_vertex_2,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2,
	CGAL::Verbose_ostream& vout
) {
	vout << "\nRunning 'sanitize_events_around_vv_overlaps'..." << endl;
	size_t debug_counter = 0;
	size_t debug_counter_2 = 0;

	for (auto kvp : out_intersection_events.intersection_events_1) {
		size_t h1_ind = kvp.first;
		std::vector<IntersectionEvent>* events = &out_intersection_events.intersection_events_1.at(h1_ind);
		for (IntersectionEvent ie : *events) {
			// TODO: this code should be fixed, because of 'iterator invalidation' on erase,
			// but it should work just fine assuming there's at most one IE_ENDS_ON_VERT,
			// which is the case in my current executions.
			if (ie.tag == IntersectionEvent::IET_ENDS_ON_VERT) {
				size_t v_to_purge = ie.other_vert_index();
				for (size_t i = events->size() - 1; ; i--) {
					IntersectionEvent ie2 = (*events)[i];
					if (ie2.tag == IntersectionEvent::IET_EE || ie2.tag == IntersectionEvent::IET_ENDS_ON_EDGE) {
						size_t other_h_ind = ie2.other_edge_index();
						size_t other_h_v1 = get_start_vert(mesh_2, other_h_ind, v_index_2);
						size_t other_h_v2 = get_end_vert(mesh_2, other_h_ind, v_index_2);
						if (other_h_v1 == v_to_purge || other_h_v2 == v_to_purge) {
							// Remove the E-E event.
							events->erase(events->begin()+i);

							if (ie2.tag == IntersectionEvent::IET_ENDS_ON_EDGE) {
								// We're done with this one.
								continue;
							}
							// Otherwise, if it was an EE, then look for dual EEs and remove them too...


							// Remove the corresponding E-E event in the other hedge (need to take its opposite).
							std::vector<IntersectionEvent>* other_events = &out_intersection_events.intersection_events_2.at(get_opposite(mesh_2, other_h_ind, h_index_2));
							bool found_it = false;
							for (size_t i = 0; i < other_events->size(); i++) {
								if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
									if ((*other_events)[i].other_edge_index() == h1_ind) {
										// Found the corresponding event!
										found_it = true;
										other_events->erase(other_events->begin() + i);
										break;
									}
								}
							}
							CGAL_assertion(found_it); 
							debug_counter++;
							// If this assertion triggers, it means we failed to find the correponding E-E event.

							// Also remove the dual pair of E-E events if it exists.
							size_t h1_ind_opp = get_opposite(mesh_1, h1_ind, h_index_1);
							// search against 'other_h_ind'
							other_events = &out_intersection_events.intersection_events_2.at(other_h_ind);
							found_it = false;
							for (size_t i = 0; i < other_events->size(); i++) {
								if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
									if ((*other_events)[i].other_edge_index() == h1_ind_opp) {
										// Found the corresponding event!
										found_it = true;
										other_events->erase(other_events->begin() + i);
										break;
									}
								}
							}
							if (found_it) { // Only search for the corresponding event, if we found the original one.
								other_events = &out_intersection_events.intersection_events_1.at(h1_ind_opp);
								found_it = false;
								for (size_t i = 0; i < other_events->size(); i++) {
									if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
										if ((*other_events)[i].other_edge_index() == get_opposite(mesh_2, other_h_ind, h_index_2)) {
											// Found the corresponding event!
											found_it = true;
											other_events->erase(other_events->begin() + i);
											break;
										}
									}
								}
								if (!found_it) {
									cerr << "Didn't find dual EE event of m1.h_" << h1_ind_opp << ". There should be an EE intersection with m2.h_" << other_h_ind << "." << endl;
									CGAL_assertion(false);
								}
								debug_counter++;
								// If this assertion triggers, it means we failed to find the correponding E-E event.
							}
						}
					}

					if (i == 0)
						break;
				}
			}
		}
	}

	// do the same thing with meshes swapped
	for (auto kvp : out_intersection_events.intersection_events_2) {
		size_t h2_ind = kvp.first;
		std::vector<IntersectionEvent> events = kvp.second;
		for (IntersectionEvent ie : events) {
			if (ie.tag == IntersectionEvent::IET_ENDS_ON_VERT) {
				size_t v_to_purge = ie.other_vert_index();
				for (size_t i = events.size() - 1; ; i--) {
					IntersectionEvent ie1 = events[i];
					if (ie1.tag == IntersectionEvent::IET_EE) {
						size_t other_h_ind = ie1.other_edge_index();
						size_t other_h_v1 = get_start_vert(mesh_1, other_h_ind, v_index_1);
						size_t other_h_v2 = get_end_vert(mesh_1, other_h_ind, v_index_1);
						if (other_h_v1 == v_to_purge || other_h_v2 == v_to_purge) {
							// Remove the E-E event.
							events.erase(events.begin() + i);

							// Remove the corresponding E-E event in the other hedge (need to take its opposite).
							std::vector<IntersectionEvent>* other_events = &out_intersection_events.intersection_events_1.at(get_opposite(mesh_1, other_h_ind, h_index_1));
							bool found_it = false;
							for (size_t i = 0; i < other_events->size(); i++) {
								if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
									if ((*other_events)[i].other_edge_index() == h2_ind) {
										// Found the corresponding event!
										found_it = true;
										other_events->erase(other_events->begin() + i);
										break;
									}
								}
							}
							CGAL_assertion(found_it);
							debug_counter++;
							// If this assertion triggers, it means we failed to find the correponding E-E event.

							// Also remove the dual pair of E-E events if it exists.
							size_t h2_ind_opp = get_opposite(mesh_2, h2_ind, h_index_2);
							// search against 'other_h_ind'
							other_events = &out_intersection_events.intersection_events_1.at(other_h_ind);
							found_it = false;
							for (size_t i = 0; i < other_events->size(); i++) {
								if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
									if ((*other_events)[i].other_edge_index() == h2_ind_opp) {
										// Found the corresponding event!
										found_it = true;
										other_events->erase(other_events->begin() + i);
										break;
									}
								}
							}
							if (found_it) { // Only search for the corresponding event, if we found the original one.
								other_events = &out_intersection_events.intersection_events_2.at(h2_ind_opp);
								bool found_it = false;
								for (size_t i = 0; i < other_events->size(); i++) {
									if ((*other_events)[i].tag == IntersectionEvent::IET_EE) {
										if ((*other_events)[i].other_edge_index() == get_opposite(mesh_1, other_h_ind, h_index_1)) {
											// Found the corresponding event!
											found_it = true;
											other_events->erase(other_events->begin() + i);
											break;
										}
									}
								}
								if (!found_it) {
									cerr << "Didn't find dual EE event of m2.h_" << h2_ind_opp << ". There should be an EE intersection with m1.h_" << other_h_ind << "." << endl;
									CGAL_assertion(false);
								}
								debug_counter++;
							}
						}
					}

					if (i == 0)
						break;
				}
			}
		}
	}
	vout << "'sanitize_events_around_vv_overlaps' is done and culled " << debug_counter << " E-E event pairs (EEs nearby edge-end verts that have an overlap)!" << endl;
	return debug_counter;
}

// Ensure that a hedge (v->w) where v overlaps with v' of the other mesh has no EV intersection with v'.
size_t sanitize_ev_events_around_vv_overlaps(Polyhedron& mesh_1, Polyhedron& mesh_2,
	IntersectionEventStore& out_intersection_events,
	std::unordered_map<size_t /*mesh_1 edge*/, std::unordered_map<size_t /*mesh_2 edge*/, size_t /*into intersection_points*/>>& out_intersections_data,
	std::vector<Point>& out_intersection_points,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_2,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 hedges*/>>& out_hedges_through_vertex_1,
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 hedges*/>>& out_hedges_through_vertex_2,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2,
	CGAL::Verbose_ostream& vout
) {
	vout << "\nRunning 'sanitize_ev_events_around_vv_overlaps'..." << endl;
	size_t debug_counter = 0;

	for (auto kvp : out_intersection_events.intersection_events_1) {
		size_t h1_ind = kvp.first;
		size_t start_vert_h1 = get_start_vert(mesh_1, h1_ind, v_index_1);
		std::vector<IntersectionEvent>* events = &out_intersection_events.intersection_events_1.at(h1_ind);
		for (size_t i = events->size() - 1; ; i--) {
			IntersectionEvent ie = (*events)[i];
			if (ie.tag == IntersectionEvent::IET_EV && vertex_overlap_map_1.find(start_vert_h1) != vertex_overlap_map_1.end()) {
				size_t intersected_v2 = ie.other_vert_index();
				for (size_t v2 : vertex_overlap_map_1.at(start_vert_h1)) {
					if (intersected_v2 == v2) {
						// We need to remove this EV event.
						events->erase(events->begin() + i);
						// Erase the corresponding 'hedges_through_vertex' entry.
						bool successfully_erased_from_hedges_through_vertex = false;
						if (out_hedges_through_vertex_2.find(v2) != out_hedges_through_vertex_2.end()) {
							std::vector<size_t>& hedges_through_vertex_list = out_hedges_through_vertex_2.at(v2);
							for (size_t i = hedges_through_vertex_list.size() - 1; ; i--) {
								if (hedges_through_vertex_list[i] == h1_ind) {
									hedges_through_vertex_list.erase(hedges_through_vertex_list.begin() + i);
									successfully_erased_from_hedges_through_vertex = true;
								}
								if (i == 0) break;
							}
						}
						vout << "successfully_erased_from_hedges_through_vertex = " << successfully_erased_from_hedges_through_vertex << std::endl;

						debug_counter++;
					}
				}
			}
			if (i == 0) break;
		}
		// handle end_vertex
		size_t end_vert_h1 = get_end_vert(mesh_1, h1_ind, v_index_1);
		events = &out_intersection_events.intersection_events_1.at(h1_ind);
		for (size_t i = events->size() - 1; ; i--) {
			IntersectionEvent ie = (*events)[i];
			if (ie.tag == IntersectionEvent::IET_EV && vertex_overlap_map_1.find(end_vert_h1) != vertex_overlap_map_1.end()) {
				size_t intersected_v2 = ie.other_vert_index();
				for (size_t v2 : vertex_overlap_map_1.at(end_vert_h1)) {
					if (intersected_v2 == v2) {
						// We need to remove this EV event.
						events->erase(events->begin() + i);
						// Erase the corresponding 'hedges_through_vertex' entry.
						bool successfully_erased_from_hedges_through_vertex = false;
						if (out_hedges_through_vertex_2.find(v2) != out_hedges_through_vertex_2.end()) {
							std::vector<size_t>& hedges_through_vertex_list = out_hedges_through_vertex_2.at(v2);
							for (size_t i = hedges_through_vertex_list.size() - 1; ; i--) {
								if (hedges_through_vertex_list[i] == h1_ind) {
									hedges_through_vertex_list.erase(hedges_through_vertex_list.begin() + i);
									successfully_erased_from_hedges_through_vertex = true;
								}
								if (i == 0) break;
							}
						}
						vout << "successfully_erased_from_hedges_through_vertex = " << successfully_erased_from_hedges_through_vertex << std::endl;

						debug_counter++;
					}
				}
			}
			if (i == 0) break;
		}
	}
	// Repeat for the other mesh.
	for (auto kvp : out_intersection_events.intersection_events_2) {
		size_t h2_ind = kvp.first;
		size_t start_vert_h2 = get_start_vert(mesh_2, h2_ind, v_index_2);
		std::vector<IntersectionEvent>* events = &out_intersection_events.intersection_events_2.at(h2_ind);
		for (size_t i = events->size() - 1; ; i--) {
			IntersectionEvent ie = (*events)[i];
			if (ie.tag == IntersectionEvent::IET_EV && vertex_overlap_map_2.find(start_vert_h2) != vertex_overlap_map_2.end()) {
				size_t intersected_v1 = ie.other_vert_index();
				for (size_t v1 : vertex_overlap_map_2.at(start_vert_h2)) {
					if (intersected_v1 == v1) {
						// We need to remove this EV event.
						events->erase(events->begin() + i);
						// Erase the corresponding 'hedges_through_vertex' entry.
						bool successfully_erased_from_hedges_through_vertex = false;
						if (out_hedges_through_vertex_1.find(v1) != out_hedges_through_vertex_1.end()) {
							std::vector<size_t>& hedges_through_vertex_list = out_hedges_through_vertex_1.at(v1);
							for (size_t i = hedges_through_vertex_list.size() - 1; ; i--) {
								if (hedges_through_vertex_list[i] == h2_ind) {
									hedges_through_vertex_list.erase(hedges_through_vertex_list.begin() + i);
									successfully_erased_from_hedges_through_vertex = true;
								}
								if (i == 0) break;
							}
						}
						vout << "successfully_erased_from_hedges_through_vertex = " << successfully_erased_from_hedges_through_vertex << std::endl;

						debug_counter++;
					}
				}
			}
			if (i == 0) break;
		}
		// handle end_vertex
		size_t end_vert_h2 = get_end_vert(mesh_2, h2_ind, v_index_2);
		events = &out_intersection_events.intersection_events_2.at(h2_ind);
		for (size_t i = events->size() - 1; ; i--) {
			IntersectionEvent ie = (*events)[i];
			if (ie.tag == IntersectionEvent::IET_EV && vertex_overlap_map_2.find(end_vert_h2) != vertex_overlap_map_2.end()) {
				size_t intersected_v1 = ie.other_vert_index();
				for (size_t v1 : vertex_overlap_map_2.at(end_vert_h2)) {
					if (intersected_v1 == v1) {
						// We need to remove this EV event.
						events->erase(events->begin() + i);
						// Erase the corresponding 'hedges_through_vertex' entry.
						bool successfully_erased_from_hedges_through_vertex = false;
						if (out_hedges_through_vertex_1.find(v1) != out_hedges_through_vertex_1.end()) {
							std::vector<size_t>& hedges_through_vertex_list = out_hedges_through_vertex_1.at(v1);
							for (size_t i = hedges_through_vertex_list.size() - 1; ; i--) {
								if (hedges_through_vertex_list[i] == h2_ind) {
									hedges_through_vertex_list.erase(hedges_through_vertex_list.begin() + i);
									successfully_erased_from_hedges_through_vertex = true;
								}
								if (i == 0) break;
							}
						}
						vout << "successfully_erased_from_hedges_through_vertex = " << successfully_erased_from_hedges_through_vertex << std::endl;
						
						debug_counter++;
					}
				}
			}
			if (i == 0) break;
		}
	}
	vout << "'sanitize_ev_events_around_vv_overlaps' is done and culled " << debug_counter << " EV events." << endl;
	return debug_counter;
}

void sort_intersection_events_helper(Polyhedron& mesh, const bool doing_mesh_1, IntersectionEventStore& intersection_events, CGAL::Verbose_ostream& vout) {
	for (auto kvp : (doing_mesh_1 ? intersection_events.intersection_events_1 : intersection_events.intersection_events_2)) {
		size_t h = kvp.first;
		std::vector<IntersectionEvent> events = kvp.second;
		Point start_point_of_h = get_halfedge(mesh, h)->opposite()->vertex()->point();
		Point end_point_of_h = get_halfedge(mesh, h)->vertex()->point();

		std::vector<Vector> event_vecs{};
		for (IntersectionEvent event : events) {
			Vector event_pos = event.pos;
			event_vecs.push_back(event_pos);
		}

		// Choose which dimension we will look at to sort the intersection points.
		// Pick the one along which the edge h is the longest.
		int xyz = 0;
		double max_diff = max(max(abs(start_point_of_h.x() - end_point_of_h.x()), abs(start_point_of_h.y() - end_point_of_h.y())), abs(start_point_of_h.z() - end_point_of_h.z()));
		if (max_diff == abs(start_point_of_h.x() - end_point_of_h.x())) xyz = 0;
		else if (max_diff == abs(start_point_of_h.y() - end_point_of_h.y())) xyz = 1;
		else xyz = 2;

		std::vector<std::pair<double, size_t>> sorter{};
		for (size_t i = 0; i < event_vecs.size(); i++) {
			// To sort points on an arc, we assume that the arc is shorter than half the circumference
			// and equivalently sort projections onto the line going through the endpoints of the arc
			Vector projected_event_vec = line_projection(start_point_of_h - CGAL::ORIGIN, end_point_of_h - CGAL::ORIGIN, event_vecs[i]);
			double event_coord = projected_event_vec.x();
			if (xyz == 1) event_coord = projected_event_vec.y();
			else if (xyz == 2) event_coord = projected_event_vec.z();
			sorter.push_back(make_pair(event_coord, i));
		}
		sort(sorter.begin(), sorter.end());

		double start_coord = start_point_of_h.x();
		if (xyz == 1) start_coord = start_point_of_h.y();
		else if (xyz == 2) start_coord = start_point_of_h.z();
		double end_coord = end_point_of_h.x();
		if (xyz == 1) end_coord = end_point_of_h.y();
		else if (xyz == 2) end_coord = end_point_of_h.z();

		if (start_coord > end_coord) {
			std::reverse(sorter.begin(), sorter.end());
		}

		// Insert intersections in the right order this time.
		std::vector<IntersectionEvent> correct_order_events_for_h{};
		for (auto pair : sorter) {
			size_t ind = pair.second;
			correct_order_events_for_h.push_back(events[ind]);
		}

		// Overwrite intersections with the correctly ordered list.
		CGAL_assertion(correct_order_events_for_h.size() == events.size());
		(doing_mesh_1 ? intersection_events.intersection_events_1 : intersection_events.intersection_events_2)[h] = correct_order_events_for_h;
	}
}
void sort_intersection_events(Polyhedron& mesh_1, Polyhedron& mesh_2, IntersectionEventStore& intersection_events, CGAL::Verbose_ostream& vout) {
	sort_intersection_events_helper(mesh_1, true, intersection_events, vout);
	sort_intersection_events_helper(mesh_2, false, intersection_events, vout);
}


// Collect all IntersectionEvents, taking into account vertex-vertex overlap but not accounting for any vertex-edge intersections. Those will be found later.
void collect_intersection_events(Polyhedron& mesh_1, Polyhedron& mesh_2,
	IntersectionEventStore& out_intersection_events,
	std::unordered_map<size_t /*mesh_1 edge*/, std::unordered_map<size_t /*mesh_2 edge*/, size_t /*into intersection_points*/>>& out_intersections_data,
	std::vector<Point>& out_intersection_points,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_2,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 hedges*/>>& out_hedges_through_vertex_1,
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 hedges*/>>& out_hedges_through_vertex_2,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2,
	CGAL::Verbose_ostream& vout
) {
	// Initialize the event lists with IE_START event for every halfedge of each mesh.
	for (size_t h_1 = 0; h_1 < mesh_1.size_of_halfedges(); h_1++) {
		insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_1, h_1, IE_START{ get_halfedge(mesh_1, h_1)->opposite()->vertex()->point() });
	}
	for (size_t h_2 = 0; h_2 < mesh_2.size_of_halfedges(); h_2++) {
		insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_2, h_2, IE_START{ get_halfedge(mesh_2, h_2)->opposite()->vertex()->point() });
	}

	// Initialize the algorithm.
	Halfedge h_1 = mesh_1.halfedges_begin(); // an arbitrary halfedge of mesh_1
	Vertex v_1 = h_1->opposite()->vertex(); // the origin of h_1
	FCI f_2 = get_facet(mesh_2, containing_face(v_1->point(), get_edge_direction_normalized(h_1), mesh_2, f_index_2));
	std::deque<Halfedge> halfedge_stack{};
	halfedge_stack.push_front(h_1);
	std::deque<FCI> face_stack{};
	face_stack.push_front(f_2);
	std::unordered_set<size_t> used{};

	// This remembers which VE intersections of a vertex from mesh_1 with an edge from mesh_2 have already been found, for deduplication purposes.
	std::unordered_map<size_t /*mesh_1 vert*/, std::unordered_set<size_t /*mesh_2 edge*/>> v1_e2_intersections_found{};

	// This remembers which VE intersections of a vertex from mesh_2 with an edge from mesh_1 have already been found, for deduplication purposes.
	std::unordered_map<size_t /*mesh_2 vert*/, std::unordered_set<size_t /*mesh_1 edge*/>> v2_e1_intersections_found{};

	// bookkeeping
	used.insert(h_index_1[h_1]);
	used.insert(h_index_1[h_1->opposite()]);

	// Debug.
	size_t debug_counter1 = 0;

	while (!halfedge_stack.empty()) {
		// pop halfedge of mesh_1 and face of mesh_2 containing the origin of that halfedges from stack
		h_1 = stack_pop<Halfedge>(halfedge_stack);
		f_2 = stack_pop<FCI>(face_stack);

		std::unordered_set<size_t> last_edges_2{}; // the set of edges that we should not test for intersections with h_1
		bool last_edges_2_valid = false;
		bool came_in_through_vertex = false;
		size_t came_in_through_vertex_index = -1;

		// debug
		debug_counter1++;
		auto debug_popped = f_index_2[f_2];
		auto debug_correct = containing_face(h_1->opposite()->vertex()->point(), get_edge_direction_normalized(h_1), mesh_2, f_index_2);

		// Temporary!
		f_2 = get_facet(mesh_2, containing_face(h_1->opposite()->vertex()->point(), get_edge_direction_normalized(h_1), mesh_2, f_index_2));
		//CGAL_assertion(debug_popped == debug_correct); // TODO: Get this assertion to pass without the hack!
		//vout << "Popped halfedge " << h_index_1[h_1] << endl;


		// Check whether h_1 starts at a vertex of mesh_2.
		// Do this by checking if h_1's start vertex overlaps with a vertex of mesh_2.
		if (vertex_overlap_map_1.find(get_start_vert(mesh_1, h_index_1[h_1], v_index_1)) != vertex_overlap_map_1.end()) {
			if (!vertex_overlap_map_1.at(get_start_vert(mesh_1, h_index_1[h_1], v_index_1)).empty()) {
				CGAL_assertion(vertex_overlap_map_1.at(get_start_vert(mesh_1, h_index_1[h_1], v_index_1)).size() == 1);
				came_in_through_vertex = true;
				came_in_through_vertex_index = vertex_overlap_map_1.at(get_start_vert(mesh_1, h_index_1[h_1], v_index_1))[0];
			}
		}

		// Check whether h_1 starts at an edge of mesh_2.
		// Do this by checking if h_1's start vertex intersects any edge of f_2.
		// Only do this if we didn't find that h_1 starts at a vertex of mesh_2. (hence 'if (!came_in_through_vertex)' at the beginning)
		// Together with this, we also later check whether h_2 ends at an edge of mesh_2.
		// By doing both these things, we make sure we have found all the vertex-edge intersections for edges of mesh_2!
		if (!came_in_through_vertex) {
			HCI edge_of_f_2 = f_2->halfedge()->next();
			HCI end = f_2->halfedge();
			while (true) {
				Point intersection_point;
				if (vert_edge_intersection(h_1->opposite()->vertex()->point(), edge_of_f_2)) {
					// h_1 starts at an edge of mesh_2!

					// 1) Add an EV IntersectionEvent for mesh_2.
					// First check if this EV hasn't already been added - this is what we use 'v1_e2_intersections_found' for.
					size_t start_vertex_of_h_1_index = v_index_1[h_1->opposite()->vertex()];
					Point start_vertex_of_h_1_point = h_1->opposite()->vertex()->point();
					size_t index_of_edge_of_f_2 = h_index_2[edge_of_f_2];
					if (v1_e2_intersections_found.find(start_vertex_of_h_1_index) == v1_e2_intersections_found.end()) {
						v1_e2_intersections_found[start_vertex_of_h_1_index] = std::unordered_set<size_t>{};
					} 
					if (!set_contains<size_t>(v1_e2_intersections_found[start_vertex_of_h_1_index], index_of_edge_of_f_2))
					{
						insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_2, index_of_edge_of_f_2, IE_EV{ start_vertex_of_h_1_point, start_vertex_of_h_1_index });
						insert_into_map<size_t, size_t>(out_hedges_through_vertex_1, start_vertex_of_h_1_index, index_of_edge_of_f_2);
						v1_e2_intersections_found[start_vertex_of_h_1_index].insert(index_of_edge_of_f_2);
					}
					// Also add the EV event for the opposite halfedge!
					size_t index_of_opposite_edge_of_f_2 = get_opposite(mesh_2, index_of_edge_of_f_2, h_index_2);
					if (!set_contains<size_t>(v1_e2_intersections_found[start_vertex_of_h_1_index], index_of_opposite_edge_of_f_2))
					{
						insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_2, index_of_opposite_edge_of_f_2, IE_EV{ start_vertex_of_h_1_point, start_vertex_of_h_1_index });
						insert_into_map<size_t, size_t>(out_hedges_through_vertex_1, start_vertex_of_h_1_index, index_of_opposite_edge_of_f_2);
						v1_e2_intersections_found[start_vertex_of_h_1_index].insert(index_of_opposite_edge_of_f_2);
					}

					// 2) Figure out which face we're going into next.
					// Do this by checking the containing face an epsilon along h_1 behind the vertex we just intersected
					Vector epsilon_position = (start_vertex_of_h_1_point - CGAL::ORIGIN) + STEP_SIZE_ALONG_HEDGE * get_edge_direction_normalized(h_1); // TODO: Get the constant right.
					size_t next_face_2 = containing_face_arbitrary(epsilon_position, mesh_2, f_index_2);

					// Advance.
					came_in_through_vertex = false;
					last_edges_2_valid = true;
					last_edges_2 = std::unordered_set<size_t>{ h_index_2[edge_of_f_2], h_index_2[edge_of_f_2->opposite()] }; // We don't have to check which way we went, just put both of them in. The opposite one being there won't change anything.
					f_2 = get_facet(mesh_2, next_face_2);
					break;
				}
				else { // h_1 did not intersect around_face_2.
					if (edge_of_f_2 == end) {
						break;
					}
					edge_of_f_2 = edge_of_f_2->next(); // There's still another edge of f_2 left to try, do it.
				}
			}
		}


		// Go along h_1 and find its intersections with edges of mesh_2. (Main step.)
		bool done_with_h_1 = false;
		while (!done_with_h_1) {

			// What's going on:
			// We are following h_1 and are currently inside the triangle described by the around_face_2 circulator.
			// We have found out that h_1 intersects a side of that triangle, namely the one currently pointed to by around_face_2.
			// We need to check whether h_1 actually goes through a vertex of that triangle, or just cleanly through the side.
			// In the latter case, we get the next face of mesh_2 by flipping across the edge 'around_face_2'.
			// In the former case, we need to figure out which triangle of mesh_2 we're going into through the vertex.
			// We also need to consider the corner case of this corner case, in which we do not go into a triangle but straight into an edge of mesh_2.

			// TODO: Check if the vertex we're intersecting is not identified with our own endpoint!

			// Check whether we're intersecting a vertex of f_2 on our way out.
			bool we_are_intersecting_a_vertex_on_way_out = false;
			bool we_have_terminated_at_a_vertex_of_mesh_2 = false;
			if (came_in_through_vertex) {
				// Check whether we're intersecting a vertex of f_2 on our way out.
				// Check the two vertices that are not the vertex we came in through.

				// Identify the two vertices to check.
				HCI edge_of_f_2 = f_2->halfedge();
				size_t debug_counter = 0;
				while (v_index_2[edge_of_f_2->vertex()] == came_in_through_vertex_index || v_index_2[edge_of_f_2->opposite()->vertex()] == came_in_through_vertex_index) {
					edge_of_f_2 = edge_of_f_2->next();
					if(debug_counter++ > 100) { // we must have entered an infinite loop due to a 2-edge face I guess? 
						CGAL_assertion(false);
					}
				}
				VCI vertex_to_check_2_a = edge_of_f_2->vertex();
				VCI vertex_to_check_2_b = edge_of_f_2->opposite()->vertex();

				// Check one of the candidates, i.e. 'vertex_to_check_2_a'.
				bool vertex_to_check_2_a_was_intersected = check_if_we_intersect_a_vertex(
					mesh_1, mesh_2, h_1,
					vertex_to_check_2_a, // this arg changes between calls
					vertex_overlap_map_1,
					std::unordered_set<size_t>{ h_index_2[edge_of_f_2->opposite()], h_index_2[edge_of_f_2->next()->opposite()] },  // this arg changes between calls
					out_intersection_events, out_hedges_through_vertex_1, out_hedges_through_vertex_2, 
					v2_e1_intersections_found,
					done_with_h_1,
					we_are_intersecting_a_vertex_on_way_out, we_have_terminated_at_a_vertex_of_mesh_2,
					came_in_through_vertex, came_in_through_vertex_index, last_edges_2_valid, last_edges_2, f_2,
					v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout
				);

				// Check the other one if the previous one was not intersected.
				if (!vertex_to_check_2_a_was_intersected) {
					bool vertex_to_check_2_b_was_intersected = check_if_we_intersect_a_vertex( // return value unused
						mesh_1, mesh_2, h_1,
						vertex_to_check_2_b, // this arg changes between calls
						vertex_overlap_map_1,
						std::unordered_set<size_t>{ h_index_2[edge_of_f_2->opposite()], h_index_2[edge_of_f_2->opposite()->next()] }, // this arg changes between calls
						out_intersection_events, out_hedges_through_vertex_1, out_hedges_through_vertex_2, 
						v2_e1_intersections_found,
						done_with_h_1,
						we_are_intersecting_a_vertex_on_way_out, we_have_terminated_at_a_vertex_of_mesh_2,
						came_in_through_vertex, came_in_through_vertex_index, last_edges_2_valid, last_edges_2, f_2,
						v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout
					);
				}
			}
			else if (!came_in_through_vertex) { // We came in through an edge (or just started).
				if (last_edges_2.empty()) { // if it's empty it means we just started and so we should check all three vertices
					HCI edge_of_f_2 = f_2->halfedge();
					VCI vertex_to_check_2_a = edge_of_f_2->vertex();
					VCI vertex_to_check_2_b = edge_of_f_2->next()->vertex();
					VCI vertex_to_check_2_c = edge_of_f_2->next()->vertex();

					bool vertex_to_check_2_a_was_intersected = check_if_we_intersect_a_vertex(
						mesh_1, mesh_2, h_1,
						vertex_to_check_2_a, // this arg changes between calls
						vertex_overlap_map_1,
						std::unordered_set<size_t>{ h_index_2[edge_of_f_2->opposite()], h_index_2[edge_of_f_2->next()->opposite()] },  // this arg changes between calls
						out_intersection_events, out_hedges_through_vertex_1, out_hedges_through_vertex_2,
						v2_e1_intersections_found,
						done_with_h_1,
						we_are_intersecting_a_vertex_on_way_out, we_have_terminated_at_a_vertex_of_mesh_2,
						came_in_through_vertex, came_in_through_vertex_index, last_edges_2_valid, last_edges_2, f_2,
						v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout
					);
					// Check the other one if the previous one was not intersected.
					if (!vertex_to_check_2_a_was_intersected) {
						bool vertex_to_check_2_c_was_intersected = check_if_we_intersect_a_vertex(
							mesh_1, mesh_2, h_1,
							vertex_to_check_2_c, // this arg changes between calls
							vertex_overlap_map_1,
							std::unordered_set<size_t>{ h_index_2[edge_of_f_2->opposite()], h_index_2[edge_of_f_2->opposite()->next()] }, // this arg changes between calls
							out_intersection_events, out_hedges_through_vertex_1, out_hedges_through_vertex_2,
							v2_e1_intersections_found,
							done_with_h_1,
							we_are_intersecting_a_vertex_on_way_out, we_have_terminated_at_a_vertex_of_mesh_2,
							came_in_through_vertex, came_in_through_vertex_index, last_edges_2_valid, last_edges_2, f_2,
							v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout
						);
						
						if (!vertex_to_check_2_c_was_intersected) {
							bool vertex_to_check_2_b_was_intersected = check_if_we_intersect_a_vertex( // return value unused
								mesh_1, mesh_2, h_1,
								vertex_to_check_2_b, // this arg changes between calls
								vertex_overlap_map_1,
								std::unordered_set<size_t>{ h_index_2[edge_of_f_2->next()->opposite()], h_index_2[edge_of_f_2->next()->next()->opposite()] }, // this arg changes between calls
								out_intersection_events, out_hedges_through_vertex_1, out_hedges_through_vertex_2,
								v2_e1_intersections_found,
								done_with_h_1,
								we_are_intersecting_a_vertex_on_way_out, we_have_terminated_at_a_vertex_of_mesh_2,
								came_in_through_vertex, came_in_through_vertex_index, last_edges_2_valid, last_edges_2, f_2,
								v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout
							);
						}
					}
				}
				else { // We really did come in through an edge, so 
				// Check whether we're intersecting a vertex of f_2 on our way out.
				// Only check the vertex opposite the edge that we came in through.
					HCI edge_we_came_in_through = get_halfedge(mesh_2, *(last_edges_2.begin()));
					VCI vertex_to_check_2 = edge_we_came_in_through->next()->vertex();

					bool vertex_to_check_2_was_intersected = check_if_we_intersect_a_vertex( // return value unused
						mesh_1, mesh_2, h_1,
						vertex_to_check_2, // this arg changes between calls
						vertex_overlap_map_1,
						std::unordered_set<size_t>{ h_index_2[edge_we_came_in_through->next()->opposite()], h_index_2[edge_we_came_in_through->next()->next()->opposite()] }, // this arg changes between calls
						out_intersection_events, out_hedges_through_vertex_1, out_hedges_through_vertex_2,
						v2_e1_intersections_found,
						done_with_h_1,
						we_are_intersecting_a_vertex_on_way_out, we_have_terminated_at_a_vertex_of_mesh_2,
						came_in_through_vertex, came_in_through_vertex_index, last_edges_2_valid, last_edges_2, f_2,
						v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout
					);
				}
			}

			// Almost the common case - we're either:
			// 1) terminating inside f_2, or
			// 2) cleanly exiting f_2 through exactly one of its edges, which we will find.
			// 3) [!] terminating at an edge of f_2. - in this case we do not want to store an EE intersection, but an EV intersection for mesh_2 and an ENDS_AT_EDGE event for mesh_1.
			if (!we_have_terminated_at_a_vertex_of_mesh_2 && !we_are_intersecting_a_vertex_on_way_out) {

				bool case_3_took_place = false;

				// Check if case 3) is taking place.
				// Check whether h_1 ends at an edge of mesh_2.
				// Do this by checking if h_1's end vertex intersects any edge of f_2.
				// Together with this, we also have earlier checked whether h_2 starts at an edge of mesh_2.
				// By doing both these things, we make sure we have found all the vertex-edge intersections for edges of mesh_2
				HCI edge_of_f_2 = f_2->halfedge()->next();
				HCI end = f_2->halfedge();
				while (true) {
					Point intersection_point;
					if (vert_edge_intersection(h_1->vertex()->point(), edge_of_f_2)) {
						// h_1 ends at an edge of mesh_2!

						// 1) Add an EV IntersectionEvent for mesh_2.
						// First check if this EV hasn't already been added - this is what we use 'v1_e2_intersections_found' for.
						size_t end_vertex_of_h_1_index = v_index_1[h_1->vertex()];
						Point end_vertex_of_h_1_point = h_1->vertex()->point();
						size_t index_of_edge_of_f_2 = h_index_2[edge_of_f_2];
						if (v1_e2_intersections_found.find(end_vertex_of_h_1_index) == v1_e2_intersections_found.end()) {
							v1_e2_intersections_found[end_vertex_of_h_1_index] = std::unordered_set<size_t>{};
						}
						if (!set_contains<size_t>(v1_e2_intersections_found[end_vertex_of_h_1_index], index_of_edge_of_f_2))
						{
							insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_2, index_of_edge_of_f_2, IE_EV{ end_vertex_of_h_1_point, end_vertex_of_h_1_index });
							insert_into_map<size_t, size_t>(out_hedges_through_vertex_1, end_vertex_of_h_1_index, index_of_edge_of_f_2);
							v1_e2_intersections_found[end_vertex_of_h_1_index].insert(index_of_edge_of_f_2);
						}
						// Also add the EV event for the opposite halfedge!
						size_t index_of_opposite_edge_of_f_2 = get_opposite(mesh_2, index_of_edge_of_f_2, h_index_2);
						if (!set_contains<size_t>(v1_e2_intersections_found[end_vertex_of_h_1_index], index_of_opposite_edge_of_f_2))
						{
							insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_2, index_of_opposite_edge_of_f_2, IE_EV{ end_vertex_of_h_1_point, end_vertex_of_h_1_index });
							insert_into_map<size_t, size_t>(out_hedges_through_vertex_1, end_vertex_of_h_1_index, index_of_opposite_edge_of_f_2);
							v1_e2_intersections_found[end_vertex_of_h_1_index].insert(index_of_opposite_edge_of_f_2);
						}

						// 2) Add an END IntersectionEvent for mesh_1.
						IE_ENDS_ON_EDGE intersection_event_mesh_1{ end_vertex_of_h_1_point, index_of_edge_of_f_2 };
						insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_1, h_index_1[h_1], intersection_event_mesh_1);
						
						// 3) Exit.
						case_3_took_place = true;
						done_with_h_1 = true;
						break;
					}
					else { // h_1 did not intersect around_face_2.
						if (edge_of_f_2 == end) {
							break;
						}
						edge_of_f_2 = edge_of_f_2->next(); // There's still another edge of f_2 left to try, do it.
					}
				}

				// Set the last_edges_2 appropriately if we entered through vertex:
				if (came_in_through_vertex) {
					HCI around_face_2 = f_2->halfedge();
					HCI around_face_2_end = f_2->halfedge();
					while (true) {
						if (get_start_vert(mesh_2, h_index_2[around_face_2], v_index_2) == came_in_through_vertex_index ||
							get_end_vert(mesh_2, h_index_2[around_face_2], v_index_2) == came_in_through_vertex_index) {
							last_edges_2.insert(h_index_2[around_face_2]);
							last_edges_2_valid = true;
						}

						around_face_2 = around_face_2->next();
						if (around_face_2 == around_face_2_end) { break; }
					}
				}

				// If case 3) didn't take place, i.e. h_1 does not end on an edge of mesh_2 (or a vertex of mesh_2, as checked before).
				// The common case - we're either:
				// 1) terminating inside f_2, or
				// 2) cleanly exiting f_2 through exactly one of its edges, which we will find.
				if (!case_3_took_place) {
					HCI around_face_2 = f_2->halfedge()->next();
					HCI around_face_2_end = f_2->halfedge();
					while (true) { // check all the edges of f_2 except for the one(s) we came in through
								   // explanation: there might be 'two' edges we came in through if we came in through a vertex

						Point intersection_point;
						if ((!last_edges_2_valid || !set_contains<size_t>(last_edges_2, h_index_2[around_face_2])) && intersect(around_face_2, h_1, intersection_point)) {

							// Register the intersection.
							insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_1, h_index_1[h_1], IE_EE{ intersection_point, h_index_2[around_face_2] });
							insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_1, h_index_1[h_1->opposite()], IE_EE{ intersection_point, h_index_2[around_face_2->opposite()] });
							// We insert the opposite of h_1 here to keep the invariant that when an edge of one mesh enters a triangle of the other mesh, it stores an intersection with an outside halfedge of that triangle.
							insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_2, h_index_2[around_face_2], IE_EE{ intersection_point, h_index_1[h_1->opposite()] });
							insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_2, h_index_2[around_face_2->opposite()], IE_EE{ intersection_point, h_index_1[h_1] });

							// Officially create an intersection point.
							out_intersection_points.push_back(intersection_point);
							size_t intersection_point_index = out_intersection_points.size() - 1;
							insert_into_map<size_t, size_t, size_t>(out_intersections_data, h_index_1[h_1], h_index_2[around_face_2], intersection_point_index);

							// Advance.
							came_in_through_vertex = false;
							last_edges_2_valid = true;
							last_edges_2 = std::unordered_set<size_t>{ h_index_2[around_face_2->opposite()] };
							f_2 = around_face_2->opposite()->facet(); // we flip across the edge 'around_face_2' to the next face
							break;
						}
						else { // h_1 did not intersect around_face_2.
							if (around_face_2 == around_face_2_end) {
								// We've tried all the edges and failed, meaning that h_1 terminates inside f_2.
								done_with_h_1 = true; // This will cause a break out of the outer while loop.
								break;
							}
							around_face_2 = around_face_2->next(); // There's still another edge of f_2 left to try, do it.
						}
					}
				}
			}
			/*
			else if (we_have_terminated_at_a_vertex_of_mesh_2) {
				done_with_h_1 = true;
			}
			else { // !we_have_terminated_at_a_vertex_of_mesh_2 && we_are_intersecting_a_vertex_on_way_out

			}
			*/
		} // We are almost done with halfedge h_1. 
		  // Now just enqueue all the not-yet-seen edges leaving from its endpoint to be examined.

		VCI endpoint = h_1->vertex();
		HVCC hv = endpoint->vertex_begin();
		HVCC hv_end{ hv };

		size_t debug_counter = 0;
		// For all edges starting at the endpoint of h_1, push them to stack
		do {
			debug_counter++;
			Halfedge starting_from_endpoint = hv->opposite();
			if (!set_contains<size_t>(used, h_index_1[starting_from_endpoint])) {
				halfedge_stack.push_front(starting_from_endpoint);
				face_stack.push_front(f_2);

				// bookkeeping
				used.insert(h_index_1[starting_from_endpoint]);
				used.insert(h_index_1[starting_from_endpoint->opposite()]);
			}
			++hv;
		} while (hv != hv_end);


		// Debug progress.
		for (size_t percent = 1; percent <= 100; percent++) {
			size_t threshold = percent * mesh_1.size_of_halfedges() / 100;
			if (2*debug_counter1 == threshold || 2*debug_counter1 + 1 == threshold) {
				// We increment debug_counter1 once per popped halfedge, but we never do both a hedge and its twin, so we take the counter times two when comparing against mesh.size_of_halfedges.
				vout << "collect_intersection_events has processed ~" << (percent < 10 ? "0" : "") << percent << "% edges of mesh_1." << std::endl;
				break;
			}
		}
		// End of debug.
	}


	// We have gathered all events of type edge-edge, edge-vertex, or ends-at-edge.
	// Now we add all events of type ends-at-vertex for mesh_1.
	HCI h_it_1 = mesh_1.halfedges_begin();
	for (size_t i = 0; i < mesh_1.size_of_halfedges(); i++, h_it_1++) {
		size_t end_of_h_it_1 = v_index_1[h_it_1->vertex()];
		if (vertex_overlap_map_1.find(end_of_h_it_1) != vertex_overlap_map_1.end()) {
			for (size_t v_2 : vertex_overlap_map_1.at(end_of_h_it_1)) {
				IE_ENDS_ON_VERT intersection_event_mesh_1{ h_it_1->vertex()->point(), v_2 };
				insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_1, h_index_1[h_it_1], intersection_event_mesh_1);
			}
		}
	}
	// Now we add all events of type ends-at-vertex for mesh_2.
	HCI h_it_2 = mesh_2.halfedges_begin();
	for (size_t i = 0; i < mesh_2.size_of_halfedges(); i++, h_it_2++) {
		size_t end_of_h_it_2 = v_index_2[h_it_2->vertex()];
		if (vertex_overlap_map_2.find(end_of_h_it_2) != vertex_overlap_map_2.end()) {
			for (size_t v_1 : vertex_overlap_map_2.at(end_of_h_it_2)) {
				IE_ENDS_ON_VERT intersection_event_mesh_2{ h_it_2->vertex()->point(), v_1 };
				insert_into_map<size_t, IntersectionEvent>(out_intersection_events.intersection_events_2, h_index_2[h_it_2], intersection_event_mesh_2);
			}
		}
	}

	// Sanitize events. Return values are only used for debugging.
	size_t sanitized_pairs = sanitize_events_around_vv_overlaps(mesh_1, mesh_2, out_intersection_events, out_intersections_data, out_intersection_points,
		vertex_overlap_map_1, vertex_overlap_map_2, out_hedges_through_vertex_1, out_hedges_through_vertex_2,
		v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout);
	size_t sanitized_pairs_2 = sanitize_events_around_vv_overlaps_2(mesh_1, mesh_2, out_intersection_events, out_intersections_data, out_intersection_points,
		vertex_overlap_map_1, vertex_overlap_map_2, out_hedges_through_vertex_1, out_hedges_through_vertex_2,
		v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout);
	size_t sanitized_evs = sanitize_ev_events_around_vv_overlaps(mesh_1, mesh_2, out_intersection_events, out_intersections_data, out_intersection_points,
		vertex_overlap_map_1, vertex_overlap_map_2, out_hedges_through_vertex_1, out_hedges_through_vertex_2,
		v_index_1, h_index_1, f_index_1, v_index_2, h_index_2, f_index_2, vout);

	// The first pass should actually catch all the EE events that need to be gotten rid of.
	CGAL_assertion(sanitized_pairs_2 == 0);

	// Sort events.
	sort_intersection_events(mesh_1, mesh_2, out_intersection_events, vout);

	// Debug.
	debug_print_intersection_events(debug_counter1, out_intersection_events, vout);
}