#include "intersections_deprecated.h"

#include "geometry_utils.h"
#include "mesh_merger.h"
#include "cpp_utils.h"

using namespace std;

void intersect_asymmetric_deprecated(Polyhedron& mesh_1, Polyhedron& mesh_2,
	std::unordered_map<size_t /*HCI*/, std::vector<size_t /*HCI*/>>&       out_intersections_of_edges_1,
	std::unordered_map<size_t /*VCI*/, std::vector<size_t /*VCI*/>>&       out_verts_of_containing_face_1,
	std::unordered_map<size_t /*HCI*/, std::unordered_map<size_t /*HCI*/, size_t /*into intersection_points*/>>& out_intersections_data,
	std::vector<Point>& out_intersection_points,
	bool collect_intersections,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2,
	CGAL::Verbose_ostream& vout) {

#if 0
	Halfedge h_1 = mesh_1.halfedges_begin(); // an arbitrary halfedge
	Vertex v_1 = h_1->opposite()->vertex(); // the origin of h_1
	Facet f_2 = get_facet(mesh_2, containing_face(v_1->point(), mesh_2, f_index_2));

	out_intersections_of_edges_1.clear();
	out_verts_of_containing_face_1.clear();
	if (collect_intersections) {
		out_intersections_data.clear();
		out_intersection_points.clear();
	}

	// TODO: Make sure all halfedges of mesh_1 are _used_ by the end!

	std::deque<Halfedge> halfedge_stack{};
	halfedge_stack.push_front(h_1);
	std::deque<Facet> face_stack{};
	face_stack.push_front(f_2);
	std::unordered_set<size_t> used{};

	// bookkeeping
	used.insert(h_index_1[h_1]);
	used.insert(h_index_1[h_1->opposite()]);

	// Debug.
	size_t debug_counter1 = 0;

	size_t last_edge_2 = -1; // an arbitrary halfedge for initialization, it is invalid now anyway so the starting value won't be used
	bool last_edge_2_valid = false;
	while (!halfedge_stack.empty()) {
		// pop halfedge of mesh_1 and face of mesh_2 containing the origin of that halfedges from stack
		h_1 = stack_pop<Halfedge>(halfedge_stack);
		f_2 = stack_pop<Facet>(face_stack);

		// debug
		debug_counter1++;
		auto debug_popped = f_index_2[f_2];
		auto debug_correct = containing_face(h_1->opposite()->vertex()->point(), mesh_2, f_index_2);

		// Temporary!
		f_2 = get_facet(mesh_2, containing_face(h_1->opposite()->vertex()->point(), mesh_2, f_index_2));

		//CGAL_assertion(debug_popped == debug_correct); // TODO: Get this assertion to pass without the hack!
		//vout << "Popped halfedge " << h_index_1[h_1] << endl;

		// Write start vertex of h_1 to out_verts_of_containing_face_1 if necessary
		size_t start_vertex_of_h_1 = get_start_vert(mesh_1, h_index_1[h_1], v_index_1);
		// check to avoid duplication
		if (out_verts_of_containing_face_1.find(start_vertex_of_h_1) == out_verts_of_containing_face_1.end()) {
			HCI h_it = f_2->halfedge();
			HCI h_it_end = f_2->halfedge();
			do {
				insert_into_map<size_t, size_t>(out_verts_of_containing_face_1, start_vertex_of_h_1, v_index_2[h_it->vertex()]);
				h_it = h_it->next();
			} while (h_it != h_it_end);
		}

		// find intersections of h_1 with edges of mesh_2
		HCI around_face_2 = f_2->halfedge()->next();
		HCI around_face_2_end = f_2->halfedge();
		while (true) {
			Point intersection_point;
			if ((!last_edge_2_valid || h_index_2[around_face_2] != last_edge_2) && intersect(around_face_2, h_1, intersection_point)) {

				// Debug.
				//vout << "intersected " << h_index_1[h_1] << " with " << h_index_2[around_face_2] << std::endl;
				// CGAL_assertion(intersect(h_1, around_face_2, intersection_point));

				last_edge_2_valid = true;
				last_edge_2 = h_index_2[around_face_2->opposite()]; // right?
				insert_into_map<size_t, size_t>(out_intersections_of_edges_1, h_index_1[h_1], h_index_2[around_face_2]);

				if (collect_intersections) {
					out_intersection_points.push_back(intersection_point);
					size_t intersection_point_index = out_intersection_points.size() - 1;
					if (out_intersections_data.find(h_index_1[h_1]) == out_intersections_data.end()) {
						out_intersections_data[h_index_1[h_1]] = std::unordered_map<size_t, size_t>{};
					}
					out_intersections_data[h_index_1[h_1]][h_index_2[around_face_2]] = intersection_point_index;
				}

				// TODO Handle case where edge of mesh_1 goes through vertex of mesh_2
				f_2 = around_face_2->opposite()->facet();
				around_face_2 = f_2->halfedge()->next();
				around_face_2_end = f_2->halfedge();
			}
			else {
				if (around_face_2 == around_face_2_end)
					break;
				around_face_2 = around_face_2->next();
			}
		}

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

		// Write end vertex of h_1 to out_verts_of_containing_face_1 if necessary.
		// This is required because if we only wrote to that data structure for start vertices of halfedges, there could be vertices that are never written.
		size_t end_vertex_of_h_1 = v_index_1[endpoint];
		// check to avoid duplication
		if (out_verts_of_containing_face_1.find(end_vertex_of_h_1) == out_verts_of_containing_face_1.end()) {
			HCI h_it = f_2->halfedge();
			HCI h_it_end = f_2->halfedge();
			do {
				insert_into_map<size_t, size_t>(out_verts_of_containing_face_1, end_vertex_of_h_1, v_index_2[h_it->vertex()]);
				h_it = h_it->next();
			} while (h_it != h_it_end);
		}

		//vout << "hedge " << h_index_1[h_1] << "(" << h_index_1[h_1->opposite()] << ") has " << debug_counter << " neighbours.\n";
	}
	vout << "FINISHED intersect_asymmetric(), investigated " << debug_counter1 << " edges." << endl;
}

// To be used after a single call to intersect_assymetric that generated intersections_of_edges_1.
// This will generate intersections_of_edges_2 that are consistent with intersections_of_edges_1.
// This will also fill out the 'out_verts_of_containing_face_2'.
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
	CGAL::Verbose_ostream& vout) {

	// Do the out_verts_of_containing_face_2 for each vertex v of mesh_2.
	for (size_t v = 0; v < mesh_2.size_of_vertices(); v++) {
		size_t f_1 = containing_face(get_vertex(mesh_2, v)->point(), mesh_1, f_index_1);
		FCI facet_1 = get_facet(mesh_1, f_1);
		// check to avoid duplication
		if (out_verts_of_containing_face_2.find(v) == out_verts_of_containing_face_2.end()) {
			HCI h_it = facet_1->halfedge();
			HCI h_it_end = facet_1->halfedge();
			do {
				insert_into_map<size_t, size_t>(out_verts_of_containing_face_2, v, v_index_1[h_it->vertex()]);
				h_it = h_it->next();
			} while (h_it != h_it_end);
		}
	}

	// Remember which h_1s got their neighbours stored in out_intersections_of_edges_2 below.
	std::unordered_set<size_t> h_1s_flipped{};

	// Generate the intersections.
	for (auto kvp : intersections_of_edges_1) {
		size_t h_1 = kvp.first;
		std::vector<size_t> h_2s = kvp.second;
		for (size_t h_2 : h_2s) {
			// Make sure for a single edge all the data goes into one of the halfedges, instead of being split.
			size_t h_2_opp = h_index_2[get_halfedge(mesh_2, h_2)->opposite()];
			size_t h_2_to_use = min(h_2, h_2_opp);

			size_t h_1_to_store = h_1;
			if (h_2 == h_2_to_use) {
				// Preserve orientation properties, i.e.
				// If intersections_of_edges_1 stores that          1->0 is intersected by 6->4
				// Then intersections_of_edges_2 either stores that 6->4 is intersected by 0->1
				// Or that                                          4->6 is intersected by 1->0
				h_1_to_store = h_index_1[get_halfedge(mesh_1, h_1)->opposite()];
				h_1s_flipped.insert(h_1_to_store);
			}
			insert_into_map<size_t, size_t>(out_intersections_of_edges_2, min(h_2, h_2_opp), h_1_to_store);
		}
	}

	// Now we need to sort the intersections for each edge of mesh_2.
	for (auto kvp : out_intersections_of_edges_2) {
		size_t h_2 = kvp.first;
		std::vector<size_t> h_1s = kvp.second;
		Point start_point_of_h_2 = get_halfedge(mesh_2, h_2)->opposite()->vertex()->point();
		Point end_point_of_h_2 = get_halfedge(mesh_2, h_2)->vertex()->point();

		std::vector<Point> p_1s{}; // mimicks h_1s but stores corresponding intersection points
		for (size_t h_1 : h_1s) {

			// Flip h_1 back if necessary for lookup in intersections_data.
			if (h_1s_flipped.find(h_1) != h_1s_flipped.end()) {
				h_1 = h_index_1[get_halfedge(mesh_1, h_1)->opposite()];
			}

			// Flip h_2 if necessary for lookup in intersections_data.
			if (intersections_data[h_1].find(h_2) == intersections_data[h_1].end()) {
				h_2 = h_index_2[get_halfedge(mesh_2, h_2)->opposite()];
			}
			CGAL_assertion(intersections_data[h_1].find(h_2) != intersections_data[h_1].end());

			Point p_intersection = intersection_points[intersections_data[h_1][h_2]];
			p_1s.push_back(p_intersection);
		}

		// Choose which dimension we will look at to sort the intersection points.
		// Pick the one along which the edge h_2 is the longest.
		int xyz = 0;
		double max_diff = max(max(abs(start_point_of_h_2.x() - end_point_of_h_2.x()), abs(start_point_of_h_2.y() - end_point_of_h_2.y())), abs(start_point_of_h_2.z() - end_point_of_h_2.z()));
		if (max_diff == abs(start_point_of_h_2.x() - end_point_of_h_2.x())) xyz = 0;
		else if (max_diff == abs(start_point_of_h_2.y() - end_point_of_h_2.y())) xyz = 1;
		else xyz = 2;

		std::vector<std::pair<double, size_t>> sorter{};
		for (size_t i = 0; i < p_1s.size(); i++) {
			double p_coord = p_1s[i].x();
			if (xyz == 1) p_coord = p_1s[i].y();
			else if (xyz == 2) p_coord = p_1s[i].z();
			sorter.push_back(make_pair(p_coord, i));
		}
		sort(sorter.begin(), sorter.end());

		double start_coord = start_point_of_h_2.x();
		if (xyz == 1) start_coord = start_point_of_h_2.y();
		else if (xyz == 2) start_coord = start_point_of_h_2.z();
		double end_coord = end_point_of_h_2.x();
		if (xyz == 1) end_coord = end_point_of_h_2.y();
		else if (xyz == 2) end_coord = end_point_of_h_2.z();

		if (start_coord > end_coord) {
			std::reverse(sorter.begin(), sorter.end());
		}

		// Insert intersections in the right order this time.
		std::vector<size_t> correct_order_intersections_for_h_2{};
		size_t h_2_for_lookup = h_2;
		for (auto pair : sorter) {
			size_t ind = pair.second;

			// Flip h_2 if necessary for lookup in intersections_data.
			
			if (out_intersections_of_edges_2.find(h_2_for_lookup) == out_intersections_of_edges_2.end()) {
				h_2_for_lookup = h_index_2[get_halfedge(mesh_2, h_2_for_lookup)->opposite()];
			}
			CGAL_assertion(out_intersections_of_edges_2.find(h_2_for_lookup) != out_intersections_of_edges_2.end());

			correct_order_intersections_for_h_2.push_back(out_intersections_of_edges_2[h_2_for_lookup][ind]);
		}

		// Overwrite intersections with the correctly ordered list.
		CGAL_assertion(correct_order_intersections_for_h_2.size() == out_intersections_of_edges_2[h_2_for_lookup].size());
		out_intersections_of_edges_2[h_2_for_lookup] = correct_order_intersections_for_h_2;
	}
#endif
}