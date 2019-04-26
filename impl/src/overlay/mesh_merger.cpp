#include "mesh_merger.h"

#include "../common/constants.h"
#include "../embedding/io.h"
#include "../embedding/geometry_utils.h"
#include "cpp_utils.h"
#include "debug_utils.h"
#include "geometry_utils.h"

#include <CGAL/Timer.h>

// Assumes that ha ends at a vertex v, such that
// hb starts at a vertex v' that overlaps with v.
bool do_hedges_overlap(Polyhedron& mesh_1, Polyhedron& mesh_2, size_t ha_index, bool ha_is_from_mesh_1, // then hb is from mesh_2
	size_t hb_index, const IntersectionEventStore& intersection_events,
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	const std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>>& vertex_overlap_map_2,
	const VInvIndex& v_index_1, const VInvIndex& v_index_2) {

	Polyhedron& mesh_a = ha_is_from_mesh_1 ? mesh_1 : mesh_2;
	Polyhedron& mesh_b = ha_is_from_mesh_1 ? mesh_2 : mesh_1;
	const VInvIndex& v_index_a = ha_is_from_mesh_1 ? v_index_1 : v_index_2;
	const VInvIndex& v_index_b = ha_is_from_mesh_1 ? v_index_2 : v_index_1;

	const std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>>& vertex_overlap_map_a = ha_is_from_mesh_1 ? vertex_overlap_map_1 : vertex_overlap_map_2;
	const std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>>& vertex_overlap_map_b = ha_is_from_mesh_1 ? vertex_overlap_map_2 : vertex_overlap_map_1;

	const auto& intersection_events_a = ha_is_from_mesh_1 ? intersection_events.intersection_events_1 : intersection_events.intersection_events_2;
	const auto& intersection_events_b = ha_is_from_mesh_1 ? intersection_events.intersection_events_2 : intersection_events.intersection_events_1;

	size_t ha_start_vert_ind = get_start_vert(mesh_a, ha_index, v_index_a);
	size_t hb_end_vert_ind = get_end_vert(mesh_b, hb_index, v_index_b);

	// Step 1/3. Check if hb ends at start of ha:
	if (vertex_overlap_map_a.find(ha_start_vert_ind) != vertex_overlap_map_a.end()) {
		for (size_t vb : vertex_overlap_map_a.at(ha_start_vert_ind)) {
			if (vb == hb_end_vert_ind) {
				return true;
			}
		}
	}

	// Step 2/3. Check if hb has an EV intersection with the start of ha:
	for (IntersectionEvent ie : intersection_events_b.at(hb_index)) {
		if (ie.tag == IntersectionEvent::IET_EV) {
			if (ha_start_vert_ind == ie.other_vert_index()) {
				return true;
			}
		}
	}

	// Step 3/3. Check if ha has an EV intersection with the end of hb:
	for (IntersectionEvent ie : intersection_events_a.at(ha_index)) {
		if (ie.tag == IntersectionEvent::IET_EV) {
			if (hb_end_vert_ind == ie.other_vert_index()) {
				return true;
			}
		}
	}

	// We haven't found any events that would signal there's an overlap.
	return false;
}


bool halfedge_or_opposite_is_intersected(Polyhedron& mesh, size_t halfedge, const HInvIndex& h_index, std::unordered_map<size_t, std::vector<size_t>>& intersections_of_edges) {
	HCI h = get_halfedge(mesh, halfedge);
	size_t ind1 = h_index[h]; // == halfedge
	size_t ind2 = h_index[h->opposite()];
	if (intersections_of_edges.find(ind1) != intersections_of_edges.end() && !intersections_of_edges[ind1].empty())
		return true;
	if (intersections_of_edges.find(ind2) != intersections_of_edges.end() && !intersections_of_edges[ind2].empty())
		return true;
	return false;
}

size_t/*vert*/ read_intersections_data(Polyhedron& mesh_1, Polyhedron& mesh_2, size_t ind1, size_t ind2, const HInvIndex& h_index_1, const HInvIndex& h_index_2,
	std::unordered_map<size_t /*HCI*/, std::unordered_map<size_t /*HCI*/, size_t /*into intersection_points*/>>& intersections_data) {
	if (intersections_data.find(ind1) == intersections_data.end())
		ind1 = get_opposite(mesh_1, ind1, h_index_1);
	CGAL_assertion(intersections_data.find(ind1) != intersections_data.end());

	if (intersections_data[ind1].find(ind2) == intersections_data[ind1].end())
		ind2 = get_opposite(mesh_2, ind2, h_index_2);
	CGAL_assertion(intersections_data[ind1].find(ind2) != intersections_data[ind1].end());

	return intersections_data.at(ind1).at(ind2);
}

// Returns (true, e') if e' intersects e and that is the next intersection after the intersection of e and "intersecting_edge".
// Returns (false, v) if that was the last intersection, and v is the following endpoint of e.
std::pair<bool /*res_is_edge*/, size_t> next_intersection(Polyhedron& mesh, Polyhedron& mesh_other, size_t e, const HInvIndex& h_index, const VInvIndex& v_index, const HInvIndex& h_index_other,
	size_t intersecting_edge,
	size_t intersecting_edge_twin,
	std::unordered_map<size_t /*HCI*/, std::vector<size_t /*HCI*/>>& intersections_of_edges) {

	bool reversed = false;
	if (intersections_of_edges.find(e) == intersections_of_edges.end()) {
		reversed = true;
		e = get_opposite(mesh, e, h_index);
	}

	std::vector<size_t /*HCI*/> intersections_of_e = intersections_of_edges[e];
	size_t i = 0;
	if (intersections_of_e.empty()) {
		CGAL_assertion(false);
		// The intersection we were looking for was not present.
	}
	for (; i < intersections_of_e.size(); i++) {
		if (intersections_of_e[i] == intersecting_edge || intersections_of_e[i] == intersecting_edge_twin)
			break;
		if (i == intersections_of_e.size() - 1)
			CGAL_assertion(false);
			// The intersection we were looking for was not present.
	}
	if (!reversed) {
		if (i != intersections_of_e.size() - 1) {
			return std::make_pair(true, intersections_of_e[i + 1]);
		}
		else {
			return std::make_pair(false, get_end_vert(mesh, e, v_index));
		}
	}
	else {
		if (i != 0) {
			return std::make_pair(true, get_opposite(mesh_other, intersections_of_e[i - 1], h_index_other));
		}
		else {
			return std::make_pair(false, get_start_vert(mesh, e, v_index));
		}
	}
}

// Returns (true, e') if e' intersects e and that is the first intersection for e.
// Returns (false, v) if there are no intersections, and v is the endpoint of e.
std::pair<bool /*res_is_edge*/, size_t> first_intersection(Polyhedron& mesh, Polyhedron& mesh_other, size_t e, const HInvIndex& h_index, const VInvIndex& v_index, const HInvIndex& h_index_other,
	std::unordered_map<size_t /*HCI*/, std::vector<size_t /*HCI*/>>& intersections_of_edges) {

	bool reversed = false;
	if (intersections_of_edges.find(e) == intersections_of_edges.end()) {
		reversed = true;
		e = get_opposite(mesh, e, h_index);
		if (intersections_of_edges.find(e) == intersections_of_edges.end()) {
			// no intersections whatsoever
			return std::make_pair(false, get_start_vert(mesh, e, v_index));
		}
	}

	std::vector<size_t /*HCI*/> intersections_of_e = intersections_of_edges[e];
	if (intersections_of_e.empty())
		// no intersections whatsoever
		return std::make_pair(false, get_start_vert(mesh, e, v_index));

	if (!reversed) {
		return std::make_pair(true, intersections_of_e[0]);
	}
	else {
		return std::make_pair(true, get_opposite(mesh_other, intersections_of_e[intersections_of_e.size() - 1], h_index_other));
	}
}

// Deduplicate faces stored in a hashtable and return them in a single list of faces.
std::vector<std::vector<size_t>> deduplicate(std::unordered_map<size_t, std::vector<NormalizedVertexList>> faces, CGAL::Verbose_ostream& vout) {
	std::vector<std::vector<size_t>> result{};
	for (auto kvp : faces) {
		std::vector<NormalizedVertexList> some_faces = kvp.second;
		std::unordered_set<size_t> indices_to_ignore{};
		for (size_t i = 0; i < some_faces.size(); i++) {
			if (set_contains<size_t>(indices_to_ignore, i)) 
				continue;
			for (size_t j = i + 1; j < some_faces.size(); j++) {
				if (set_contains<size_t>(indices_to_ignore, j))
					continue;
				if (some_faces[i].equals(some_faces[j]))
					indices_to_ignore.insert(j);
			}
			result.push_back(some_faces[i].get_face());
		}
	}

	// Debug.
	vout << "Deduplicated faces:" << std::endl;
	for (std::vector<size_t> face : result) {
		for (size_t vind : face) {
			vout << " " << vind;
		}
		vout << std::endl;
	} vout << std::endl;

	return result;
}

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
	CGAL::Verbose_ostream& vout)
{
	size_t debug_overlap_counter = 0;
	size_t debug_max_overlap_size = 0;

	VCI v1 = mesh_1.vertices_begin();
	for (size_t i = 0; i < V1; i++, v1++) {
		VCI v2 = mesh_2.vertices_begin();
		for (size_t j = 0; j < V2; j++, v2++) {
			Vector v1pos = v1->point() - CGAL::ORIGIN;
			Vector v2pos = v2->point() - CGAL::ORIGIN;
			if (identify_vertices(v1pos, v2pos)) {
				insert_into_map<size_t, size_t>(out_vertex_overlap_map_1, i, j);
				insert_into_map<size_t, size_t>(out_vertex_overlap_map_2, j, i);
				debug_overlap_counter++;
				debug_max_overlap_size = std::max(debug_max_overlap_size, std::max(out_vertex_overlap_map_1[i].size(), out_vertex_overlap_map_2[j].size()));
			}
		}
	}
	vout << "The call to generate_vv_overlaps() is done. " << debug_overlap_counter << " overlaps were generated between the two meshes." << std::endl;
	vout << "A record vertex has overlapped with " << debug_max_overlap_size << " vertices of the other mesh." << std::endl;
}

// Vertices of mesh_1 will all persist.
// Vertices of mesh_2 that overlap them shall all perish. Their indices need to be mapped onto the mesh_1 indices.
// The remaining vertices of mesh_2 may shift down. That needs to be reflected in this mapping.
std::unordered_map<size_t, size_t> generate_vert_index_mapping_for_merged_mesh(
	const std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& vertex_overlap_map_1,
	const std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>>& vertex_overlap_map_2,
	const size_t V1, // original vertex count of mesh_1
	const size_t V2, // original vertex count of mesh_2
	const size_t V3  // count of intersection vertices
) {
	std::unordered_map<size_t, size_t> res{};
	for (size_t i = 0; i < V1 + V2 + V3; i++) {
		res[i] = i;
	}

	// We need to go through the pairs in the order of increasing v2!
	// So let's sort.
	std::vector<std::pair<size_t, size_t>> reversed_overlaps{};
	for (auto kvp : vertex_overlap_map_1) {
		size_t v1 = kvp.first;
		for (size_t v2 : kvp.second) {
			reversed_overlaps.push_back(std::make_pair(v2, v1));
		}
	}
	std::sort(reversed_overlaps.begin(), reversed_overlaps.end());

	for (std::pair<size_t, size_t> reversed_pair : reversed_overlaps) {
		size_t v1 = reversed_pair.second;
		size_t v2 = reversed_pair.first;

		size_t v2_offset = V1 + v2;
		res[v2_offset] = v1;
		for (size_t j = v2_offset + 1; j < V1 + V2 + V3; j++) {
			res[j] = res[j] - 1;
		}
	}

	// Debug.
	for (size_t v1 = 0; v1 < V1; v1++) {
		CGAL_assertion(res[v1] == v1);
	}
	for (auto kvp : vertex_overlap_map_1) {
		size_t v1 = kvp.first;
		CGAL_assertion(v1 < V1);
		for (size_t v2 : kvp.second) {
			CGAL_assertion(res[V1+v2] == v1);
		}
	}

	return res;
}
size_t normalize_vertex_index_based_on_overlaps(size_t v_ind,
	const std::unordered_map<size_t, size_t> vert_index_mapping_for_merged_mesh
) {
	return vert_index_mapping_for_merged_mesh.at(v_ind);
}

// Returns halfedges of a mesh leaving a vertex of that mesh.
std::vector<size_t> get_halfedges_leaving_vert(Polyhedron& mesh, size_t vert, const HInvIndex& h_index) {
	std::vector<size_t> res{};
	HVCC hv = get_vertex(mesh, vert)->vertex_begin();
	HVCC hv_end{ hv };
	do {
		HCI halfedge_leaving_vert = hv->opposite();
		res.push_back(h_index[halfedge_leaving_vert]);
	} while (++hv != hv_end);
	return res;
}

Vector get_halfedge_dir(HCI& h) {
	return normalized(h->vertex()->point() - h->opposite()->vertex()->point());
}

// Project arg onto the plane perpendicular to normal.
Vector get_perpendicular_component(const Vector& normal, const Vector& arg) {
	return arg - normal * ((arg*normal) / (normal*normal));
}

// Return the rotation that maps v to (1,0,0). TODO: Verify this works.
Transformation get_rotation_to_positive_x_axis(Vector v) {
	double around_x = asin(v.y());
	double around_y = acos(v.x()/cos(around_x));
	return (rotY(around_y)*rotX(around_x)).inverse();
}

// Map the halfedge onto the plane tangent to the unit sphere at vert_pos, and then calculate its angle (between -pi and pi) around the normal of that plane.
double get_angle_of_halfedge_relative_to_vert(Vector vert_pos, Vector halfedge_dir) {
	Transformation t = get_rotation_to_positive_x_axis(vert_pos);
	Vector proj = get_perpendicular_component(vert_pos, halfedge_dir);
	Vector yz = t.transform(proj);
	double incoming_angle = atan2(yz.z(), yz.y());
	return incoming_angle;
}



// Core of the method. Projects the candidate_dirs and the reference_dir onto a plane tangent to the unit sphere at the 'normal' vector, and sorts them by angle.
#define GO_THE_OTHER_WAY false
size_t get_neighbouring_halfedge_by_angle(const Vector& normal, const std::vector<std::pair<Vector, size_t>>& candidate_dirs, const Vector& reference_dir, bool go_the_other_way) {

	double incoming_angle = get_angle_of_halfedge_relative_to_vert(normal, reference_dir);

	// Get the angles of the candidate halfedges.
	std::vector<std::pair<double, size_t>> candidate_angles{};
	for (auto kvp : candidate_dirs) {
		Vector candidate_dir = kvp.first;
		size_t candidate_index = kvp.second;
		double candidate_angle = get_angle_of_halfedge_relative_to_vert(normal, candidate_dir);
		candidate_angles.push_back(std::make_pair(candidate_angle, candidate_index));
	}

	std::sort(candidate_angles.begin(), candidate_angles.end());

	// Find the candidate halfedge which should be turned to after following the incoming edge.
	if (incoming_angle <= candidate_angles[0].first || incoming_angle > candidate_angles[candidate_angles.size() - 1].first) {
		if(go_the_other_way)
			return candidate_angles[0].second;
		else
			return candidate_angles[candidate_angles.size() - 1].second;
	}
	else {
		for (size_t i = 1; i < candidate_angles.size(); i++) {
			if (incoming_angle <= candidate_angles[i].first) {
				if(go_the_other_way)
					return candidate_angles[i].second;
				else
					return candidate_angles[(i - 1) % candidate_angles.size()].second;
			}
		}
	}

	// Unreachable.
	CGAL_assertion(false);
	return -1;
}


// Core of the method. Projects the candidate_dirs and the reference_dir onto a plane tangent to the unit sphere at the 'normal' vector, and sorts them by angle.
// This version stores an extra bool with each candidate.
std::pair<size_t, bool> get_neighbouring_halfedge_by_angle(const Vector& normal, const std::vector<std::pair<Vector, std::pair<size_t, bool>>>& candidate_dirs, const Vector& reference_dir,
	bool go_the_other_way) {

	double incoming_angle = get_angle_of_halfedge_relative_to_vert(normal, reference_dir);

	// Get the angles of the candidate halfedges.
	std::vector<std::pair<double, std::pair<size_t, bool>>> candidate_angles{};
	for (auto kvp : candidate_dirs) {
		Vector candidate_dir = kvp.first;
		std::pair<size_t, bool> candidate_index = kvp.second;
		double candidate_angle = get_angle_of_halfedge_relative_to_vert(normal, candidate_dir);
		candidate_angles.push_back(std::make_pair(candidate_angle, candidate_index));
	}

	std::sort(candidate_angles.begin(), candidate_angles.end());

	// Find the candidate halfedge which should be turned to after following the incoming edge.
	if (incoming_angle <= candidate_angles[0].first || incoming_angle > candidate_angles[candidate_angles.size() - 1].first) {
		if (go_the_other_way)
			return candidate_angles[0].second;
		else
			return candidate_angles[candidate_angles.size() - 1].second;
	}
	else {
		for (size_t i = 1; i < candidate_angles.size(); i++) {
			if (incoming_angle <= candidate_angles[i].first) {
				if(go_the_other_way)
					return candidate_angles[i].second;
				else
					return candidate_angles[(i - 1) % candidate_angles.size()].second;
			}
		}
	}

	// Unreachable.
	CGAL_assertion(false);
	return std::make_pair(-1, false);
}

// Wrapper.
// This method is used to handle EV intersections.
size_t get_neighbouring_halfedge_by_angle(Polyhedron& mesh, Polyhedron& mesh_of_incoming_halfedge, size_t vert, std::vector<size_t> candidate_halfedges, size_t incoming_halfedge,
	bool go_the_other_way) {

	// Get the direction of the incoming edge's twin.
	Vector normal = get_vertex(mesh, vert)->point() - CGAL::ORIGIN;
	Vector incoming_opposite_direction = get_halfedge_dir(get_halfedge(mesh_of_incoming_halfedge, incoming_halfedge)->opposite());
	
	// Get the candidate dirs.
	std::vector<std::pair<Vector, size_t>> candidate_dirs{};
	for (size_t h : candidate_halfedges) {
		Vector candidate_dir = get_halfedge_dir(get_halfedge(mesh, h));
		candidate_dirs.push_back(std::make_pair(candidate_dir, h));
	}

	// Call core method.
	return get_neighbouring_halfedge_by_angle(normal, candidate_dirs, incoming_opposite_direction, go_the_other_way);
}

// Wrapper.
// This method is used to handle ENDS_AT_VERT intersection events.
std::pair<size_t, bool/*is_from_mesh_1*/> get_neighbouring_halfedge_by_angle(Polyhedron& mesh_1, Polyhedron& mesh_2, const Point& vert, 
	std::vector<std::pair<size_t, bool/*is_from_mesh_1*/>> candidate_halfedges, const HCI& incoming_halfedge, bool go_the_other_way) {

	// Get the direction of the incoming edge's twin.
	Vector normal = vert - CGAL::ORIGIN;
	Vector incoming_opposite_direction = get_halfedge_dir(incoming_halfedge->opposite());

	// Get the candidate dirs.
	std::vector<std::pair<Vector, std::pair<size_t, bool>>> candidate_dirs{};
	for (std::pair<size_t, bool/*is_from_mesh_1*/> kvp : candidate_halfedges) {
		bool h_is_from_mesh_1 = kvp.second;
		Vector candidate_dir = get_halfedge_dir(get_halfedge(h_is_from_mesh_1 ? mesh_1 : mesh_2, kvp.first));
		candidate_dirs.push_back(std::make_pair(candidate_dir, kvp));
	}

	// Call core method.
	return get_neighbouring_halfedge_by_angle(normal, candidate_dirs, incoming_opposite_direction, go_the_other_way);
}

// Returns the index of the next vertex to be added to the face.
// Updates the seven variables passed as arguments immediately after 'intersection_events'
size_t collect_edge_starting_at_event(Polyhedron& mesh_1, Polyhedron& mesh_2,
	const IntersectionEventStore& intersection_events,

	size_t& last_added_vert,
	bool& last_added_vert_is_from_mesh_1, // this variable is to be meaningless when 'last_added_vert_is_an_intersection' is true
	bool& last_added_vert_is_an_intersection,
	size_t& halfedge_leading_to_last_added_vertex_index,
	bool& halfedge_leading_to_last_added_vertex_is_from_mesh_1,
	size_t& last_added_vert_event_index, // this variable is to be meaningless when 'last_added_vert_was_an_endpoint' is true
	bool& last_added_vert_was_an_endpoint,

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
	CGAL::Verbose_ostream& vout,
	bool go_the_other_way
) {
	if (!last_added_vert_was_an_endpoint) {
		IntersectionEvent last_event = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? intersection_events.intersection_events_1 : intersection_events.intersection_events_2).
			at(halfedge_leading_to_last_added_vertex_index)[last_added_vert_event_index];
		switch (last_event.tag) {
		case IntersectionEvent::IET_EE:
		{
			// The last added vertex was an edge-edge intersection.
			// This means we need to locate that e-e event on the other edge, and then investigate the successor event.
			// If there's no successor, then we investigate the endpoint of the other edge.
			size_t halfedge_leading_to_last_added_vertex_opposite_index = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? h_index_1 : h_index_2)[
				get_halfedge(halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? mesh_1 : mesh_2, halfedge_leading_to_last_added_vertex_index)->opposite()
			];
			size_t other_edge_index = last_event.other_edge_index();
			HCI other_edge = get_halfedge(halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? mesh_2 : mesh_1, other_edge_index);
			auto other_edge_events = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? intersection_events.intersection_events_2 : intersection_events.intersection_events_1).
				at(other_edge_index);
			bool found_the_event_on_the_other_edge = false;
			size_t i = 0;
			for ( ; i < other_edge_events.size(); i++) {
				IntersectionEvent other_edge_event = other_edge_events[i];
				if (other_edge_event.tag == IntersectionEvent::IET_EE && other_edge_event.other_edge_index() == halfedge_leading_to_last_added_vertex_opposite_index) {
					found_the_event_on_the_other_edge = true;
					break;
				}
			}
			CGAL_assertion(found_the_event_on_the_other_edge);
			if (i == other_edge_events.size() - 1 || other_edge_events[i + 1].tag == IntersectionEvent::IET_ENDS_ON_EDGE || other_edge_events[i + 1].tag == IntersectionEvent::IET_ENDS_ON_VERT) {
				// The next vertex will be the endpoint of the other edge.
				size_t last_added_vert_no_offset = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? v_index_2 : v_index_1)[other_edge->vertex()];
				last_added_vert = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? V1 : 0) + last_added_vert_no_offset;
				last_added_vert_is_from_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
				last_added_vert_is_an_intersection = false;
				halfedge_leading_to_last_added_vertex_index = other_edge_index;
				halfedge_leading_to_last_added_vertex_is_from_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
				last_added_vert_event_index = -1;
				last_added_vert_was_an_endpoint = true;
				return last_added_vert;
			}
			else {
				// The next vertex will be the (i+1)th event on the other_edge.
				// In order to update the state correctly we need to investigate that event closer.
				switch (other_edge_events[i + 1].tag) {
				case IntersectionEvent::IET_EE:
				{
					// The next vertex will be an intersection vertex again.
					size_t index1 = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? other_edge_events[i + 1].other_edge_index() : other_edge_index;
					size_t index2 = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? other_edge_index : other_edge_events[i + 1].other_edge_index();
					last_added_vert = V1 + V2 + read_intersections_data(mesh_1, mesh_2, index1, index2, h_index_1, h_index_2, intersections_data);
					last_added_vert_is_from_mesh_1 = (bool)(-1); // just to emphasize this value is meaningless now because 'last_added_vert_is_an_intersection' is true
					last_added_vert_is_an_intersection = true;
					halfedge_leading_to_last_added_vertex_index = other_edge_index;
					halfedge_leading_to_last_added_vertex_is_from_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
					last_added_vert_event_index = i + 1;
					last_added_vert_was_an_endpoint = false;
					return last_added_vert;
				}
					break;
				case IntersectionEvent::IET_EV:
				{
					// The next vertex will be a vertex of the mesh whose edge we read the previous event from.
					size_t last_added_vert_no_offset = other_edge_events[i + 1].other_vert_index();
					last_added_vert = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? 0 : V1) + last_added_vert_no_offset;
					last_added_vert_is_from_mesh_1 = last_added_vert_is_from_mesh_1;
					last_added_vert_is_an_intersection = false;
					halfedge_leading_to_last_added_vertex_index = other_edge_index;
					halfedge_leading_to_last_added_vertex_is_from_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
					last_added_vert_event_index = i + 1;
					last_added_vert_was_an_endpoint = false;
					return last_added_vert;
				}
					break;
				default:
					CGAL_assertion(false); return -1;
				}
			}
		}
			break;
		case IntersectionEvent::IET_EV:
		{
			// The last added vertex was an edge-vertex intersection.
			// This means we need to consider all halfedges of the other mesh leaving from that vertex and pick the right one to turn into 
			//   (i.e. the one counterclockwise from our incoming edge by angle).
			
			size_t intersected_vert = last_event.other_vert_index();
			bool intersected_vert_in_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
			Polyhedron& mesh_of_intersected_vert = intersected_vert_in_mesh_1 ? mesh_1 : mesh_2;
			Polyhedron& mesh_of_intersected_halfedge = intersected_vert_in_mesh_1 ? mesh_2 : mesh_1;
			const HInvIndex& h_index_of_intersected_vert = intersected_vert_in_mesh_1 ? h_index_1 : h_index_2;

			std::vector<size_t> halfedges_leaving_intersected_vert = get_halfedges_leaving_vert(mesh_of_intersected_vert, intersected_vert, h_index_of_intersected_vert);
			size_t next_halfedge_index = get_neighbouring_halfedge_by_angle(
				mesh_of_intersected_vert, mesh_of_intersected_halfedge, intersected_vert, halfedges_leaving_intersected_vert, halfedge_leading_to_last_added_vertex_index, go_the_other_way);
			HCI next_halfedge = get_halfedge(mesh_of_intersected_vert, next_halfedge_index);
			auto next_halfedge_events = (intersected_vert_in_mesh_1 ? intersection_events.intersection_events_1 : intersection_events.intersection_events_2).at(next_halfedge_index);


			// Note! The code below is copy-pasted (and adjusted) from above!!! Could be encapsulated into a helper function.

			// The last added vert was the start vert of 'next_halfedge'.
			if (next_halfedge_events.size() <= 1 || next_halfedge_events[1].tag == IntersectionEvent::IET_ENDS_ON_EDGE || next_halfedge_events[1].tag == IntersectionEvent::IET_ENDS_ON_VERT) {
				// The next vertex will be the endpoint of the other edge.
				size_t last_added_vert_no_offset = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? v_index_2 : v_index_1)[next_halfedge->vertex()];
				last_added_vert = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? V1 : 0) + last_added_vert_no_offset;
				last_added_vert_is_from_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
				last_added_vert_is_an_intersection = false;
				halfedge_leading_to_last_added_vertex_index = next_halfedge_index;
				halfedge_leading_to_last_added_vertex_is_from_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
				last_added_vert_event_index = -1;
				last_added_vert_was_an_endpoint = true;
				return last_added_vert;
			}
			else {
				// The next vertex will be the #1 event on the other_edge (#0 is the IE_START event).
				// In order to update the state correctly we need to investigate that event closer.
				switch (next_halfedge_events[1].tag) {
				case IntersectionEvent::IET_EE:
				{
					// The next vertex will be an intersection vertex again.
					size_t index1 = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? next_halfedge_events[1].other_edge_index() : next_halfedge_index;
					size_t index2 = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? next_halfedge_index : next_halfedge_events[1].other_edge_index();
					last_added_vert = V1 + V2 + read_intersections_data(mesh_1, mesh_2, index1, index2, h_index_1, h_index_2, intersections_data);
					last_added_vert_is_from_mesh_1 = (bool)(-1); // just to emphasize this value is meaningless now because 'last_added_vert_is_an_intersection' is true
					last_added_vert_is_an_intersection = true;
					halfedge_leading_to_last_added_vertex_index = next_halfedge_index;
					halfedge_leading_to_last_added_vertex_is_from_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
					last_added_vert_event_index = 1;
					last_added_vert_was_an_endpoint = false;
					return last_added_vert;
				}
				break;
				case IntersectionEvent::IET_EV:
				{
					// The next vertex will be a vertex of the mesh whose edge we read the previous event from.
					size_t last_added_vert_no_offset = next_halfedge_events[1].other_vert_index();
					last_added_vert = (halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? 0 : V1) + last_added_vert_no_offset;
					last_added_vert_is_from_mesh_1 = last_added_vert_is_from_mesh_1;
					last_added_vert_is_an_intersection = false;
					halfedge_leading_to_last_added_vertex_index = next_halfedge_index;
					halfedge_leading_to_last_added_vertex_is_from_mesh_1 = !halfedge_leading_to_last_added_vertex_is_from_mesh_1;
					last_added_vert_event_index = 1;
					last_added_vert_was_an_endpoint = false;
					return last_added_vert;
				}
				break;
				default:
					CGAL_assertion(false); return -1;
				}
			}
		}
			break;
		default:
			// The other cases should be impossible.
			CGAL_assertion(false); return -1;
			break;
		}
	}
	else { //if (last_added_vert_was_an_endpoint)

		// Last added vertex was an endpoint. That means we need to investigate all hedges going through this endpoint and choose the right one to follow.
		// Then we need to investigate its first non-start event.

		// This will store pairs (hedge_index, is_the_hedge_from_mesh_1) of hedges of the other mesh that go through the vertex (either terminating there or just passing through)
		// as well as the next() halfedge of the halfedge that brought us there.
		// It would be useful to have a lookup datastructure that's a map vertex_index->mesh_number->hedges_from_that_mesh_passing_through_vertex.
		std::vector<std::pair<size_t, bool>> hedges_through_vertex{};
		const HInvIndex& last_halfedge_h_index = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? h_index_1 : h_index_2;
		const HInvIndex& other_h_index = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? h_index_2 : h_index_1;
		Polyhedron& mesh_of_last_halfedge = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? mesh_1 : mesh_2;
		Polyhedron& other_mesh = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? mesh_2 : mesh_1;
		size_t next_halfedge_ind = last_halfedge_h_index[get_halfedge(mesh_of_last_halfedge, halfedge_leading_to_last_added_vertex_index)->next()];
		hedges_through_vertex.push_back(std::make_pair(next_halfedge_ind, halfedge_leading_to_last_added_vertex_is_from_mesh_1));
		
		// Cancel the vertex index offset:
		size_t last_added_vert_no_offset = (last_added_vert < V1) ? last_added_vert : last_added_vert - V1;

		// For every 'vert_of_other_mesh' vertex index of the other mesh that corresponds to our endpoint.
		const std::unordered_map<size_t, std::vector<size_t>>& overlap_map = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? vertex_overlap_map_1 : vertex_overlap_map_2;
		if (overlap_map.find(last_added_vert_no_offset) != overlap_map.end()) {
			for (size_t vert_of_other_mesh : overlap_map.at(last_added_vert_no_offset)) { // overlap map stores no-offset vert indices, but the last_added_vert has offset if it is from mesh_2
				// Add every edge leaving the 'vert_of_other_mesh' into the 'hedges_through_vertex' candidate set.
				// Do not add edges overlapping with the current edge - we wouldn't want a degenerate face with an angle 0 in it.
				for (size_t candidate_ind : get_halfedges_leaving_vert(other_mesh, vert_of_other_mesh, other_h_index)) {
					if (!do_hedges_overlap(mesh_1, mesh_2, halfedge_leading_to_last_added_vertex_index, halfedge_leading_to_last_added_vertex_is_from_mesh_1,
						candidate_ind, intersection_events, vertex_overlap_map_1, vertex_overlap_map_2, v_index_1, v_index_2)) 
					{
						hedges_through_vertex.push_back(std::make_pair(candidate_ind, !halfedge_leading_to_last_added_vertex_is_from_mesh_1));
						std::cerr << "[1] hedges_through_vertex.push_back: " << candidate_ind << ", " << !halfedge_leading_to_last_added_vertex_is_from_mesh_1 << std::endl;
					}
				}
			}
		}

		// Also include edges of the other mesh that are just passing through that vertex in the 'hedges_through_vertex' candidate set. 
		auto corresponding_hedges_through_vertex = halfedge_leading_to_last_added_vertex_is_from_mesh_1 ? hedges_through_vertex_1 : hedges_through_vertex_2;
		if (corresponding_hedges_through_vertex.find(last_added_vert_no_offset) != corresponding_hedges_through_vertex.end())
		{
			std::vector<size_t> hedges_passing_through = corresponding_hedges_through_vertex.at(last_added_vert_no_offset);
			for (size_t hedge_passing_through_vertex : hedges_passing_through) {
				hedges_through_vertex.push_back(std::make_pair(hedge_passing_through_vertex, !halfedge_leading_to_last_added_vertex_is_from_mesh_1));
				std::cerr << "[2] hedges_through_vertex.push_back: " << hedge_passing_through_vertex << ", " << !halfedge_leading_to_last_added_vertex_is_from_mesh_1 << std::endl;
			}
		}
	
		// Get the next edge to be followed...
		Point last_added_vert_pos = get_vertex(last_added_vert_is_from_mesh_1 ? mesh_1 : mesh_2, last_added_vert_no_offset)->point();
		HCI halfedge_leading_to_last_added_vertex = get_halfedge(mesh_of_last_halfedge, halfedge_leading_to_last_added_vertex_index);
		auto kvp = get_neighbouring_halfedge_by_angle(mesh_1, mesh_2, last_added_vert_pos, hedges_through_vertex, halfedge_leading_to_last_added_vertex, go_the_other_way);
		size_t next_halfedge_index = kvp.first;
		bool next_halfedge_is_from_mesh_1 = kvp.second;
		HCI next_halfedge = get_halfedge(next_halfedge_is_from_mesh_1 ? mesh_1 : mesh_2, next_halfedge_index);
		auto next_halfedge_events = (next_halfedge_is_from_mesh_1 ? intersection_events.intersection_events_1 : intersection_events.intersection_events_2).at(next_halfedge_index);

		// ...and determine which event on it is the next one.


		// Note! The code below is copy-pasted (and adjusted) from above!!! Could be encapsulated into a helper function.
		// Event #0 is IE_START, we want to look at event #1.
		if (next_halfedge_events.size() <= 1 || next_halfedge_events[1].tag == IntersectionEvent::IET_ENDS_ON_EDGE || next_halfedge_events[1].tag == IntersectionEvent::IET_ENDS_ON_VERT) {
			size_t last_added_vert_no_offset = (next_halfedge_is_from_mesh_1 ? v_index_1 : v_index_2)[next_halfedge->vertex()];
			last_added_vert = (next_halfedge_is_from_mesh_1 ? 0 : V1) + last_added_vert_no_offset;
			last_added_vert_is_from_mesh_1 = next_halfedge_is_from_mesh_1;
			last_added_vert_is_an_intersection = false;
			halfedge_leading_to_last_added_vertex_index = next_halfedge_index;
			halfedge_leading_to_last_added_vertex_is_from_mesh_1 = next_halfedge_is_from_mesh_1;
			last_added_vert_event_index = -1;
			last_added_vert_was_an_endpoint = true;
			return last_added_vert;
		}
		else {
			switch (next_halfedge_events[1].tag) {
			case IntersectionEvent::IET_EE:
			{
				// The next vertex will be an intersection vertex.
				size_t index1 = next_halfedge_is_from_mesh_1 ? next_halfedge_index : next_halfedge_events[1].other_edge_index(); 
				size_t index2 = next_halfedge_is_from_mesh_1 ? next_halfedge_events[1].other_edge_index() : next_halfedge_index;
				last_added_vert = V1 + V2 + read_intersections_data(mesh_1, mesh_2, index1, index2, h_index_1, h_index_2, intersections_data);
				last_added_vert_is_from_mesh_1 = (bool)(-1); // just to emphasize this value is meaningless now because 'last_added_vert_is_an_intersection' is true
				last_added_vert_is_an_intersection = true;
				halfedge_leading_to_last_added_vertex_index = next_halfedge_index;
				halfedge_leading_to_last_added_vertex_is_from_mesh_1 = next_halfedge_is_from_mesh_1;
				last_added_vert_event_index = 1;
				last_added_vert_was_an_endpoint = false;
				return last_added_vert;
			}
			break;
			case IntersectionEvent::IET_EV:
			{
				// The next vertex will be a vertex of the other mesh than the one whose edge we are investigating.
				size_t last_added_vert_no_offset = next_halfedge_events[1].other_vert_index();
				last_added_vert = (next_halfedge_is_from_mesh_1 ? V1 : 0) + last_added_vert_no_offset;
				last_added_vert_is_from_mesh_1 = !next_halfedge_is_from_mesh_1;
				last_added_vert_is_an_intersection = false;
				halfedge_leading_to_last_added_vertex_index = next_halfedge_index;
				halfedge_leading_to_last_added_vertex_is_from_mesh_1 = next_halfedge_is_from_mesh_1;
				last_added_vert_event_index = 1;
				last_added_vert_was_an_endpoint = false;
				return last_added_vert;
			}
			break;
			default:
				CGAL_assertion(false); return -1;
			}
		}
	}
}

std::vector<size_t> collect_face_starting_at_event(Polyhedron& mesh_1, Polyhedron& mesh_2,
	const IntersectionEventStore& intersection_events,

	bool event_is_from_mesh_1,
	size_t event_halfedge_index,
	size_t event_list_index,

	const std::unordered_map<size_t, size_t>& vert_index_mapping_for_merged_mesh,
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
	CGAL::Verbose_ostream& vout,
	bool go_the_other_way
) {
	std::vector<size_t> face{};
	
	// Determine the starting vertex of this face, i.e. the vertex index corresponding to the starting event.
	size_t starting_vertex = -1;
	std::vector<IntersectionEvent> event_list = (event_is_from_mesh_1 ? intersection_events.intersection_events_1 : intersection_events.intersection_events_2).at(event_halfedge_index);
	IntersectionEvent starting_event = event_list[event_list_index];
	if (starting_event.tag == IntersectionEvent::IET_START) { // IE_START
		size_t no_offset_starting_vertex = (event_is_from_mesh_1 ? v_index_1 : v_index_2)[get_halfedge((event_is_from_mesh_1 ? mesh_1 : mesh_2), event_halfedge_index)->opposite()->vertex()];
		starting_vertex = (event_is_from_mesh_1 ? 0 : V1) + no_offset_starting_vertex;
	}
	else if (starting_event.tag == IntersectionEvent::IET_EE) {
		size_t halfedge_index_1 = event_is_from_mesh_1 ? event_halfedge_index : starting_event.other_edge_index();
		size_t halfedge_index_2 = event_is_from_mesh_1 ? starting_event.other_edge_index() : event_halfedge_index;
		starting_vertex = V1 + V2 + read_intersections_data(mesh_1, mesh_2, halfedge_index_1, halfedge_index_2, h_index_1, h_index_2, intersections_data);
	}
	else if (starting_event.tag == IntersectionEvent::IET_EV) {
		starting_vertex = (event_is_from_mesh_1 ? V1 : 0) + starting_event.other_vert_index();
	}
	else {
		// This situation should be disallowed in the parent (caller) function.
		CGAL_assertion(false);
	}
	// Insert the starting_vertex into the face.
	face.push_back(normalize_vertex_index_based_on_overlaps(starting_vertex, vert_index_mapping_for_merged_mesh));

	// Determine the second vertex of this face. It will necessarily lie along the given halfedge.
	size_t second_vertex = -1;
	bool second_vertex_is_an_endpoint = false;
	if (event_list_index == event_list.size() - 1 || event_list[event_list_index + 1].tag == IntersectionEvent::IET_ENDS_ON_EDGE
		                                          || event_list[event_list_index + 1].tag == IntersectionEvent::IET_ENDS_ON_VERT) {
		// If this event is the last non-end event on this halfedge, then the next vertex is the endpoint of this halfedge.
		// Recall that endpoints are not explicitly stored as events. The ENDS_ON events are there as extras when face collection might not amount to following the ->next() pointer.
		second_vertex_is_an_endpoint = true;
		size_t no_offset_second_vertex = (event_is_from_mesh_1 ? v_index_1 : v_index_2)[get_halfedge((event_is_from_mesh_1 ? mesh_1 : mesh_2), event_halfedge_index)->vertex()];
		second_vertex = (event_is_from_mesh_1 ? 0 : V1) + no_offset_second_vertex;
	}
	else {
		IntersectionEvent successor_event = event_list[event_list_index + 1];
		if (successor_event.tag == IntersectionEvent::IET_EE) {
			size_t halfedge_index_1 = event_is_from_mesh_1 ? event_halfedge_index : successor_event.other_edge_index();
			size_t halfedge_index_2 = event_is_from_mesh_1 ? successor_event.other_edge_index() : event_halfedge_index;
			second_vertex = V1 + V2 + read_intersections_data(mesh_1, mesh_2, halfedge_index_1, halfedge_index_2, h_index_1, h_index_2, intersections_data);
		}
		else if (successor_event.tag == IntersectionEvent::IET_EV) {
			second_vertex = (event_is_from_mesh_1 ? V1 : 0) + successor_event.other_vert_index();
		}
		else {
			// The ENDS_ON cases dealt with in the previous if scope.
			// The START case should be impossible, as this is not the first event on event_list.
			CGAL_assertion(false);
		}
	}
	// Insert the second_vertex into the face.
	face.push_back(normalize_vertex_index_based_on_overlaps(second_vertex, vert_index_mapping_for_merged_mesh));
	CGAL_assertion(face.size() == 2);
	CGAL_assertion(face[0] != face[1]);


	// Enter the main loop...
	size_t last_added_vert = second_vertex;
	bool last_added_vert_is_from_mesh_1 = (second_vertex < V1); // this variable is to be meaningless when 'last_added_vert_is_an_intersection' is true
	bool last_added_vert_is_an_intersection = (second_vertex >= V1 + V2);
	size_t halfedge_leading_to_last_added_vertex = event_halfedge_index;
	bool halfedge_leading_to_last_added_vertex_is_from_mesh_1 = event_is_from_mesh_1;
	size_t last_added_vert_event_index = event_list_index + 1; // this variable is to be meaningless when 'last_added_vert_was_an_endpoint' is true
	bool last_added_vert_was_an_endpoint = second_vertex_is_an_endpoint;

	while (true) {
		size_t next_vert_of_face = collect_edge_starting_at_event(mesh_1, mesh_2, intersection_events,
			last_added_vert, last_added_vert_is_from_mesh_1, last_added_vert_is_an_intersection, halfedge_leading_to_last_added_vertex,
			halfedge_leading_to_last_added_vertex_is_from_mesh_1, last_added_vert_event_index, last_added_vert_was_an_endpoint,
			intersections_data, vertex_overlap_map_1, vertex_overlap_map_2, hedges_through_vertex_1, hedges_through_vertex_2, V1, V2, V3, v_index_1, h_index_1, v_index_2, h_index_2, vout,
			go_the_other_way);
		
		size_t next_vert_of_face_normalized = normalize_vertex_index_based_on_overlaps(next_vert_of_face, vert_index_mapping_for_merged_mesh);
		if (face[0] == next_vert_of_face_normalized) {
			// We have looped back to the beginning of the face.
			break;
		}
		// Debug.
		bool face_is_malformed = false;
		for (size_t i = 1; i < face.size(); i++) {
			// Assert that we haven't looped back to the 'middle' of the face.
			//CGAL_assertion(face[i] != next_vert_of_face_normalized);
			if (face[i] == next_vert_of_face_normalized) {
				face_is_malformed = true;
				std::cerr << "Malformed face using GO_THE_OTHER_WAY=" + std::to_string(go_the_other_way) + "... Vertices:";
				for (size_t malformed_face_v : face) {
					std::cerr << " " << malformed_face_v << ",";
				}
				std::cerr << " [" << next_vert_of_face_normalized << "]" << std::endl;
				return std::vector<size_t>{}; // signal a problem by returning an empty face
			}
		}
		face.push_back(next_vert_of_face_normalized);
	}
	
	// Debug.
	for (size_t i = 1; i < face.size(); i++) {
		// Assert that all vertex indices are in the right range.
		CGAL_assertion(face[i] < V1 + V2 + V3);
	}
	return face;
}

// Collect faces of the merged mesh - new implementation based on IntersectionEvents.
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
	CGAL::Verbose_ostream& vout
) {
	std::unordered_map<size_t, std::vector<NormalizedVertexList>> result_faces{};
	CGAL::Timer debug_timer{}; debug_timer.start();
	double debug_last_timer_time = debug_timer.time();

	// Collect all faces containing at least a piece of some halfedge of mesh_1.
	for (size_t h_1 = 0; h_1 < mesh_1.size_of_halfedges(); h_1++) {

		// Debug progress.
		for (size_t percent = 2; percent <= 100; percent+=2) {
			if (h_1 == percent * mesh_1.size_of_halfedges() / 100) {
				// Found the percentage that's done.
				vout << "generate_merged_faces_from_events ~" << (percent/2 < 10 ? "0" : "") << percent/2 << "% done...";
				double debug_last_percent_time = debug_timer.time() - debug_last_timer_time; // How long it took to process the last 1% of total work.
				double debug_seconds_left_d = (200 - percent) * debug_last_percent_time;
				size_t debug_seconds_left = (size_t)debug_seconds_left_d;
				size_t debug_minutes_left = debug_seconds_left / 60;
				debug_seconds_left = debug_seconds_left - 60 * debug_minutes_left;
				vout << " Est. time left: " + std::to_string(debug_minutes_left) + "m " + std::to_string(debug_seconds_left) + "s." << std::endl;
				debug_last_timer_time = debug_timer.time();
				break;
			}
		}
		// End of debug.

		CGAL_assertion(intersection_events.intersection_events_1.find(h_1) != intersection_events.intersection_events_1.end());
		CGAL_assertion(intersection_events.intersection_events_1.at(h_1).size() >= 1); // At least the IE_START event should be present.
		CGAL_assertion(intersection_events.intersection_events_1.at(h_1)[0].tag == 4); // Check that the first event is the IE_START event.
		std::vector<IntersectionEvent> events = intersection_events.intersection_events_1.at(h_1);
		size_t event_index = 0;
		for (; event_index < events.size() && events[event_index].tag != 2 && events[event_index].tag != 3; event_index++) {

			// Get a face containing (a piece of) h_1.
			bool go_the_other_way = false;
			std::vector<size_t> face = collect_face_starting_at_event(mesh_1, mesh_2, intersection_events,
				true, h_1, event_index,
				vert_index_mapping_for_merged_mesh,
				intersections_data, vertex_overlap_map_1, vertex_overlap_map_2, hedges_through_vertex_1, hedges_through_vertex_2,
				V1, V2, V3, v_index_1, h_index_1, v_index_2, h_index_2, vout, go_the_other_way);

			if (face.empty()) {
				// The face was malformed (i.e. it was a lasso-shape) with go_the_other_way=false. Let's try with go_the_other_way=true.
				std::cerr << "Attempting to collect a face using go_the_other_way=true, as we got a malformed face with 'false'." << std::endl;
				face = collect_face_starting_at_event(mesh_1, mesh_2, intersection_events,
					true, h_1, event_index,
					vert_index_mapping_for_merged_mesh,
					intersections_data, vertex_overlap_map_1, vertex_overlap_map_2, hedges_through_vertex_1, hedges_through_vertex_2,
					V1, V2, V3, v_index_1, h_index_1, v_index_2, h_index_2, vout, !go_the_other_way);
			}
			CGAL_assertion(!face.empty()); // If this assertion fails, the face came out malformed with both go_the_other_way=true and go_the_other_way=false. :(

			// Insert the face into the deduplication data structure.
			NormalizedVertexList norm_face{ face };
			insert_into_map<size_t, NormalizedVertexList>(result_faces, norm_face.hash(), norm_face);
		}

		// Debug. Make sure that after the first 'ENDS_AT' event, there were no events of other type than 'ENDS_AT'.
		for (; event_index < events.size(); event_index++) {
			CGAL_assertion(events[event_index].tag == 2 || events[event_index].tag == 3);
		}
	}
	// Collect all faces containing at least a piece of some halfedge of mesh_2. We deduplicate at the end.
	for (size_t h_2 = 0; h_2 < mesh_2.size_of_halfedges(); h_2++) {

		// Debug progress.
		for (size_t percent = 2; percent <= 100; percent += 2) {
			if (h_2 == percent * mesh_2.size_of_halfedges() / 100) {
				// Found the percentage that's done.
				vout << "generate_merged_faces_from_events ~" << 50 + percent / 2 << "% done...";
				double debug_last_percent_time = debug_timer.time() - debug_last_timer_time; // How long it took to process the last 1% of total work.
				double debug_seconds_left_d = (100 - percent) * debug_last_percent_time;
				size_t debug_seconds_left = (size_t)debug_seconds_left_d;
				size_t debug_minutes_left = debug_seconds_left / 60;
				debug_seconds_left = debug_seconds_left - 60 * debug_minutes_left;
				vout << " Est. time left: " + std::to_string(debug_minutes_left) + "m " + std::to_string(debug_seconds_left) + "s." << std::endl;
				debug_last_timer_time = debug_timer.time();
				break;
			}
		}
		// End of debug.

		CGAL_assertion(intersection_events.intersection_events_2.find(h_2) != intersection_events.intersection_events_2.end());
		CGAL_assertion(intersection_events.intersection_events_2.at(h_2).size() >= 1);
		CGAL_assertion(intersection_events.intersection_events_2.at(h_2)[0].tag == 4); // Check that the first event is the IE_START event.
		std::vector<IntersectionEvent> events = intersection_events.intersection_events_2.at(h_2);
		size_t event_index = 0;
		for (; event_index < events.size() && events[event_index].tag != 2 && events[event_index].tag != 3; event_index++) {

			bool go_the_other_way = false;
			// Get a face containing (a piece of) h_1.
			std::vector<size_t> face = collect_face_starting_at_event(mesh_1, mesh_2, intersection_events,
				false, h_2, event_index,
				vert_index_mapping_for_merged_mesh,
				intersections_data, vertex_overlap_map_1, vertex_overlap_map_2, hedges_through_vertex_1, hedges_through_vertex_2,
				V1, V2, V3, v_index_1, h_index_1, v_index_2, h_index_2, vout, go_the_other_way);

			if (face.empty()) {
				// The face was malformed (i.e. it was a lasso-shape) with go_the_other_way=false. Let's try with go_the_other_way=true.
				std::cerr << "Attempting to collect a face using go_the_other_way=true, as we got a malformed face with 'false'." << std::endl;
				face = collect_face_starting_at_event(mesh_1, mesh_2, intersection_events,
					false, h_2, event_index,
					vert_index_mapping_for_merged_mesh,
					intersections_data, vertex_overlap_map_1, vertex_overlap_map_2, hedges_through_vertex_1, hedges_through_vertex_2,
					V1, V2, V3, v_index_1, h_index_1, v_index_2, h_index_2, vout, !go_the_other_way);
			}
			CGAL_assertion(!face.empty()); // If this assertion fails, the face came out malformed with both go_the_other_way=true and go_the_other_way=false. :(

			// Insert the face into the deduplication data structure.
			NormalizedVertexList norm_face{ face };
			insert_into_map<size_t, NormalizedVertexList>(result_faces, norm_face.hash(), norm_face);
		}

		// Debug. Make sure that after the first 'ENDS_AT' event, there were no events of other type than 'ENDS_AT'.
		for (; event_index < events.size(); event_index++) {
			CGAL_assertion(events[event_index].tag == 2 || events[event_index].tag == 3);
		}
	}

	auto ret = deduplicate(result_faces, vout);
	return ret;
}

// Make the points of the merged_mesh consistent with the faces -- because some vertices may have overlapped,
// some vertices of mesh_2 should possibly be discarded -- otherwise we'd get loose isolated points hanging around.
// Example: suppose that mesh_1 has verts 0..3 and mesh_2 has verts 4..8 and there are overlaps 1=5, 3=7.
// Then the 'normalize' function has the mapping as follows: [0,1,2,3,4,5,6,7,8]->
//                                                         ->[0,1,2,3,4,1,5,3,6].
// And thus the vector of Points should skip the original points number 5 and 7 which will happen with the implementation below.
std::vector<Point> generate_merged_points(const std::vector<Point>& points, const std::unordered_map<size_t, size_t>& vert_index_mapping_for_merged_mesh) {
	std::vector<Point> res{};
	size_t j = 0;
	for (size_t i = 0; i < points.size(); i++) {
		size_t mapped_i = vert_index_mapping_for_merged_mesh.at(i);
		if (mapped_i == j) {
			res.push_back(points[i]);
			j++;
		}
	}
	return res;
}


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
	std::unordered_map<size_t /*HCI*/, std::vector<size_t /*HCI*/>>& intersections_of_edges_1,
	std::unordered_map<size_t /*HCI*/, std::vector<size_t /*HCI*/>>& intersections_of_edges_2,
	std::unordered_map<size_t /*HCI*/, std::unordered_map<size_t /*HCI*/, size_t /*into intersection_points*/>>& intersections_data,
	const size_t V1,
	const size_t V2,
	const size_t V3,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	CGAL::Verbose_ostream& vout)
{
	std::unordered_map<size_t, std::vector<NormalizedVertexList>> result_faces{};

	// Debug. 
	debug_print_hedges(mesh_1, mesh_2, v_index_1, h_index_1, v_index_1, h_index_2, V1, vout);

	for (auto kvp : intersections_of_edges_1) {
		// Both for the halfedge kvp.first of mesh_1, and its twin, collect all faces adjacent to (a piece of) that halfedge.
		for (size_t e_1 : std::vector<size_t>{ kvp.first, get_opposite(mesh_1, kvp.first, h_index_1) }) {
			std::vector<size_t> e_2s = kvp.second;
			// If we're doing the twin we need to adjust the e_2s intersection list:
			if (e_1 != kvp.first) {
				std::vector<size_t> e_2s_adjusted{};
				for (int reverse_it = (int)e_2s.size() - 1; reverse_it >= 0; reverse_it--) {
					// reverse list and take opposites
					e_2s_adjusted.push_back(get_opposite(mesh_2, e_2s[reverse_it], h_index_2));
				}
				// replace sneakily
				e_2s = e_2s_adjusted;
			}

			size_t v_1 = get_start_vert(mesh_1, e_1, v_index_1);

			for (size_t e_2 : e_2s) {
				std::vector<size_t> face{};
				face.push_back(0 + v_1);
				//size_t v = get_intersection(e_1, e_2, V1, V2, intersections_data);
				size_t v = read_intersections_data(mesh_1, mesh_2, e_1, e_2, h_index_1, h_index_2, intersections_data);
				size_t v_global = V1 + V2 + v;
				size_t prev_e = e_1;
				size_t e = e_2;
				bool e_in_mesh_1 = false;
				bool v_is_an_original_vertex = false; // i.e. not an intersection vertex
													  // while we don't complete the loop around the face


													  // Debug.
				std::cerr << std::endl << 0 + v_1 << "x(" << v_is_an_original_vertex << "," << e_in_mesh_1 << ")";

				while (v_global != 0 + v_1) {
					face.push_back(v_global);

					// Debug.
					std::cerr << " " << v_global << "x(" << v_is_an_original_vertex << "," << e_in_mesh_1 << ")";

					if (e_in_mesh_1) {
						std::pair<bool, size_t> p;
						if (!v_is_an_original_vertex) {
							size_t prev_e_twin = get_opposite(mesh_2, prev_e, h_index_2);
							p = next_intersection(mesh_1, mesh_2, e, h_index_1, v_index_1, h_index_2, prev_e, prev_e_twin, intersections_of_edges_1);
						}
						else
							p = first_intersection(mesh_1, mesh_2, e, h_index_1, v_index_1, h_index_2, intersections_of_edges_1);
						if (p.first) {
							size_t e_next = p.second;
							//v = intersections_data[e][e_next];
							v = read_intersections_data(mesh_1, mesh_2, e, e_next, h_index_1, h_index_2, intersections_data);
							v_global = V1 + V2 + v;
							prev_e = e;
							e = e_next;
							e_in_mesh_1 = !e_in_mesh_1;
							v_is_an_original_vertex = false;
						}
						else {
							v = p.second;
							v_global = 0 + v;
							prev_e = e;
							e = get_next(mesh_1, e, h_index_1);
							v_is_an_original_vertex = true;
						}
					}
					else {
						std::pair<bool, size_t> p;
						if (!v_is_an_original_vertex) {
							size_t prev_e_twin = get_opposite(mesh_1, prev_e, h_index_1);
							p = next_intersection(mesh_2, mesh_1, e, h_index_2, v_index_2, h_index_1, prev_e, prev_e_twin, intersections_of_edges_2);
						}
						else
							p = first_intersection(mesh_2, mesh_1, e, h_index_2, v_index_2, h_index_1, intersections_of_edges_2);
						if (p.first) {
							size_t e_next = p.second;
							//v = intersections_data[e_next][e];
							v = read_intersections_data(mesh_1, mesh_2, e_next, e, h_index_1, h_index_2, intersections_data);
							v_global = V1 + V2 + v;
							prev_e = e;
							e = e_next;
							e_in_mesh_1 = !e_in_mesh_1;
							v_is_an_original_vertex = false;
						}
						else {
							v = p.second;
							v_global = V1 + v;
							prev_e = e;
							e = get_next(mesh_2, e, h_index_2);
							v_is_an_original_vertex = true;
						}
					}
				} // We have completed the face.

				NormalizedVertexList norm_face{ face };
				insert_into_map<size_t, NormalizedVertexList>(result_faces, norm_face.hash(), norm_face);
				//std::cout << "FACE!" << std::endl;

				// Move on to the next segment of e_1.
				//v_1 = V1 + V2 + get_intersection(e_1, e_2, V1, V2, intersections_data);
				v_1 = V1 + V2 + read_intersections_data(mesh_1, mesh_2, e_1, e_2, h_index_1, h_index_2, intersections_data);
			}

			// Do the last segment of e_1.
			{
				std::vector<size_t> face{};
				face.push_back(0 + v_1);
				size_t v = get_end_vert(mesh_1, e_1, v_index_1);
				size_t v_global = 0 + v;
				size_t prev_e = e_1;
				size_t e = get_next(mesh_1, e_1, h_index_1);
				bool e_in_mesh_1 = true;
				bool v_is_an_original_vertex = true; // i.e. not an intersection vertex
													 // while we don't complete the loop around the face

				HCI h1 = get_halfedge(mesh_1, e_1);

				// Debug.
				//size_t v1a = v_index_1[h1->opposite()->vertex()];
				//size_t v1b = v_index_1[h1->vertex()];
				//vout << "e_1: (" << v1a << "->" << v1b << ") " << "[m1]" << std::endl;
				std::cerr << std::endl << 0 + v_1 << "(" << v_is_an_original_vertex << "," << e_in_mesh_1 << ")";

				while (v_global != 0 + v_1) {
					face.push_back(v_global);

					// Debug.
					std::cerr << " " << v_global << "(" << v_is_an_original_vertex << "," << e_in_mesh_1 << ")";

					if (e_in_mesh_1) {
						std::pair<bool, size_t> p;
						if (!v_is_an_original_vertex) {
							size_t prev_e_twin = get_opposite(mesh_2, prev_e, h_index_2);
							p = next_intersection(mesh_1, mesh_2, e, h_index_1, v_index_1, h_index_2, prev_e, prev_e_twin, intersections_of_edges_1);
						}
						else
							p = first_intersection(mesh_1, mesh_2, e, h_index_1, v_index_1, h_index_2, intersections_of_edges_1);
						if (p.first) {
							size_t e_next = p.second;
							//v = intersections_data[e][e_next];
							v = read_intersections_data(mesh_1, mesh_2, e, e_next, h_index_1, h_index_2, intersections_data);
							v_global = V1 + V2 + v;
							prev_e = e;
							e = e_next;
							e_in_mesh_1 = !e_in_mesh_1;
							v_is_an_original_vertex = false;
						}
						else {
							v = p.second;
							v_global = 0 + v;
							prev_e = e;
							e = get_next(mesh_1, e, h_index_1);
							v_is_an_original_vertex = true;
						}
					}
					else {
						std::pair<bool, size_t> p;
						if (!v_is_an_original_vertex) {
							size_t prev_e_twin = get_opposite(mesh_1, prev_e, h_index_1);
							p = next_intersection(mesh_2, mesh_1, e, h_index_2, v_index_2, h_index_1, prev_e, prev_e_twin, intersections_of_edges_2);
						}
						else
							p = first_intersection(mesh_2, mesh_1, e, h_index_2, v_index_2, h_index_1, intersections_of_edges_2);
						if (p.first) {
							size_t e_next = p.second;
							//v = intersections_data[e_next][e];
							v = read_intersections_data(mesh_1, mesh_2, e_next, e, h_index_1, h_index_2, intersections_data);
							v_global = V1 + V2 + v;
							prev_e = e;
							e = e_next;
							e_in_mesh_1 = !e_in_mesh_1;
							v_is_an_original_vertex = false;
						}
						else {
							v = p.second;
							v_global = V1 + v;
							prev_e = e;
							e = get_next(mesh_2, e, h_index_2);
							v_is_an_original_vertex = true;
						}
					}
				} // We have completed the face.
				  //result_faces.push_back(face);

				NormalizedVertexList norm_face{ face };
				insert_into_map<size_t, NormalizedVertexList>(result_faces, norm_face.hash(), norm_face);
				//std::cout << "FACE!" << std::endl;
			}
			// We are done with e_1.
		}

	}
	// We are now done with all faces which contain at least one edge that's a segment of some edge of mesh_1.
	// All that's left are faces made entirely out of edges of mesh_2. 
	// Those are exactly the original faces of mesh_2, such that none of the face's edges was intersected.

	// So...
	// Iterate through each face of mesh_2. For each of them, get all the edges making them up. Check the intersection info of each edge for emptiness.

	for (FCI fit = mesh_2.facets_begin(); fit != mesh_2.facets_end(); fit++) {
		HCI eit = fit->halfedge();
		HCI eit_start = fit->halfedge();
		bool include_this_face = true;
		do {
			if (halfedge_or_opposite_is_intersected(mesh_2, h_index_2[eit], h_index_2, intersections_of_edges_2)) {
				include_this_face = false;
				break;
			}
			eit = eit->next();
		} while (eit != eit_start);

		// the face's edges had no intersections, so we should include the face.
		if (include_this_face) {
			std::vector<size_t> face{};
			eit = fit->halfedge();
			eit_start = fit->halfedge();
			do {
				face.push_back(V1 + v_index_2[eit->vertex()]);
				eit = eit->next();
			} while (eit != eit_start);

			NormalizedVertexList norm_face{ face };
			insert_into_map<size_t, NormalizedVertexList>(result_faces, norm_face.hash(), norm_face);
			//std::cout << "FACE!" << std::endl;
		}
	}

	auto ret = deduplicate(result_faces, vout);
	return ret;
}