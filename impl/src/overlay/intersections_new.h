#pragma once

#include "../common/typedefs.h"
#include "cpp_utils.h"

#include <unordered_map>

struct IntersectionEventPointer {
public:
	bool valid; // is this pointer meant to be read
	bool pointing_to_mesh_1;
	size_t halfedge_index;
	size_t event_index_in_halfedge;
};
class IntersectionEvent {
public:
	static const int IET_EE = 0;
	static const int IET_EV = 1;
	static const int IET_ENDS_ON_EDGE = 2;
	static const int IET_ENDS_ON_VERT = 3;
	static const int IET_START = 4;
//	// True iff this is the last event in this halfedge.
//	bool is_last_event;
//
//	// Points to the next IntersectionEvent on this halfedge, if one exists.
//	IntersectionEventPointer next;

	// Points to the IntersectionEvent that should be examined next if you are currently collecting faces of the merged mesh.
	IntersectionEventPointer turn;

	// Used for ordering events in a single halfedge.
	Vector pos; 

	// tag
	size_t tag; // IE_EE = 0, IE_EV = 1, IE_ENDS_ON_EDGE = 2, IE_ENDS_ON_VERT = 3, IE_START = 4
protected:
	size_t other_mesh_index; // edge_index for IE_EE, IE_ENDS_ON_EDGE, vert_index for IE_EV, IE_ENDS_ON_VERT, unused for IE_START
private:
	std::string pt_to_string(const Point& p) {
		return "(" + std::to_string(p.x()) + ", " + std::to_string(p.y()) + ", " + std::to_string(p.z()) + ")";
	}
public:
	IntersectionEvent(Point pos_, size_t index_) : pos(pos_ - CGAL::ORIGIN), other_mesh_index(index_) {}
	size_t other_vert_index() { if (tag == 0 || tag == 2) CGAL_assertion(false); return other_mesh_index; }
	size_t other_edge_index() { if (tag == 1 || tag == 3) CGAL_assertion(false); return other_mesh_index; }

	std::string to_string() {
		std::string s1 = ((tag == 0) ? "EDGE-EDGE" : (tag == 1) ? "EDGE-VERT" : (tag == 2) ? "ENDS-EDGE" : (tag == 3) ? "ENDS-VERT" : "ERROR-TAG");
		std::string s2 = " event at location (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " + std::to_string(pos.z()) + ") with other mesh's ";
		std::string s3 = ((tag == 0 || tag == 2) ? "edge #" : "vert #") + std::to_string(other_mesh_index);
		return s1 + s2 + s3;
	}
	std::string to_string(Polyhedron& mesh_1, Polyhedron& mesh_2, bool this_is_mesh_1, size_t this_hedge_index, const HInvIndex& h_index_1, const HInvIndex& h_index_2) { 
		if (tag != IET_EE) return "this to_string() was not_implemented";
		// edge-edge case:
		Polyhedron& this_mesh = this_is_mesh_1 ? mesh_1 : mesh_2;
		Polyhedron& other_mesh = this_is_mesh_1 ? mesh_2 : mesh_1;
		const HInvIndex& this_h_index = this_is_mesh_1 ? h_index_1 : h_index_2;
		const HInvIndex& other_h_index = this_is_mesh_1 ? h_index_2 : h_index_1;
		HCI this_hedge = get_halfedge(this_mesh, this_hedge_index);
		HCI other_hedge = get_halfedge(other_mesh, this->other_mesh_index);
		size_t this_hedge_ind = this_h_index[this_hedge];
		size_t other_hedge_ind = other_h_index[other_hedge];
		CGAL_assertion(other_hedge_ind == this->other_mesh_index);
		Point u1 = this_hedge->opposite()->vertex()->point();
		Point u2 = this_hedge->vertex()->point();
		Point v1 = other_hedge->opposite()->vertex()->point();
		Point v2 = other_hedge->vertex()->point();
		std::string this_h_str = std::string("mesh_") + (this_is_mesh_1 ? "1" : "2") + ".h#" + std::to_string(this_hedge_ind) + "[" + pt_to_string(u1) + "->" + pt_to_string(u2) + "]";
		std::string other_h_str = std::string("mesh_") + (this_is_mesh_1 ? "2" : "1") + ".h#" + std::to_string(other_hedge_ind) + "[" + pt_to_string(v1) + "->" + pt_to_string(v2) + "]";
		return ((IntersectionEvent*)(this))->to_string() + "; This hedge: " + this_h_str + ", other hedge: " + other_h_str;
	}
};

// edge-edge IntersectionEvent (i.e. this halfedge intersects some halfedge of the other mesh here)
class IE_EE : public IntersectionEvent {
public:
	IE_EE(Point pos_, size_t other_edge_index_) : IntersectionEvent(pos_, other_edge_index_) { tag = IET_EE; }
};

// edge-vertex IntersectionEvent (i.e. this halfedge goes right through some vertex of the other mesh here)
class IE_EV : public IntersectionEvent {
public:
	IE_EV(Point pos_, size_t other_vert_index_) : IntersectionEvent(pos_, other_vert_index_) { tag = IET_EV; }
};

// There can be multiple END events for a single halfedge, all of them helping determine where to turn from its endpoint during face collection.
// endpoint IntersectionEvent (i.e. this halfedge ends here at an edge of the other mesh - depending on the angles we might want to turn into that other edge)
class IE_ENDS_ON_EDGE : public IntersectionEvent {
public:
	IE_ENDS_ON_EDGE(Point pos_, size_t other_edge_index_) : IntersectionEvent(pos_, other_edge_index_) { tag = IET_ENDS_ON_EDGE; }
};

// endpoint IntersectionEvent (i.e. this halfedge ends here at a vertex of _the other_ mesh - depending on the angles we might want to turn into an edge of the other mesh)
class IE_ENDS_ON_VERT : public IntersectionEvent {
public:
	IE_ENDS_ON_VERT(Point pos_, size_t other_vert_index_) : IntersectionEvent(pos_, other_vert_index_) { tag = IET_ENDS_ON_VERT; }
};

// endpoint IntersectionEvent (i.e. this halfedge ends here at a vertex of _the other_ mesh - depending on the angles we might want to turn into an edge of the other mesh)
class IE_START : public IntersectionEvent {
public:
	IE_START(Point pos_) : IntersectionEvent(pos_, 0) { tag = IET_START; }
};

// holds all intersection events for a pair of meshes
class IntersectionEventStore {
public:
	// Maps indices of halfedges of mesh_1 to ordered lists of IntersectionEvents.
	std::unordered_map<size_t, std::vector<IntersectionEvent>> intersection_events_1;
	// Maps indices of halfedges of mesh_2 to ordered lists of IntersectionEvents.
	std::unordered_map<size_t, std::vector<IntersectionEvent>> intersection_events_2;

	IntersectionEvent retrieve(IntersectionEventPointer ptr) {
		return (ptr.pointing_to_mesh_1 ? intersection_events_1 : intersection_events_2)[ptr.halfedge_index][ptr.event_index_in_halfedge];
	}
};


void generate_verts_of_containing_face(
	Polyhedron& mesh_1, Polyhedron& mesh_2,
	std::unordered_map<size_t /*mesh_1 verts*/, std::vector<size_t /*mesh_2 verts*/>>& out_verts_of_containing_face_1,
	std::unordered_map<size_t /*mesh_2 verts*/, std::vector<size_t /*mesh_1 verts*/>>& out_verts_of_containing_face_2,
	const VInvIndex& v_index_1,
	const HInvIndex& h_index_1,
	const FInvIndex& f_index_1,
	const VInvIndex& v_index_2,
	const HInvIndex& h_index_2,
	const FInvIndex& f_index_2);

void collect_intersection_events(Polyhedron& mesh_1, Polyhedron& mesh_2,
	IntersectionEventStore& out_intersection_events,
	std::unordered_map<size_t /*HCI*/, std::unordered_map<size_t /*HCI*/, size_t /*into intersection_points*/>>& out_intersections_data,
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
	CGAL::Verbose_ostream& vout);
