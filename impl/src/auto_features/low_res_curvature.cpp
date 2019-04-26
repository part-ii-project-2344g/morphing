#include "low_res_curvature.h"

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

#include "../measurement/measurement.h"

namespace SMS = CGAL::Surface_mesh_simplification;

double weight(double dist_sq, double radius, double sqrt_radius) {
	double dist = sqrt(dist_sq);
	if (dist > radius)
		return 0.0;
	
	double sqrt_dist = sqrt(dist);
	return sqrt_radius - sqrt_dist; // 0.0 at dist=radius, sqrt_radius at dist=0.0, falls off quickly at the beginning
}

double low_res_curvature_at_vert(const VCI& vert, const std::vector<std::pair<Point, double>>& low_res_curvatures, const double radius, const double sqrt_radius, CGAL::Verbose_ostream& vout) {
	Point p = vert->point();
	double res = 0.0;
	double total_w = 0.0;
	double closest_curv = 0.0;
	double closest_distSq = CGAL_IA_MAX_DOUBLE;
	for (auto pair : low_res_curvatures) {
		double distSq = (pair.first - p).squared_length();
		if (distSq < closest_distSq) {
			closest_distSq = distSq;
			closest_curv = pair.second; // Used in case no contributions with nonzero weight are found.
		}
		double w = weight(distSq, radius, sqrt_radius);
		if (w > 0.0) {
			res += w * pair.second;
			total_w += w;
		}	
	}
	if (total_w == 0.0) {
		// Counted no contributions.
		return closest_curv;
	}
	return res / total_w;
}

Polyhedron low_res_mesh(const Polyhedron& mesh, size_t n_edges, CGAL::Verbose_ostream& vout) {
	// This is a stop predicate (defines when the algorithm terminates).
	// In this example, the simplification stops when the number of undirected edges
	// left in the surface mesh drops below the specified number (1000)
	SMS::Count_stop_predicate<Polyhedron> stop(n_edges);

	Polyhedron result{ mesh };

	// This the actual call to the simplification algorithm.
	// The surface mesh and stop conditions are mandatory arguments.
	// The index maps are needed because the vertices and edges
	// of this surface mesh lack an "id()" field.
	int r = SMS::edge_collapse(result,
		stop,
		CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, result))
		.halfedge_index_map(get(CGAL::halfedge_external_index, result))
	);

	vout << "\nFinished decimation, " << r << " edges removed. Final edge count: " << result.size_of_halfedges() / 2 << " whereas n_edges was " << n_edges << "." << std::endl;

	return result;
}

// Extract the curvature information from a low_res mesh, that will be passed to the 'low_res_curvature_at_vert' method, so we don't have to recompute it all the time.
std::vector<std::pair<Point, double>> get_low_res_curvatures(const Polyhedron& low_res) {
	std::vector<std::pair<Point, double>> res{};
	for (VCI it = low_res.vertices_begin(); it != low_res.vertices_end(); ++it) {
		res.push_back(std::make_pair(it->point(), curvature_at_vert(low_res, it)));
	}
	return res;
}