#include "relaxation.h"
#include "geometry_utils.h"
#include "io.h"

using namespace std;

#define COLLAPSE_THRESHOLD 0.2

size_t debug_relaxation_counter = 1;
// Assumes that the ordering in points corresponds to ordering in mesh.vertices.
// Returns the largest movement of any point during this relaxation step.
double relaxation_step(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::cpp11::tuple<size_t, size_t, size_t, size_t> fixed_vertices,
	CGAL::Verbose_ostream& vout, std::string& input_filename, bool export_intermediate_mesh) {
	std::vector<Point> new_points{};

	Vertex_iterator vertex_it = mesh.vertices_begin();
	CGAL_assertion(mesh.size_of_vertices() == points.size());

	double largest_movement_sq = 0.0;

	for (int i = 0; i < points.size(); i++, vertex_it++) {
		Point p = vertex_it->point();
		CGAL_assertion(p == points[i]);
		if (i == CGAL::cpp11::get<0>(fixed_vertices) ||
			i == CGAL::cpp11::get<1>(fixed_vertices) ||
			i == CGAL::cpp11::get<2>(fixed_vertices) ||
			i == CGAL::cpp11::get<3>(fixed_vertices)) {
			new_points.push_back(normalized(p)); // The fixed points get normalized at first iteration. This was not described explicitly in the paper.
			continue;
		}

		HVCC hv = vertex_it->vertex_begin();
		HVCC hv_end = hv;

		// we can add up vectors, but this wouldn't be possible with points
		Vector sum_of_neighbours{ 0.0,0.0,0.0 };
		do {
			// point to vector conversion
			sum_of_neighbours += hv->opposite()->vertex()->point() - CGAL::ORIGIN;
			++hv;
		} while (hv != hv_end);

		sum_of_neighbours = normalized(sum_of_neighbours);
		// vector to point conversion
		new_points.push_back(CGAL::ORIGIN + sum_of_neighbours);
		largest_movement_sq = max(largest_movement_sq, (points[i] - new_points[i]).squared_length());
	}

	CGAL_assertion(vertex_it == mesh.vertices_end());
	CGAL_assertion(points.size() == new_points.size());
	points = new_points;
#define SLOW 0
#if SLOW
	make_mesh_from_points_and_faces(points, faces, mesh);
#else
	int i = 0;
	for (vertex_it = mesh.vertices_begin(); vertex_it != mesh.vertices_end(); vertex_it++, i++) {
		(*vertex_it).point() = new_points[i];
	}
#endif
	vout << "Relaxation step [" << debug_relaxation_counter << "] had largest movement: " << sqrt(largest_movement_sq) << endl;

	//debug
	if (export_intermediate_mesh) {
		std::string S = input_filename; // for brevity in next line
		std::string output_filename = S.substr(0, S.length() - 4) + "[after_" + std::to_string(debug_relaxation_counter) + "_relaxations]" + S.substr(S.find_last_of("."));
		output_mesh(output_filename, mesh, vout, false);
	}

	debug_relaxation_counter++;
	return sqrt(largest_movement_sq);
}

bool is_embedding_collapsed(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, 
	CGAL::cpp11::tuple<Point, Point, Point, Point> tetrahedron_vertices, CGAL::Verbose_ostream& vout) {
	double tetrahedron_edge_length = sqrt((CGAL::cpp11::get<0>(tetrahedron_vertices) - CGAL::cpp11::get<1>(tetrahedron_vertices)).squared_length());

	Point diametric_of_tetrahedron_a = diametric(CGAL::cpp11::get<0>(tetrahedron_vertices));
	Point diametric_of_tetrahedron_b = diametric(CGAL::cpp11::get<1>(tetrahedron_vertices));
	Point diametric_of_tetrahedron_c = diametric(CGAL::cpp11::get<2>(tetrahedron_vertices));
	Point diametric_of_tetrahedron_d = diametric(CGAL::cpp11::get<3>(tetrahedron_vertices));
	Point neighbour_a = points[find_nearest_neighbour_of_point(points, diametric_of_tetrahedron_a)];
	Point neighbour_b = points[find_nearest_neighbour_of_point(points, diametric_of_tetrahedron_b)];
	Point neighbour_c = points[find_nearest_neighbour_of_point(points, diametric_of_tetrahedron_c)];
	Point neighbour_d = points[find_nearest_neighbour_of_point(points, diametric_of_tetrahedron_d)];
	double dist_a_sq = (diametric_of_tetrahedron_a - neighbour_a).squared_length();
	double dist_b_sq = (diametric_of_tetrahedron_b - neighbour_b).squared_length();
	double dist_c_sq = (diametric_of_tetrahedron_c - neighbour_c).squared_length();
	double dist_d_sq = (diametric_of_tetrahedron_d - neighbour_d).squared_length();
	double dist = sqrt(max(max(max(dist_a_sq, dist_b_sq), dist_c_sq), dist_d_sq));

	vout << "COLLAPSE_THRESHOLD <" << dist / tetrahedron_edge_length << " would be necessary to classify this embedding as collapsed." << endl;
	return dist > COLLAPSE_THRESHOLD * tetrahedron_edge_length;
}