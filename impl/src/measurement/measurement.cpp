#include"measurement.h"
#include "../embedding/io.h"
#include"../interpolation/interpolate.h" // for 'generate_path'
#include"../overlay/cpp_utils.h" // set_contains
#include <fstream>
#include <iostream>
#include <ostream>
#include <unordered_set>

bool is_angle_degenerate(const Vector& v, const Vector& a, const Vector& b) {
	Vector va = a - v;
	Vector vb = b - v;
	return sqrt((va*va)*(vb*vb)) == 0;
}

double angle(const Vector& v, const Vector& a, const Vector& b) {
	Vector va = a - v;
	Vector vb = b - v;
	double cos_angle = (va)*(vb) / sqrt((va*va)*(vb*vb));

	if (std::abs(cos_angle) > 1.0 || cos_angle != std::max(-1.0, std::min(1.0, cos_angle))) {
		std::cerr << "cos_angle outside [-1,1]: " << cos_angle << " will be clamped." << std::endl;
		CGAL_assertion(std::abs(cos_angle) < 1.1);
	}

	cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
	return std::acos( cos_angle ); // TODO Verify formula
}

double curvature_at_vert(const Polyhedron& mesh, const VCI& vert) {
	double angle_deficit = 2 * CGAL_PI;
	HVCC around_vert_circulator = vert->vertex_begin();
	HVCC end = around_vert_circulator;
	Vector this_vert = vert->point() - CGAL::ORIGIN;
	Vector prev_nei = vert->vertex_begin()->opposite()->vertex()->point() - CGAL::ORIGIN;
	do {
		++around_vert_circulator;
		Vector next_nei = around_vert_circulator->opposite()->vertex()->point() - CGAL::ORIGIN;

		// Compute an angle adjacent to 'vert'.
		double angle_at_vert = angle(this_vert, prev_nei, next_nei);
		angle_deficit -= angle_at_vert;
		prev_nei = next_nei;
	} while (around_vert_circulator != end);
	return angle_deficit;
}

// O(|V|)
double curvature_at_vert(const Polyhedron& mesh, const VCI& vert, double averaging_radius) {
	double curvature = 0.0;
	size_t verts_counted = 0;
	Point p = vert->point();
	double dist_sq = averaging_radius * averaging_radius;
	for (VCI v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v) {
		if ((v->point() - p).squared_length() <= dist_sq) {
			verts_counted++;
			curvature += curvature_at_vert(mesh, v);
		}
	}
	return curvature / (double)verts_counted;
}

double measure_step(const Polyhedron& mesh_1, const Polyhedron& mesh_2, const std::unordered_set<size_t>& verts_to_exclude, CGAL::Verbose_ostream& vout) {
	double result = 0.0;

	CGAL_assertion(mesh_1.size_of_vertices() == mesh_2.size_of_vertices());
	size_t V = mesh_1.size_of_vertices();

	size_t debug_strongest_i = 0;
	double debug_strongest_diff = 0.0;
	double debug_strongest_curv_1 = 0.0;
	double debug_strongest_curv_2 = 0.0;

	VCI vertices_1 = mesh_1.vertices_begin();
	VCI vertices_2 = mesh_2.vertices_begin();
	for (size_t i = 0; i < V; i++, vertices_1++, vertices_2++) {
		if (set_contains<size_t>(verts_to_exclude, i))
			continue;
		double curv_1 = curvature_at_vert(mesh_1, vertices_1);
		double curv_2 = curvature_at_vert(mesh_2, vertices_2);
		double diff = abs(curv_1 - curv_2);
		result += diff;

		// Debug.
		if (diff > debug_strongest_diff) {
			debug_strongest_diff = diff;
			debug_strongest_curv_1 = curv_1;
			debug_strongest_curv_2 = curv_2;
			debug_strongest_i = i;
		}
	}

	// Debug.
	vout << "Vertex " << debug_strongest_i << " moved from (" << get_vertex(mesh_1, debug_strongest_i)->point() << ") to (" << get_vertex(mesh_2, debug_strongest_i)->point() << ") and changed curvature from " << debug_strongest_curv_1 << " to " << debug_strongest_curv_2 << " by a maximal " << debug_strongest_curv_2 - debug_strongest_curv_1 << std::endl;

	return result;
}

std::unordered_set<size_t> verts_to_exclude_from_measurement_helper(const Polyhedron& mesh) {
	std::unordered_set<size_t> res{};
	VCI vert = mesh.vertices_begin();
	for (size_t i = 0; i < mesh.size_of_vertices(); i++, ++vert) {
		HVCC around_vert_circulator = vert->vertex_begin();
		HVCC end = around_vert_circulator;
		Vector this_vert = vert->point() - CGAL::ORIGIN;
		Vector prev_nei = vert->vertex_begin()->opposite()->vertex()->point() - CGAL::ORIGIN;
		do {
			++around_vert_circulator;
			Vector next_nei = around_vert_circulator->opposite()->vertex()->point() - CGAL::ORIGIN;

			// Compute an angle adjacent to 'vert'.
			bool degeneracy = is_angle_degenerate(this_vert, prev_nei, next_nei);
			if (degeneracy) {
				res.insert(i);
				break;
			}
		} while (around_vert_circulator != end);
	}
	return res;
}

std::unordered_set<size_t> verts_to_exclude_from_measurement(const std::vector<Polyhedron>& mesh_frames, CGAL::Verbose_ostream& vout) {
	std::unordered_set<size_t> res{};
	size_t debug_counter = 0;
	for (const Polyhedron& mesh : mesh_frames) {
		auto verts = verts_to_exclude_from_measurement_helper(mesh);
		for (size_t v : verts)
			res.insert(v);
		vout << "verts_to_exclude_from_measurement done with " << ++debug_counter << "/" << mesh_frames.size() << " frames." << std::endl;
	}
	return res;
}

// Returns true iff the function ran successfully.
bool measure(double& out_result, size_t frames, std::string base_path, CGAL::Verbose_ostream& vout) {
	
	std::vector<Polyhedron> mesh_frames{};
	for (size_t i = 0; i < frames; i++) {
		std::string path_i = generate_path(base_path, i);
		std::vector<Point> points_tmp;
		std::vector<std::vector<size_t>> faces_tmp;
		Polyhedron mesh;
		if (!input_mesh(points_tmp, faces_tmp, mesh, vout, path_i)) {
			std::cerr << "Error reading input meshes." << std::endl;
			return false;
		} // mesh loaded in
		mesh_frames.push_back(mesh);
	}

	std::unordered_set<size_t> verts_to_exclude = verts_to_exclude_from_measurement(mesh_frames, vout);

	// Debug.
	vout << "verts_to_exclude.size = " << verts_to_exclude.size() << std::endl;

	double start_to_end_difference = measure_step(mesh_frames[0], mesh_frames[frames-1], verts_to_exclude, vout);
	double total_difference = 0.0;

	// Debug.
	vout << "Difference between the extremal frames: " << start_to_end_difference << std::endl;

	for (size_t i = 1; i < frames; i++) {
		double this_step = measure_step(mesh_frames[i-1], mesh_frames[i], verts_to_exclude, vout);
		total_difference += this_step;

		// Debug.
		vout << "Difference between frames #" << i << " and #" << i+1 << ": " << this_step << std::endl;
	}

	double result = total_difference / start_to_end_difference;

	// Debug.
	vout << "\n\nFinal measurement result: " << total_difference << "/" << start_to_end_difference << " = " << result << "." << std::endl;
	
	out_result = result;
	return true;
}

bool measure(double& out_result, size_t frames, std::string base_path, CGAL::Verbose_ostream& vout, const std::string& output_path) {

	std::vector<Polyhedron> mesh_frames{};
	for (size_t i = 0; i < frames; i++) {
		std::string path_i = generate_path(base_path, i);
		std::vector<Point> points_tmp;
		std::vector<std::vector<size_t>> faces_tmp;
		Polyhedron mesh;
		if (!input_mesh(points_tmp, faces_tmp, mesh, vout, path_i)) {
			std::cerr << "Error reading input meshes." << std::endl;
			return false;
		} // mesh loaded in
		mesh_frames.push_back(mesh);
	}

	std::unordered_set<size_t> verts_to_exclude = verts_to_exclude_from_measurement(mesh_frames, vout);

	// Debug.
	vout << "verts_to_exclude.size = " << verts_to_exclude.size() << std::endl;

	double start_to_end_difference = measure_step(mesh_frames[0], mesh_frames[frames - 1], verts_to_exclude, vout);
	double total_difference = 0.0;

	// Debug.
	vout << "Difference between the extremal frames: " << start_to_end_difference << std::endl;

	for (size_t i = 1; i < frames; i++) {
		double this_step = measure_step(mesh_frames[i - 1], mesh_frames[i], verts_to_exclude, vout);
		total_difference += this_step;

		// Debug.
		vout << "Difference between frames #" << i << " and #" << i + 1 << ": " << this_step << std::endl;
	}

	double result = total_difference / start_to_end_difference;
	out_result = result;

	std::ostream* p_out = &std::cout;
	std::ofstream out;
	out.open(output_path);
	p_out = &out;
	if (!*p_out) {
		std::cerr << "ERROR: cannot open file for writing: " << output_path << std::endl;
		return false;
	}

	vout << "\nWriting a measurement result to " << output_path.substr(output_path.find_last_of("\\/") + 1) << " ..." << std::endl;
	*p_out << "Final measurement result: " << total_difference << "/" << start_to_end_difference << " = " << result << "." << std::endl;
	vout << "\n\nFinal measurement result: " << total_difference << "/" << start_to_end_difference << " = " << result << "." << std::endl;
	
	out.close();
	return true;
}