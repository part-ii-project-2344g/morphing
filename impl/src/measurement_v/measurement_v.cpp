#include "measurement_v.h"

#include"../interpolation/interpolate.h" // for 'generate_path'

#include <fstream>
#include <iostream>
#include <ostream>
#include <unordered_set>


AABB get_bounding_box(const ProcessedPolyhedron& p1, const ProcessedPolyhedron& p2) {
	AABB bb{};
	bb.add_mesh(p1);
	bb.add_mesh(p2);
	return bb;
}

std::unordered_set<size_t> get_faces_through_yz(double y, double z, const ProcessedPolyhedron& p) {
	std::unordered_set<size_t> res{};
	for (size_t i = 0; i < p.faces_yz.size(); i++) {
		double y_l = p.faces_yz[i][0];
		double y_h = p.faces_yz[i][1];
		double z_l = p.faces_yz[i][2];
		double z_h = p.faces_yz[i][3];
		if (y_l <= y && y <= y_h && z_l <= z && z <= z_h) {
			res.insert(i);
		}
	}
	return res;
}

// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm#cite_note-1
bool ray_triangle_intersection(const Point& ray_origin, const Vector& ray_dir, const Point& a, const Point& b, const Point& c) {
	const double EPSILON = 0.0000001;
	Vector edge1, edge2, h, s, q;
	double aa, f, u, v;
	edge1 = b - a;
	edge2 = c - a;
	h = CGAL::cross_product(ray_dir, edge2);
	aa = edge1 * h;
	if (aa > -EPSILON && aa < EPSILON)
		return false;    // This ray is parallel to this triangle.
	f = 1.0 / aa;
	s = ray_origin - a;
	u = f * (s * h);
	if (u < 0.0 || u > 1.0)
		return false;
	q = CGAL::cross_product(s, edge1);
	v = f * (ray_dir * q);
	if (v < 0.0 || u + v > 1.0)
		return false;
	// At this stage we can compute t to find out where the intersection point is on the line.
	double t = f * (edge2 * q);
	if (t > EPSILON) {
		return true;
	}
	else // This means that there is a line intersection but not a ray intersection.
		return false;
}

// ray straight in the positive x direction, starts at pt.
bool ray_face_intersection(const Point& pt, const ProcessedPolyhedron& p, size_t face_ind) {
	CGAL_assertion(p.faces[face_ind].size() == 3);
	Point a = (p.points[p.faces[face_ind][0]]);
	Point b = (p.points[p.faces[face_ind][1]]);
	Point c = (p.points[p.faces[face_ind][2]]);
	return ray_triangle_intersection(pt, Vector{ 1.0, 0.0, 0.0 }, a, b, c);
}

bool is_point_inside_mesh(const Point& pt, const ProcessedPolyhedron& p) {
	double y = pt.y();
	double z = pt.z();

	std::unordered_set<size_t> faces_through_yz = get_faces_through_yz(y, z, p);
	
	size_t hits = 0;

	for(size_t f_ind : faces_through_yz){
		if (ray_face_intersection(pt, p, f_ind)) {
			hits++;
		}
	}
	
	return (hits % 2) == 1;
}

Point sample_uniform_point_in_bb(const AABB& bb, int seed, int seq) {
	double x01 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	double y01 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	double z01 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	double x = bb.x_l + x01 * (bb.x_h - bb.x_l);
	double y = bb.y_l + y01 * (bb.y_h - bb.y_l);
	double z = bb.z_l + z01 * (bb.z_h - bb.z_l);
	return Point{x, y, z};
}

#define SEED 0
double sym_diff_volume(const ProcessedPolyhedron& p1, const ProcessedPolyhedron& p2, size_t n_samples) {
	srand(static_cast <unsigned> (SEED));
	AABB bb = get_bounding_box(p1, p2);
	size_t hits = 0;
	for (size_t i = 0; i < n_samples; i++) {
		Point pt = sample_uniform_point_in_bb(bb, (int)SEED, (int)i);
		bool pt_inside_symmetric_difference = is_point_inside_mesh(pt, p1) ^ is_point_inside_mesh(pt, p2);
		if (pt_inside_symmetric_difference)
			hits++;
	}
	double hit_percentage = ((double)hits) / ((double)n_samples);
	return bb.get_volume() * hit_percentage;
}


// Returns true iff the function ran successfully.
bool measure_v(double& out_result, size_t n_samples, size_t frames, std::string base_path, CGAL::Verbose_ostream& vout) {

	std::vector<ProcessedPolyhedron> mesh_frames{};
	for (size_t i = 0; i < frames; i++) {
		std::string path_i = generate_path(base_path, i);
		std::vector<Point> points_tmp;
		std::vector<std::vector<size_t>> faces_tmp;
		Polyhedron mesh;
		if (!input_mesh(points_tmp, faces_tmp, mesh, vout, path_i)) {
			std::cerr << "Error reading input meshes." << std::endl;
			return false;
		} // mesh loaded in
		mesh_frames.push_back(ProcessedPolyhedron{ mesh });
	}

	double start_to_end_difference = sym_diff_volume(mesh_frames[0], mesh_frames[frames - 1], n_samples);
	double total_difference = 0.0;

	// Debug.
	vout << "Volume difference between the extremal frames: " << start_to_end_difference << std::endl;

	for (size_t i = 1; i < frames; i++) {
		double this_step = sym_diff_volume(mesh_frames[i - 1], mesh_frames[i], n_samples);
		total_difference += this_step;

		// Debug.
		vout << "Volume difference between frames #" << i << " and #" << i + 1 << ": " << this_step << std::endl;
	}

	double result = total_difference / start_to_end_difference;

	// Debug.
	vout << "\n\nFinal measurement_v result: " << total_difference << "/" << start_to_end_difference << " = " << result << "." << std::endl;

	out_result = result;
	return true;
}


// Returns true iff the function ran successfully.
bool measure_v(double& out_result, size_t n_samples, size_t frames, std::string base_path, CGAL::Verbose_ostream& vout, const std::string& output_path) {

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

	double start_to_end_difference = sym_diff_volume(mesh_frames[0], mesh_frames[frames - 1], n_samples);
	double total_difference = 0.0;

	// Debug.
	vout << "Volume difference between the extremal frames: " << start_to_end_difference << std::endl;

	for (size_t i = 1; i < frames; i++) {
		double this_step = sym_diff_volume(mesh_frames[i - 1], mesh_frames[i], n_samples);
		total_difference += this_step;

		// Debug.
		vout << "Volume difference between frames #" << i << " and #" << i + 1 << ": " << this_step << std::endl;
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

	vout << "\nWriting a measurement_v result to " << output_path.substr(output_path.find_last_of("\\/") + 1) << " ..." << std::endl;
	*p_out << "Final measurement_v result: " << total_difference << "/" << start_to_end_difference << " = " << result << "." << std::endl;
	vout << "\n\nFinal measurement_v result: " << total_difference << "/" << start_to_end_difference << " = " << result << "." << std::endl;

	out.close();
	return true;
}