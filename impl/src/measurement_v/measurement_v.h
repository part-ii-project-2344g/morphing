#pragma once
#include "../common/typedefs.h"
#include "../embedding/io.h"

bool measure_v(double& out_result, size_t n_samples, size_t frames, std::string base_path, CGAL::Verbose_ostream& vout);
bool measure_v(double& out_result, size_t n_samples, size_t frames, std::string base_path, CGAL::Verbose_ostream& vout, const std::string& output_path);

struct ProcessedPolyhedron {
public:
	std::vector<Point> points;
	std::vector<std::vector<size_t>> faces;
	std::vector<std::vector<double>> faces_yz; // [0] - y low, [1] - y high, [2] - z low, [3] - z high
	ProcessedPolyhedron(const Polyhedron& p) {
		get_points(points, p);
		get_faces(faces, p);
		for (size_t i = 0; i < faces.size(); i++) {
			double y_l = CGAL_IA_MAX_DOUBLE;
			double y_h = CGAL_IA_MIN_DOUBLE;
			double z_l = CGAL_IA_MAX_DOUBLE;
			double z_h = CGAL_IA_MIN_DOUBLE;
			for (size_t v_ind : faces[i]) {
				Point pt = points[v_ind];
				y_l = std::min(y_l, pt.y());
				y_h = std::max(y_h, pt.y());
				z_l = std::min(z_l, pt.z());
				z_h = std::max(z_h, pt.z());
			}
			faces_yz.push_back(std::vector<double>{y_l, y_h, z_l, z_h});
		}
	}
};

struct AABB {
public:
	double x_l, x_h, y_l, y_h, z_l, z_h;
	double get_volume() {
		if (x_h <= x_l || y_h <= y_l || z_h <= z_l) return 0.0;
		return (x_h - x_l)*(y_h - y_l)*(z_h - z_l);
	}
	AABB() {
		x_l = CGAL_IA_MAX_DOUBLE;
		y_l = CGAL_IA_MAX_DOUBLE;
		z_l = CGAL_IA_MAX_DOUBLE;
		x_h = CGAL_IA_MIN_DOUBLE;
		y_h = CGAL_IA_MIN_DOUBLE;
		z_h = CGAL_IA_MIN_DOUBLE;
	}
	void add_mesh(const ProcessedPolyhedron& p) {
		for (size_t i = 0; i < p.points.size(); i++) {
			Point pt = p.points[i];
			double x = pt.x();
			double y = pt.y();
			double z = pt.z();
			x_l = std::min(x_l, x);
			y_l = std::min(y_l, y);
			z_l = std::min(z_l, z);
			x_h = std::max(x_h, x);
			y_h = std::max(y_h, y);
			z_h = std::max(z_h, z);
		}
	}
};

// Exposed for testing only.
double sym_diff_volume(const ProcessedPolyhedron& p1, const ProcessedPolyhedron& p2, size_t n_samples);
