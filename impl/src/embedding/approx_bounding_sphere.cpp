#include "approx_bounding_sphere.h"

using namespace std;


Sphere compute_initial_sphere(std::vector<Point>& points) {
	double minx, miny, minz, maxx, maxy, maxz;
	minx = miny = minz = std::numeric_limits<double>::max();
	maxx = maxy = maxz = std::numeric_limits<double>::min();
	Point* pts[6]{};
	for (Point p : points) {
		minx = min(minx, p.x());
		if (minx == p.x()) pts[0] = &p;
		miny = min(miny, p.y());
		if (miny == p.y()) pts[1] = &p;
		minz = min(minz, p.z());
		if (minz == p.z()) pts[2] = &p;
		maxx = max(minx, p.x());
		if (maxx == p.x()) pts[3] = &p;
		maxy = max(maxy, p.y());
		if (maxy == p.y()) pts[4] = &p;
		maxz = max(maxz, p.z());
		if (maxz == p.z()) pts[5] = &p;
	}

	Sphere res{};
	double dist_x = ((*(pts[0])) - (*(pts[3]))).squared_length();
	double dist_y = ((*(pts[1])) - (*(pts[4]))).squared_length();
	double dist_z = ((*(pts[2])) - (*(pts[5]))).squared_length();
	double max_dist = max(max(dist_x, dist_y), dist_z);
	res.radius = sqrt(max_dist) / 2.0;
	if (dist_x == max_dist) {
		res.centre = *(pts[3]) + (*(pts[0]) - *(pts[3])) / 2.0;
	}
	else if (dist_y == max_dist) {
		res.centre = *(pts[4]) + (*(pts[1]) - *(pts[4])) / 2.0;
	}
	else if (dist_z == max_dist) {
		res.centre = *(pts[5]) + (*(pts[2]) - *(pts[5])) / 2.0;
	}

	return res;
}

Sphere compute_approx_bounding_sphere(std::vector<Point>& points) {
	Sphere sph = compute_initial_sphere(points);

	for (Point p : points) {
		double dx = p.x() - sph.centre.x();
		double dy = p.y() - sph.centre.y();
		double dz = p.z() - sph.centre.z();
		double old_to_p_sq = dx * dx + dy * dy + dz * dz;

		if (old_to_p_sq > sph.radius*sph.radius) {
			// Point p is outside current sphere sph. Update sph.
			double old_to_p = sqrt(old_to_p_sq);
			sph.radius = (sph.radius + old_to_p) / 2.0;
			double old_to_new = old_to_p - sph.radius;
			double new_x = (sph.radius * sph.centre.x() + old_to_new * p.x()) / old_to_p;
			double new_y = (sph.radius * sph.centre.y() + old_to_new * p.y()) / old_to_p;
			double new_z = (sph.radius * sph.centre.z() + old_to_new * p.z()) / old_to_p;
			sph.centre = Point{ new_x, new_y, new_z };
		}
	}

	return sph;
}

void normalize_bounding_sphere(std::vector<Point>& points, Sphere sph) {
	for (int i = 0; i < points.size(); i++) {
		Point p = points[i];
		double x = (p.x() - sph.centre.x()) / sph.radius;
		double y = (p.y() - sph.centre.y()) / sph.radius;
		double z = (p.z() - sph.centre.z()) / sph.radius;
		Point new_p{ x,y,z };
		points[i] = new_p;
	}
}