#include "geometry_utils.h"
#include <CGAL/basic.h>

#include "../common/constants.h"
#include "../embedding/geometry_utils.h"
#include "../overlay/cpp_utils.h"

using namespace std;

// TODO: Get the threshold right.
bool identify_vertices(Vector v1, Vector v2) {
	return sqrt((v1 - v2).squared_length()) < VERTICES_IDENTIFICATION_THRESHOLD;
}

bool intersect_h(const Point& p1, const Point& p2, const Point& q1, const Point& q2, Point& output);

// returns true if p1--p2 and q1--q2 intersect and the result was written to the output variable
bool intersect(const Point& p1, const Point& p2, const Point& q1, const Point& q2, Point& output) {
	// We don't want a situation in which we'd say that intersect(e,f) succeeds while intersect(f,e) fails, etc.
	// That would break the code in 'intersect_asymmetric' which uses next_interesection().
	// To protect ourselves we do this:
	if (intersect_h(p1, p2, q1, q2, output))
		return true;
	if (intersect_h(p1, p2, q2, q1, output))
		return true;
	if (intersect_h(p2, p1, q1, q2, output))
		return true;
	if (intersect_h(p2, p1, q2, q1, output))
		return true;
	if (intersect_h(q1, q2, p1, p2, output))
		return true;
	if (intersect_h(q1, q2, p2, p1, output))
		return true;
	if (intersect_h(q2, q1, p1, p2, output))
		return true;
	if (intersect_h(q2, q1, p2, p1, output))
		return true;
	return false;
}

// wrapped this in intersect() for robustness
bool intersect_h(const Point& p1, const Point& p2, const Point& q1, const Point& q2, Point& output) {
	Vector v1 = p1 - CGAL::ORIGIN;
	Vector v2 = p2 - CGAL::ORIGIN;
	Vector u1 = q1 - CGAL::ORIGIN;
	Vector u2 = q2 - CGAL::ORIGIN;

	Vector r = CGAL::cross_product(CGAL::cross_product(v1, v2), CGAL::cross_product(u1, u2));
	r = normalized(r); // to make values later on easier to understand

	double p1x = p1.x();
	double p2x = p2.x();
	double p1y = p1.y();
	double p2y = p2.y();
	double p1z = p1.z();
	double p2z = p2.z();

	double q1x = q1.x();
	double q2x = q2.x();
	double q1y = q1.y();
	double q2y = q2.y();
	double q1z = q1.z();
	double q2z = q2.z();

	double rx = r.x();
	double ry = r.y();
	double rz = r.z();

	double sp = 0, sq = 0;

	bool succ_p_positive = false;
	bool succ_p_negative = false;
	bool succ_q_positive = false;
	bool succ_q_negative = false;


	if (rx == 0.0 && ry == 0.0 && rz == 0.0) { // r is the zero vector iff p1p2 and q1q2 are parallel.
											   // in that case the arcs overlap or not.
											   // if they do, their midpoints conincide and we return that.
											   // To achieve that, we set r to point to that point.
											   // If the arcs do not intersect, the code will return false as is appropriate.
		r = (v1 + v2) / 2.0;
		r = normalized(r);
		rx = r.x();
		ry = r.y();
		rz = r.z();

		// TODO: what if that r happens to be zero too (i.e. v1 and v2 are diametric)
	}

	{
		// We're dividing by p2x-p1x (or similar), we want to choose one that is as large as possible to minimize loss of precision.
		double diff_x = abs(p2x - p1x);
		double diff_y = abs(p2y - p1y);
		double diff_z = abs(p2z - p1z);
		double diff_max = max(max(diff_x, diff_y), diff_z);
		bool all_diffs_zero = (diff_x == 0.0 && diff_y == 0.0 && diff_z == 0.0);

		if (!all_diffs_zero && diff_x == diff_max) {
			double nominator = p1y * p2x - p1x * p2y;
			double denominator = ry * (p2x - p1x) - rx * (p2y - p1y);
			double tp = nominator / denominator;
			sp = (tp*rx - p1x) / (p2x - p1x);

			// Debug.
			if (abs(tp * rx - (p1x + sp * (p2x - p1x))) > 0.01) {
				cerr << "1Large error " << abs(tp * rx - (p1x + sp * (p2x - p1x))) << " in intersect(), due to small denom " << p2x - p1x << endl;
			}
			if (abs(tp * ry - (p1y + sp * (p2y - p1y))) > 0.01) {
				cerr << "2Large error " << abs(tp * ry - (p1y + sp * (p2y - p1y))) << " in intersect(), due to small denom " << p2y - p1y << endl;
			}

			double error = tp * rz - (p1z + sp * (p2z - p1z));
			if (abs(error) > 0.0001) // TODO choose a good value for catching logic errors
				return false;

			if (tp > 0) succ_p_positive = true;
			else        succ_p_negative = true;
		}
		else if (!all_diffs_zero && diff_y == diff_max) {
			double nominator = p1x * p2y - p1y * p2x;
			double denominator = rx * (p2y - p1y) - ry * (p2x - p1x);
			double tp = nominator / denominator;
			sp = (tp*ry - p1y) / (p2y - p1y);

			// Debug.
			if (abs(tp * ry - (p1y + sp * (p2y - p1y))) > 0.01) {
				cerr << "3Large error " << abs(tp * ry - (p1y + sp * (p2y - p1y))) << " in intersect(), due to small denom " << p2y - p1y << endl;
			}
			if (abs(tp * rx - (p1x + sp * (p2x - p1x))) > 0.01) {
				cerr << "4Large error " << abs(tp * rx - (p1x + sp * (p2x - p1x))) << " in intersect(), due to small denom " << p2x - p1x << endl;
			}

			double error = tp * rz - (p1z + sp * (p2z - p1z));
			if (abs(error) > 0.0001) // TODO choose a good value for catching logic errors
				return false;

			if (tp > 0) succ_p_positive = true;
			else        succ_p_negative = true;
		}
		else if (!all_diffs_zero && diff_z == diff_max) {
			double nominator = p1x * p2z - p1z * p2x;
			double denominator = rx * (p2z - p1z) - rz * (p2x - p1x);
			double tp = nominator / denominator;
			sp = (tp*rz - p1z) / (p2z - p1z);

			// Debug.
			if (abs(tp * rz - (p1z + sp * (p2z - p1z))) > 0.01) {
				cerr << "5Large error " << abs(tp * rz - (p1z + sp * (p2z - p1z))) << " in intersect(), due to small denom " << p2z - p1z << endl;
			}
			if (abs(tp * rx - (p1x + sp * (p2x - p1x))) > 0.01) {
				cerr << "6Large error " << abs(tp * rx - (p1x + sp * (p2x - p1x))) << " in intersect(), due to small denom " << p2x - p1x << endl;
			}

			double error = tp * ry - (p1y + sp * (p2y - p1y));
			if (abs(error) > 0.0001) // TODO choose a good value for catching logic errors
				return false;

			if (tp > 0) succ_p_positive = true;
			else        succ_p_negative = true;
		}
		else { // all_diffs_zero, i.e. p1 == p2
			double tp;
			if (rx != 0.0)
				tp = p1x / rx;
			else if (ry != 0.0)
				tp = p1y / ry;
			else if (rz != 0.0)
				tp = p1z / rz;
			else // that case should have been caught at the beginning
				CGAL_assertion(false);

			if (tp > 0) succ_p_positive = true;
			else        succ_p_negative = true;
		}
	}
	{
		// We're dividing by q2x-q1x (or similar), we want to choose one that is as large as qossible to minimize loss of qrecision.
		double diff_x = abs(q2x - q1x);
		double diff_y = abs(q2y - q1y);
		double diff_z = abs(q2z - q1z);
		double diff_max = max(max(diff_x, diff_y), diff_z);
		bool all_diffs_zero = (diff_x == 0.0 && diff_y == 0.0 && diff_z == 0.0);

		if (!all_diffs_zero && diff_x == diff_max) {
			double nominator = q1y * q2x - q1x * q2y;
			double denominator = ry * (q2x - q1x) - rx * (q2y - q1y);
			double tq = nominator / denominator;
			sq = (tq*rx - q1x) / (q2x - q1x);

			// Debug.
			if (abs(tq * rx - (q1x + sq * (q2x - q1x))) > 0.01) {
				cerr << "7Large error " << abs(tq * rx - (q1x + sq * (q2x - q1x))) << " in intersect(), due to small denom " << q2x - q1x << endl;
			}
			if (abs(tq * ry - (q1y + sq * (q2y - q1y))) > 0.01) {
				cerr << "8Large error " << abs(tq * ry - (q1y + sq * (q2y - q1y))) << " in intersect(), due to small denom " << q2y - q1y << endl;
			}

			double error = tq * rz - (q1z + sq * (q2z - q1z));
			if (abs(error) > 0.0001) // TODO choose a good value for catching logic errors
				return false;

			if (tq > 0) succ_q_positive = true;
			else        succ_q_negative = true;
		}
		else if (!all_diffs_zero && diff_y == diff_max) {
			double nominator = q1x * q2y - q1y * q2x;
			double denominator = rx * (q2y - q1y) - ry * (q2x - q1x);
			double tq = nominator / denominator;
			sq = (tq*ry - q1y) / (q2y - q1y);

			// Debug.
			if (abs(tq * ry - (q1y + sq * (q2y - q1y))) > 0.01) {
				cerr << "9Large error " << abs(tq * ry - (q1y + sq * (q2y - q1y))) << " in intersect(), due to small denom " << q2y - q1y << endl;
			}
			if (abs(tq * rx - (q1x + sq * (q2x - q1x))) > 0.01) {
				cerr << "10Large error " << abs(tq * rx - (q1x + sq * (q2x - q1x))) << " in intersect(), due to small denom " << q2x - q1x << endl;
			}

			double error = tq * rz - (q1z + sq * (q2z - q1z));
			if (abs(error) > 0.0001) // TODO choose a good value for catching logic errors
				return false;

			if (tq > 0) succ_q_positive = true;
			else        succ_q_negative = true;
		}
		else if (!all_diffs_zero && diff_z == diff_max) {
			double nominator = q1x * q2z - q1z * q2x;
			double denominator = rx * (q2z - q1z) - rz * (q2x - q1x);
			double tq = nominator / denominator;
			sq = (tq*rz - q1z) / (q2z - q1z);

			// Debug.
			if (abs(tq * rz - (q1z + sq * (q2z - q1z))) > 0.01) {
				cerr << "11Large error " << abs(tq * rz - (q1z + sq * (q2z - q1z))) << " in intersect(), due to small denom " << q2z - q1z << endl;
			}
			if (abs(tq * rx - (q1x + sq * (q2x - q1x))) > 0.01) {
				cerr << "12Large error " << abs(tq * rx - (q1x + sq * (q2x - q1x))) << " in intersect(), due to small denom " << q2x - q1x << endl;
			}

			double error = tq * ry - (q1y + sq * (q2y - q1y));
			if (abs(error) > 0.0001) // TODO choose a good value for catching logic errors
				return false;

			if (tq > 0) succ_q_positive = true;
			else        succ_q_negative = true;
		}
		else { // all_diffs_zero, i.e. q1 == q2
			double tq;
			if (rx != 0.0)
				tq = q1x / rx;
			else if (ry != 0.0)
				tq = q1y / ry;
			else if (rz != 0.0)
				tq = q1z / rz;
			else // that case should have been caught at the beginning
				CGAL_assertion(false);

			if (tq > 0) succ_q_positive = true;
			else        succ_q_negative = true;
		}
	}

	if (sp < 0 || sp > 1 || sq < 0 || sq > 1)
		return false;

	// if tp and tq had different signs, then the shortest arcs do not intersect, so we return false
	if (succ_p_negative && succ_q_negative)
		output = CGAL::ORIGIN - normalized(r);
	else if (succ_p_positive && succ_q_positive)
		output = CGAL::ORIGIN + normalized(r);
	else {
		return false;
	}
	return true;
}

bool intersect(const HCI& h1, const HCI& h2, Point& output) {
	Point p1 = h1->vertex()->point();
	Point p2 = h1->opposite()->vertex()->point();
	Point q1 = h2->vertex()->point();
	Point q2 = h2->opposite()->vertex()->point();
	return intersect(p1, p2, q1, q2, output);
}

double point_segment_distance(const Vector& v, const Vector& s1, const Vector& s2) {
	Vector seg = s2 - s1;
	Vector s1v = v - s1;
	Vector s2v = v - s2;

	if (seg * s1v <= 0.0)
		return sqrt(s1v.squared_length());
	if (seg * s2v >= 0.0)
		return sqrt(s2v.squared_length());
	return sqrt(CGAL::cross_product(seg, s1v).squared_length() / seg.squared_length());
}

bool does_vert_edge_intersection_qualify(double dist) {
	return dist < VERTEX_EDGE_INTERSECTION_THRESHOLD; // TODO: Make this work.
}

bool vert_edge_intersection(const Point& vert, const Point& edge_start, const Point& edge_end) {
	double dist = point_segment_distance(vert - CGAL::ORIGIN, edge_start - CGAL::ORIGIN, edge_end - CGAL::ORIGIN);
	if (does_vert_edge_intersection_qualify(dist))
		return true;
	double dist_2 = point_segment_distance(vert - CGAL::ORIGIN, edge_end - CGAL::ORIGIN, edge_start - CGAL::ORIGIN); // for robustness
	if (does_vert_edge_intersection_qualify(dist_2))
		return true;
	return false;
}

bool vert_edge_intersection(const Point & vert, const HCI & edge) {
	return vert_edge_intersection(vert, edge->opposite()->vertex()->point(), edge->vertex()->point());
}

// Does fi contain v?
bool containing_face_helper(const FCI& fi, const Vector& v) {
	Vector a = fi->halfedge()->vertex()->point() - CGAL::ORIGIN;
	Vector b = fi->halfedge()->next()->vertex()->point() - CGAL::ORIGIN;
	Vector c = fi->halfedge()->next()->next()->vertex()->point() - CGAL::ORIGIN;

	Point pa = fi->halfedge()->vertex()->point();
	Point pb = fi->halfedge()->next()->vertex()->point();
	Point pc = fi->halfedge()->next()->next()->vertex()->point();
	Point pv = CGAL::ORIGIN + v;

	if (vert_edge_intersection(pv, pa, pb) || vert_edge_intersection(pv, pb, pc) || vert_edge_intersection(pv, pc, pa)) {
		return true;
	}

	if (identify_vertices(a, v) || identify_vertices(b, v) || identify_vertices(c, v)) {
		return true;
	}

	return CGAL::cross_product(a, b) * v >= 0.0 && CGAL::cross_product(b, c) * v >= 0.0 && CGAL::cross_product(c, a) * v >= 0.0;
}

// Returns the face of 'mesh' that contains (spherically) 'v'.
// Does not support multiple results.
size_t containing_face(const Point& p, const Vector& offset, const Polyhedron& mesh, const FInvIndex& f_index) {
	Vector v = p - CGAL::ORIGIN;
	return containing_face(v, offset, mesh, f_index);
}

// Returns the face of 'mesh' that contains (spherically) 'v'.
// Returns an arbitrary one in case of multiple hits.
size_t containing_face_arbitrary(const Point& p, const Polyhedron& mesh, const FInvIndex& f_index) {
	Vector v = p - CGAL::ORIGIN;
	return containing_face_arbitrary(v, mesh, f_index);
}


// Returns the face of 'mesh' that contains (spherically) 'v'.
// Returns an arbitrary one in case of multiple hits.
size_t containing_face_arbitrary(const Vector& v, const Polyhedron& mesh, const FInvIndex& f_index) {
	size_t debug_attempts = 0;
	std::vector<size_t> solutions{};

	// Somehow it seems that a==b==c for >600 faces of the deer embedding... 
	// I think they are just too small for the precision we're using. Blender does show that neighbouring verts have the same exact position. That is problematic.
	// Yes, the deer embedding has many verts with coords -0.819375 -0.109848 0.562635
	// But they do look different in Blender, weird. Does Blender cheat and move them around slightly based on connectivity? 
	// I should try to write my own OBJ writer with higher precision.
	for (FCI fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi) {
		debug_attempts++;
		if (containing_face_helper(fi, v)) {
			solutions.push_back(f_index[fi]);
		}
	}
	
	if (solutions.size() > 0)
		return solutions[0];

	cerr << "containing_face_arbitrary() failed to find a containing face after attempting (all) " << debug_attempts << " faces." << endl;
	CGAL_assertion(false);
	return -1;
}

// Returns the face of 'mesh' that contains (spherically) 'v'.
// Does not support multiple results (containing faces).
size_t containing_face(const Vector& v, const Vector& offset, const Polyhedron& mesh, const FInvIndex& f_index) {
	size_t debug_attempts = 0;
	std::vector<size_t> solutions{};

	for (FCI fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi) {
		debug_attempts++;
		if (containing_face_helper(fi, v)) {
			solutions.push_back(f_index[fi]);
		}
	}
	bool found_2 = false;
	size_t found_2_sol = -1;
	double found_2_eps = 0.0;
	if (solutions.size() >= 1) {
		if (solutions.size() >= 2) {
			cerr << "containing_face() found " << solutions.size() << " different faces. Will use an epsilon offset to find the correct face." << endl;

			// in the case of conflict, containing_face should return the face containing point+epsilon*dir
			// (what if that is still ambiguous between two faces? - Hopefully the intersection collection algorithm should handle that.)

			std::vector<size_t> refined_solutions{};
			for (double epsilon = 1e-16; epsilon < 1; epsilon *= 10.0) {
				refined_solutions.clear();
				for (size_t sol : solutions) {
					FCI candidate_face = get_facet(mesh, sol);
					Vector v_offset = v + epsilon * offset;
					if (containing_face_helper(candidate_face, v_offset)) {
						refined_solutions.push_back(sol);
					}
				}

				if (refined_solutions.size() == 1)
					return refined_solutions[0];
				if (refined_solutions.size() == 2 && !found_2) {
					found_2 = true;
					found_2_sol = refined_solutions[0];
					found_2_eps = epsilon;
				}
			}
			if (found_2) {
				cerr << "containing_face with largest epsilon found " << 2 << " faces " << "at epsilon=" << found_2_eps << ", probably because of an edge overlap, will return an arbitrary one of the two now..." << endl;
				return found_2_sol;
			}
			cerr << "containing_face() failed to find a unique containing face (or exactly two) after attempting a refined search. This is an error." << endl;
			CGAL_assertion(false);
		}
		return solutions[0];
	}

	cerr << "containing_face() failed to find a containing face after attempting (all) " << debug_attempts << " faces." << endl;
	CGAL_assertion(false);
	return -1;
}

std::vector<double> barycentric_coords_planar(Vector p, Vector a, Vector b, Vector c) {
	// gamedev.stackexchange.com/a/23745/76255

	Vector v0 = b - a, v1 = c - a, v2 = p - a;
	double d00 = v0 * v0;
	double d01 = v0 * v1;
	double d11 = v1 * v1;
	double d20 = v2 * v0;
	double d21 = v2 * v1;
	double denom = d00 * d11 - d01 * d01;
	double v = (d11 * d20 - d01 * d21) / denom;
	double w = (d00 * d21 - d01 * d20) / denom;
	double u = 1.0f - v - w;

	return std::vector<double>{u, v, w};
}

Vector plane_line_intersection(Vector planePoint, Vector planeNormal, Vector linePoint, Vector lineDirection) {
	if (planeNormal * lineDirection == 0) {
		// TODO: Handle corner case.
		CGAL_assertion(false);
		return Vector{ 0,0,0 };
	}

	double t = (planeNormal*planePoint - planeNormal*linePoint) / (planeNormal*lineDirection);
	return linePoint + (lineDirection * t);
}

std::vector<double> barycentric_coords_spherical(Vector p, Vector a, Vector b, Vector c) {
	// gamedev.stackexchange.com/a/23745/76255

	// Cast p onto the abc plane along the vector from origin to p.
	Vector p_cast = plane_line_intersection(a, CGAL::cross_product(a - b, c - b), Vector{ 0,0,0 }, p);
	return barycentric_coords_planar(p_cast, a, b, c);
}

Vector get_edge_direction_normalized(const HCI& edge) {
	return normalized(edge->vertex()->point() - edge->opposite()->vertex()->point());
}

Vector line_projection(const Vector& a, const Vector&b, const Vector& p) {
	return a + (b - a)*(((p - a)*(b - a)) / ((b - a)*(b - a)));
}