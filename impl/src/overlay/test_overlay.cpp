#include "test_overlay.h"

#include "../embedding/io.h"
#include "geometry_utils.h"
#include "mesh_merger.h"
#include "../common/constants.h"

using namespace std;

class BasicTest {
public:
	BasicTest(bool(*_run)(void), int _index, std::string _d) : run(_run), description(_d), index(_index) {}
	bool(*run)(void);
	std::string description;
	int index;

	bool execute() {
		std::cout << "\nRunning test " << index << "..." << std::endl;
		std::cout << description << std::endl;
		bool succ = run();
		if (succ) {
			std::cout << "test " << index << " succeeded!" << std::endl;
		}
		else {
			std::cout << "test " << index << " failed!" << std::endl;
		}
		return succ;
	}
};


bool test_1() {
	Vector v1{ -0.026382, 0.809377, -0.477577 };  v1 = v1 / sqrt(v1.squared_length());
	Vector v2{ 0.499957,  0.786233, -0.124121 };  v2 = v2 / sqrt(v2.squared_length());
	Vector u1{ 0.441978,  0.62212,  -0.551479 };  u1 = u1 / sqrt(u1.squared_length());
	Vector u2{ 0.028292,  0.936964, -0.060968 };  u2 = u2 / sqrt(u2.squared_length());
	Point p1 = CGAL::ORIGIN + v1;
	Point p2 = CGAL::ORIGIN + v2;
	Point q1 = CGAL::ORIGIN + u1;
	Point q2 = CGAL::ORIGIN + u2;

	Point inter;
	intersect(p1, p2, q1, q2, inter);

	Point target{ 0.267337, 0.850778, -0.288422 };
	std::cout << "inter = " << inter << std::endl;

	bool condition = (target - inter).squared_length() < 0.01;
	CGAL_assertion(condition);
	return condition;
}


bool test_1b() {
	Vector v1{ -0.809830, 0.580352, 0.085834 }; // v1 = v1 / sqrt(v1.squared_length());
	Vector v2{ -0.810061, 0.575684, 0.111309 }; // v2 = v2 / sqrt(v2.squared_length());
	Vector u1{ 0.094696, 0.510592, 0.854593 };  //u1 = u1 / sqrt(u1.squared_length());
	Vector u2{ 0.125951, 0.502869, 0.855137 };  //u2 = u2 / sqrt(u2.squared_length());
	Point p1 = CGAL::ORIGIN + v1;
	Point p2 = CGAL::ORIGIN + v2;
	Point q1 = CGAL::ORIGIN + u1;
	Point q2 = CGAL::ORIGIN + u2;

	Point inter;
	bool b = intersect(p1, p2, q1, q2, inter);

	Point target{ 0.267337, 0.850778, -0.288422 };
	if (b)
		std::cout << "inter = " << inter << std::endl;
	if (!b)
		std::cout << "no intersection" << std::endl;

	//bool condition = (target - inter).squared_length() < 0.01;
	//CGAL_assertion(condition);
	//return condition;
	return true;
}


bool test_2() {
	Point p1 = CGAL::ORIGIN + Vector{ 0.215608, 0.764863, -0.640975 };
	Point p2 = CGAL::ORIGIN + Vector{ 0.915364, -0.059952, 0.291854 };
	Point p3 = CGAL::ORIGIN + Vector{ -0.341627, -0.796450, -0.509966 };
	Point p4 = CGAL::ORIGIN + Vector{ -0.664379, 0.297928, 0.699801 };

	Point q1 = CGAL::ORIGIN + Vector{ 0.013548, 0.527524, 0.826625 };
	Point q2 = CGAL::ORIGIN + Vector{ 0.400053, -0.792088, 0.440120 };
	Point q3 = CGAL::ORIGIN + Vector{ 0.695644, 0.396849, -0.684347 };
	Point q4 = CGAL::ORIGIN + Vector{ -0.860919, -0.132282, -0.434346 };

	std::vector<Point> ps{ p1, p2, p3, p4 };
	std::vector<Point> qs{ q1, q2, q3, q4 };

	size_t counter = 0;
	bool inconsistent = false;
	Point dummy;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int i2 = i+1; i2 < 4; i2++) {
				for (int j2 = j + 1; j2 < 4; j2++) {
					Point p = ps[i];
					Point p2 = ps[i2];
					Point q = qs[j];
					Point q2 = qs[j2];
					if (intersect(p, p2, q, q2, dummy)) {
						counter++;
					}
					if (intersect(p, p2, q, q2, dummy)) {
						if (!intersect(p2, p, q, q2, dummy)) {
							std::cout << "Inconsistency." << endl;
							inconsistent = true;
						}
						else if (!intersect(p2, p, q2, q, dummy)) {
							std::cout << "Inconsistency." << endl;
							inconsistent = true;
						}
						else if (!intersect(p, p2, q2, q, dummy)) {
							std::cout << "Inconsistency." << endl;
							inconsistent = true;
						}
					}
					else if (!intersect(p, p2, q, q2, dummy)) {
						if (intersect(p2, p, q, q2, dummy)) {
							std::cout << "Inconsistency." << endl;
							inconsistent = true;
						}
						else if (intersect(p2, p, q2, q, dummy)) {
							std::cout << "Inconsistency." << endl;
							inconsistent = true;
						}
						else if (intersect(p, p2, q2, q, dummy)) {
							std::cout << "Inconsistency." << endl;
							inconsistent = true;
						}
					}
				}
			}
		}
	}

	std::cout << "Number of distinct intersections: " << counter << endl;
	bool condition = counter == 6 && !inconsistent;
	CGAL_assertion(condition);
	return condition;
}


bool test_3() {
	Point p = CGAL::ORIGIN + Vector{ 0.915364, -0.059952, 0.291854 };

	Point q1 = CGAL::ORIGIN + Vector{ 0.013548, 0.527524, 0.826625 };
	Point q2 = CGAL::ORIGIN + Vector{ 0.400053, -0.792088, 0.440120 };
	Point q3 = CGAL::ORIGIN + Vector{ 0.695644, 0.396849, -0.684347 };
	Point q4 = CGAL::ORIGIN + Vector{ -0.860919, -0.132282, -0.434346 };
	std::vector<Point> qs{ q1, q2, q3, q4 };
	std::vector<std::vector<size_t>> q_faces{
		std::vector<size_t>{0, 1, 2},
		std::vector<size_t>{1, 0, 3},
		std::vector<size_t>{0, 2, 3},
		std::vector<size_t>{2, 1, 3},
	};
	Polyhedron poly_q;
	make_mesh_from_points_and_faces(qs, q_faces, poly_q);

	FInvIndex ind{ poly_q.facets_begin(), poly_q.facets_end() };
	size_t returned = containing_face_arbitrary(p, poly_q, ind);

	std::vector<size_t> face = q_faces[returned];
	std::cout << "Face from 'containing_face': [ind=" << returned << "], verts: ";
	for (size_t v_ind : face) {
		std::cout << v_ind << " ";
	}
	std::cout << endl;

	bool condition = true;
	for (size_t v_ind : face) {
		if (v_ind == 3) {
			condition = false;
			break;
		}
	}
	CGAL_assertion(condition);
	return condition;
}


bool test_4() {
	bool condition = true;
	auto f1 = NormalizedVertexList{ std::vector<size_t>{1,2,3,4} };
	auto f2 = NormalizedVertexList{ std::vector<size_t>{4,1,2,3} };
	auto g1 = NormalizedVertexList{ std::vector<size_t>{4,3,2,1} };
	auto g2 = NormalizedVertexList{ std::vector<size_t>{2,1,4,3} };

	std::cout << "Hash of {1,2,3,4}: " << f1.hash() << std::endl;
	std::cout << "Hash of {4,1,2,3}: " << f2.hash() << std::endl;
	std::cout << "Hash of {4,3,2,1}: " << g1.hash() << std::endl;
	std::cout << "Hash of {2,1,4,3}: " << g2.hash() << std::endl;

	condition = condition && (f1.equals(f2));
	CGAL_assertion(condition);
	condition = condition && (f2.equals(f1));
	CGAL_assertion(condition);
	condition = condition && (g1.equals(g2));
	CGAL_assertion(condition);
	condition = condition && (g2.equals(g1));
	CGAL_assertion(condition);
	condition = condition && (!f1.equals(g1));
	CGAL_assertion(condition);
	condition = condition && (!g1.equals(f1));
	CGAL_assertion(condition);
	condition = condition && (f1.hash() == f2.hash());
	CGAL_assertion(condition);
	condition = condition && (g1.hash() == g2.hash());
	CGAL_assertion(condition);
	condition = condition && (f1.hash() != g2.hash());
	CGAL_assertion(condition);
	condition = condition && (f2.hash() != g1.hash());
	CGAL_assertion(condition);
	return condition;
}


bool test_5() {
	Vector a{ 0.5130507, 0.6076790, 0.5424445 };
	Vector b{ 0.4480082, 0.4601987, 0.1618462 };
	Vector c{ 0.7646843, 0.9640288, 0.4111764 };
	double u = 0.1863127;
	double v = 0.2826308;
	double w = 1 - u - v;

	Vector p = u * a + v * b + w * c;

	std::vector<double> coords = barycentric_coords_planar(p, a, b, c);
	bool condition = true;

	condition = condition && (abs(coords[0] - u) < 0.001);
	std::cout << "coord[0] error: " << coords[0] - u << std::endl;
	CGAL_assertion(condition);
	condition = condition && (abs(coords[1] - v) < 0.001);
	std::cout << "coord[1] error: " << coords[1] - v << std::endl;
	CGAL_assertion(condition);
	condition = condition && (abs(coords[2] - w) < 0.001);
	std::cout << "coord[2] error: " << coords[2] - w << std::endl;
	CGAL_assertion(condition);
	return condition;
}


bool test_6() {
	Vector a{ 0.5130507, 0.6076790, 0.5424445 };
	Vector b{ 0.4480082, 0.4601987, 0.1618462 };
	Vector c{ 0.7646843, 0.9640288, 0.4111764 };
	double u = 0.1863127;
	double v = 0.2826308;
	double w = 1 - u - v;

	Vector p = u * a + v * b + w * c;

	p *= 10.2354;

	std::vector<double> coords = barycentric_coords_spherical(p, a, b, c);
	bool condition = true;

	condition = condition && (abs(coords[0] - u) < 0.001);
	std::cout << "coord[0] error: " << coords[0] - u << std::endl;
	CGAL_assertion(condition);
	condition = condition && (abs(coords[1] - v) < 0.001);
	std::cout << "coord[1] error: " << coords[1] - v << std::endl;
	CGAL_assertion(condition);
	condition = condition && (abs(coords[2] - w) < 0.001);
	std::cout << "coord[2] error: " << coords[2] - w << std::endl;
	CGAL_assertion(condition);
	return condition;
}


bool test_7() {
	bool success = true;

	// With a 'star' of n directions, test the function n times.
	
	Vector normal{0.75344, 0.00000, 0.65751};
	std::vector<Vector> neighbors{ // in order!
		Vector{0.83822, -0.33459, 0.41145},
		Vector{0.67557, -0.27627, 0.66528},
		Vector{0.52573, 0.0, 0.85065},
		Vector{0.63819, 0.26286, 0.72361},
		Vector{0.83105, 0.23885, 0.50230},
		Vector{0.87161, -0.15220, 0.49307}
	};
	std::vector<std::pair<Vector, size_t>> all_dirs{};
	size_t i = 0;
	for (Vector nei : neighbors) {
		all_dirs.push_back(std::make_pair(nei - normal, i++));
	}

	for (size_t i = 0; i < all_dirs.size(); i++) {
		Vector& reference_dir = all_dirs[i].first;

		std::vector<std::pair<Vector, size_t>> candidate_dirs{};
		for (size_t j = 0; j < all_dirs.size(); j++) {
			if (j != i)
				candidate_dirs.push_back(all_dirs[j]);
		}
		size_t correct_answer = (i + 1) % all_dirs.size();
		bool GO_THE_OTHER_WAY = false;
		size_t answer = get_neighbouring_halfedge_by_angle(normal, candidate_dirs, reference_dir, GO_THE_OTHER_WAY);

		if (correct_answer != answer) {
			success = false;
			std::cerr << "Wrong answer: " << "[" << answer << "] (" << all_dirs[answer].first << ") instead of: [" << correct_answer << "] (" << all_dirs[correct_answer].first << ")" << std::endl;
		}
	}
	return success;
}


void run_all_tests() {

	std::vector<BasicTest> tests{
		BasicTest(test_1, 1,
		"Check a single call to 'intersection'."),

		BasicTest(test_1b, 1,
		"Check a single call to 'intersection'."),

		BasicTest(test_2, 2,
		"Check the number of intersections found for two tetrahedrons."),

		BasicTest(test_3, 3,
		"Check a single 'containing_face'."),

		BasicTest(test_4, 4,
		"Check the NormalizedVertexList hashing."),

		BasicTest(test_5, 5,
		"Check the barycentric_coords_planar code."),

		BasicTest(test_6, 6,
		"Check the barycentric_coords_spherical code."),

		BasicTest(test_7, 7,
		"Check the get_neighbouring_halfedge_by_angle code."),
	};

	size_t failed = 0;

	for (BasicTest t : tests) {
		if (!t.execute()) {
			failed++;
		}
	}

	if (failed == 0) {
		std::cout << "All " << tests.size() << " tests successful." << std::endl;
	}
	else {
		std::cout << failed << "/" << tests.size() << " tests failed..." << std::endl;
	}
	ask_before_termination(DISABLE_USER_INTERACTION);
}