#include "simple_features.h"
#include "features_io.h"
#include "../overlay/cpp_utils.h" // get_vertex
#include "low_res_curvature.h" // low_res_curvature_at_vert
#include "../measurement/measurement.h" // curvature_at_vert
#include "../embedding/io.h" // mesh io for test function at the end

#include <vector>
#include <tuple>
#include <algorithm>


///////////////////////////////// 'poisson' mode /////////////////////////////////

#define LOW_RES_NUMBER_OF_EDGES 1110

void poisson_features_to_file(const Polyhedron& mesh_1, const Polyhedron& mesh_2, const Polyhedron& low_res_1, const Polyhedron& low_res_2, 
	double PERCENTAGE_OF_VERTS_TO_CHECK, size_t N_FEATURES, double POISSON_RADIUS,double DIST_WEIGHT, double AVERAGING_RADIUS, std::string output_path, 
	CGAL::Verbose_ostream& vout) {
	std::vector<std::pair<size_t, size_t>> features = get_poisson_features(mesh_1, mesh_2, low_res_1, low_res_2, PERCENTAGE_OF_VERTS_TO_CHECK, N_FEATURES, POISSON_RADIUS, DIST_WEIGHT, AVERAGING_RADIUS, vout);
	output_features(features, output_path, mesh_1, mesh_2);
}

double correspondence_strength(double curv_1, double curv_2, double dist, double dist_scale) {
	return (curv_1 + curv_2) - dist * dist_scale;
}

std::vector<std::pair<size_t, size_t>> get_poisson_features(const Polyhedron& mesh_1, const Polyhedron& mesh_2, const Polyhedron& low_res_1, const Polyhedron& low_res_2,
	double PERCENTAGE_OF_VERTS_TO_CHECK /* 0.1 */, size_t N_FEATURES /* 3 */, double POISSON_RADIUS /* 0.3 */, double DIST_WEIGHT /* 1.0 */, double AVERAGING_RADIUS /* 0.3 */,
	CGAL::Verbose_ostream& vout) {

	vout << "Getting poisson_features with params:";
	vout << "\nPERCENTAGE_OF_VERTS_TO_CHECK = " << PERCENTAGE_OF_VERTS_TO_CHECK;
	vout << "\nN_FEATURES = " << N_FEATURES;
	vout << "\nPOISSON_RADIUS = " << POISSON_RADIUS;
	vout << "\nDIST_WEIGHT = " << DIST_WEIGHT;
	vout << "\nAVERAGING_RADIUS = " << AVERAGING_RADIUS;
	vout << std::endl << std::endl;

	double SQRT_AVERAGING_RADIUS = sqrt(AVERAGING_RADIUS);
	auto low_res_curvatures_1 = get_low_res_curvatures(low_res_1);
	auto low_res_curvatures_2 = get_low_res_curvatures(low_res_2);

	VCI vertices_1 = mesh_1.vertices_begin();
	std::vector<std::tuple<double, Point, size_t>> features_1{};
	for (size_t i = 0; i < mesh_1.size_of_vertices(); i++, vertices_1++) {
		double curv_1 = low_res_curvature_at_vert(vertices_1, low_res_curvatures_1, AVERAGING_RADIUS, SQRT_AVERAGING_RADIUS, vout);
		Point pos = vertices_1->point();
		features_1.push_back(std::make_tuple(curv_1, pos, i));

		// Debug.
		for (size_t percent = 0; percent <= 100; percent += 10) {
			if (i == mesh_1.size_of_vertices() - 1) {
				vout << 100 << "% done with mesh_1 curvature computation..." << std::endl << std::endl;
				break;
			}
			if (i == percent * mesh_1.size_of_vertices() / 100) {
				vout << percent << "% done with mesh_1 curvature computation..." << std::endl;
				break;
			}
		}
	}
	VCI vertices_2 = mesh_2.vertices_begin();
	std::vector<std::tuple<double, Point, size_t>> features_2{};
	for (size_t i = 0; i < mesh_2.size_of_vertices(); i++, vertices_2++) {
		double curv_2 = low_res_curvature_at_vert(vertices_2, low_res_curvatures_2, AVERAGING_RADIUS, SQRT_AVERAGING_RADIUS, vout);
		Point pos = vertices_2->point();
		features_2.push_back(std::make_tuple(curv_2, pos, i));

		// Debug.
		for (size_t percent = 0; percent <= 100; percent += 10) {
			if (i == mesh_1.size_of_vertices() - 1) {
				vout << 100 << "% done with mesh_2 curvature computation..." << std::endl << std::endl;
				break;
			}
			if (i == percent * mesh_2.size_of_vertices() / 100) {
				vout << percent << "% done with mesh_2 curvature computation..." << std::endl;
				break;
			}
		}
	}
	std::sort(features_1.begin(), features_1.end());
	std::sort(features_2.begin(), features_2.end());

	// Debug.
	vout << "Top curvature verts of mesh_1:\n";
	for (size_t i = 0; i < 3; i++) {
		std::tuple<double, Point, size_t> tup = features_1[features_1.size() - 1 - i];
		vout << std::get<1>(tup) << " [v_" << std::get<2>(tup) << "]:  curv = " << std::get<0>(tup) << " = 2pi - " << (2 * CGAL_PI - std::get<0>(tup)) * (360.0 / (2.0 * CGAL_PI)) << "deg \n";
	}
	vout << "Top curvature verts of mesh_2:\n";
	for (size_t i = 0; i < 3; i++) {
		std::tuple<double, Point, size_t> tup = features_2[features_2.size() - 1 - i];
		vout << std::get<1>(tup) << " [v_" << std::get<2>(tup) << "]:  curv = " << std::get<0>(tup) << " = 2pi - " << (2 * CGAL_PI - std::get<0>(tup)) * (360.0 / (2.0 * CGAL_PI)) << "deg \n";
	}

	const size_t N = (int)(PERCENTAGE_OF_VERTS_TO_CHECK * (float)(std::min(features_1.size(), features_2.size())));
	double high_curvature = std::min(std::get<0>(features_1[features_1.size() - 1]), std::get<0>(features_2[features_2.size() - 1]));
	// Let's say we want a distance of 1.0 to count as much as the curv difference of 'high_curvature'
	const double DIST_SCALE = DIST_WEIGHT * high_curvature / 1.0;

	std::vector<std::tuple<double, size_t, size_t>> candidate_feature_pairs{};
	for (size_t i = features_1.size() - N; i < features_1.size(); i++) {
		for (size_t j = features_2.size() - N; j < features_2.size(); j++) {
			size_t feature_1_ind = std::get<2>(features_1[i]);
			size_t feature_2_ind = std::get<2>(features_2[j]);
			double feature_1_curv = std::get<0>(features_1[i]);
			double feature_2_curv = std::get<0>(features_2[j]);
			Point feature_1_pos = std::get<1>(features_1[i]);
			Point feature_2_pos = std::get<1>(features_2[j]);
			double dist = std::sqrt((feature_1_pos - feature_2_pos).squared_length());
			candidate_feature_pairs.push_back(std::make_tuple(correspondence_strength(feature_1_curv, feature_2_curv, dist, DIST_SCALE), feature_1_ind, feature_2_ind));
		}
	}
	std::sort(candidate_feature_pairs.begin(), candidate_feature_pairs.end());

	std::vector<std::pair<size_t, size_t>> result{};
	std::unordered_set<size_t> indices_used_1{};
	std::unordered_set<size_t> indices_used_2{};
	std::vector<Point> points_used_1{};
	std::vector<Point> points_used_2{};

	size_t j = candidate_feature_pairs.size();
	for (size_t i = 0; i < N_FEATURES; i++) {

		while (true) {
			if (j == 0) {
				std::cerr << "Failed to poisson-assign features." << std::endl;
				//CGAL_assertion(j != 0);
				break;
			}
			j--;
			size_t ind_1 = std::get<1>(candidate_feature_pairs[j]);
			size_t ind_2 = std::get<2>(candidate_feature_pairs[j]);

			// Check that neither of the verts has been used.
			bool valid_pair = true;
			for (size_t used_ind : indices_used_1) {
				if (ind_1 == used_ind) {
					valid_pair = false;
					break;
				}
			}
			for (size_t used_ind : indices_used_2) {
				if (ind_2 == used_ind) {
					valid_pair = false;
					break;
				}
			}
			if (!valid_pair)
				continue;

			// Check that neither of the verts lies close to a vert that was already used.
			Point p_1 = get_vertex(mesh_1, ind_1)->point();
			Point p_2 = get_vertex(mesh_2, ind_2)->point();
			bool poisson_test_passed = true;
			for (Point p : points_used_1) {
				if (std::sqrt((p - p_1).squared_length()) < POISSON_RADIUS) {
					poisson_test_passed = false;
					break;
				}
			}
			for (Point p : points_used_2) {
				if (std::sqrt((p - p_2).squared_length()) < POISSON_RADIUS) {
					poisson_test_passed = false;
					break;
				}
			}
			if (!poisson_test_passed)
				continue;


			// Success!
			result.push_back(std::make_pair(ind_1, ind_2));
			indices_used_1.insert(ind_1);
			indices_used_2.insert(ind_2);
			points_used_1.push_back(p_1);
			points_used_2.push_back(p_2);

			// Debug.
			double curv_1 = AVERAGING_RADIUS <= 0.0 ? curvature_at_vert(mesh_1, get_vertex(mesh_1, ind_1)) : curvature_at_vert(mesh_1, get_vertex(mesh_1, ind_1), AVERAGING_RADIUS);
			double curv_2 = AVERAGING_RADIUS <= 0.0 ? curvature_at_vert(mesh_2, get_vertex(mesh_2, ind_2)) : curvature_at_vert(mesh_2, get_vertex(mesh_2, ind_2), AVERAGING_RADIUS);
			double dist = std::sqrt((p_1 - p_2).squared_length());
			vout << "Selected a feature pair: v_" << ind_1 << "(" << p_1 << ") : v_" << mesh_1.size_of_vertices() + ind_2 << " (" << p_2 << ") with curvatures " << curv_1 << ", " << curv_2 <<
				" and dist " << dist << ", resulting in correspondence " << std::get<0>(candidate_feature_pairs[j]) << std::endl;

			break;
		}
	}
	return result;
}


///////////////////////////////// 'simple' mode /////////////////////////////////

void simple_features_to_file(const Polyhedron& mesh_1, const Polyhedron& mesh_2, std::string output_path) {
	std::vector<std::pair<size_t, size_t>> features = get_simple_features(mesh_1, mesh_2);
	output_features(features, output_path, mesh_1, mesh_2);
}

std::vector<std::pair<size_t, size_t>> get_simple_features(const Polyhedron& mesh_1, const Polyhedron& mesh_2) {
	VCI vertices_1 = mesh_1.vertices_begin();
	std::vector<std::pair<double, size_t>> features_1;
	for (size_t i = 0; i < mesh_1.size_of_vertices(); i++, vertices_1++) {
		double curv_1 = curvature_at_vert(mesh_1, vertices_1);
		features_1.push_back(std::make_pair(curv_1, i));
	}
	VCI vertices_2 = mesh_2.vertices_begin();
	std::vector<std::pair<double, size_t>> features_2;
	for (size_t i = 0; i < mesh_2.size_of_vertices(); i++, vertices_2++) {
		double curv_2 = curvature_at_vert(mesh_2, vertices_2);
		features_2.push_back(std::make_pair(curv_2, i));
	}
	std::sort(features_1.begin(), features_1.end());
	std::sort(features_2.begin(), features_2.end());

	std::vector<std::pair<size_t, size_t>> features{};
	for (int index : std::vector<int>{ 0, -1, -2 }) { // smallest curvature, largest curvature, 2nd largest curvature
		size_t ind1 = index >= 0 ? index : (size_t)((int)features_1.size() + index);
		size_t ind2 = index >= 0 ? index : (size_t)((int)features_2.size() + index);
		features.push_back(std::make_pair(features_1[ind1].second, features_2[ind2].second));
	}
	return features;
}

///////////////////////////////// 'spatial' mode /////////////////////////////////

void spatial_simple_features_to_file(const Polyhedron& mesh_1, const Polyhedron& mesh_2, std::string output_path) {
	std::vector<std::pair<size_t, size_t>> features = get_simple_features(mesh_1, mesh_2);
	output_features(features, output_path, mesh_1, mesh_2);
}

double compute_triple_distance(size_t i, size_t j, size_t k, const Polyhedron& mesh) {
	double res = 0.0;
	Point pi = get_vertex(mesh, i)->point();
	Point pj = get_vertex(mesh, j)->point();
	Point pk = get_vertex(mesh, k)->point();
	res += std::sqrt((pi - pj).squared_length());
	res += std::sqrt((pj - pk).squared_length());
	res += std::sqrt((pk - pi).squared_length());
	return res;
}

std::vector<size_t> choose_three_to_maximize_distance(const std::vector<size_t>& candidates, const Polyhedron& mesh) {
	size_t best_i = 0;
	size_t best_j = 1;
	size_t best_k = 2;
	double best_dist = CGAL_IA_MAX_DOUBLE;
	for (size_t i = 0; i < candidates.size(); i++) {
		for (size_t j = i + 1; j < candidates.size(); j++) {
			for (size_t k = j + 1; k < candidates.size(); k++) {
				double dist = compute_triple_distance(i, j, k, mesh);
				if (dist < best_dist) {
					best_dist = dist;
					best_i = i;
					best_j = j;
					best_k = k;
				}
			}
		}
	}
	return std::vector<size_t>{best_i, best_j, best_k};
}

double compute_in_pair_distances(std::vector<std::pair<size_t, size_t>> feature_pairs, const Polyhedron& mesh_1, const Polyhedron& mesh_2) {
	double res = 0.0;
	for (std::pair<size_t, size_t> pair : feature_pairs) {
		res += std::sqrt((get_vertex(mesh_1, pair.first)->point() - get_vertex(mesh_2, pair.second)->point()).squared_length());
	}
	return res;
}

// To minimize sum of in-pair distances.
std::vector<std::pair<size_t, size_t>> order_feature_triples(const std::vector<size_t>& features_1, const std::vector<size_t>& features_2, const Polyhedron& mesh_1, const Polyhedron& mesh_2) {
	std::vector<std::pair<size_t, size_t>> res{};
	CGAL_assertion(features_1.size() == 3 && features_2.size() == 3);
	size_t a0 = features_1[0];
	size_t a1 = features_1[1];
	size_t a2 = features_1[2];
	size_t b0 = features_2[0];
	size_t b1 = features_2[1];
	size_t b2 = features_2[2];

	auto feature_pairs_0 = std::vector<std::pair<size_t, size_t>>{ std::make_pair(a0, b0), std::make_pair(a1, b1), std::make_pair(a2, b2) };
	auto feature_pairs_1 = std::vector<std::pair<size_t, size_t>>{ std::make_pair(a0, b0), std::make_pair(a1, b2), std::make_pair(a2, b1) };
	auto feature_pairs_2 = std::vector<std::pair<size_t, size_t>>{ std::make_pair(a0, b1), std::make_pair(a1, b0), std::make_pair(a2, b2) };
	auto feature_pairs_3 = std::vector<std::pair<size_t, size_t>>{ std::make_pair(a0, b1), std::make_pair(a1, b2), std::make_pair(a2, b0) };
	auto feature_pairs_4 = std::vector<std::pair<size_t, size_t>>{ std::make_pair(a0, b2), std::make_pair(a1, b0), std::make_pair(a2, b1) };
	auto feature_pairs_5 = std::vector<std::pair<size_t, size_t>>{ std::make_pair(a0, b2), std::make_pair(a1, b1), std::make_pair(a2, b0) };
	auto feature_pairs_collection = std::vector<std::vector<std::pair<size_t, size_t>>>{ feature_pairs_0, feature_pairs_1, feature_pairs_2, feature_pairs_3, feature_pairs_4, feature_pairs_5 };
	double best_result = CGAL_IA_MAX_DOUBLE;
	for (auto& feature_pairs : feature_pairs_collection) {
		best_result = std::min(best_result, compute_in_pair_distances(feature_pairs, mesh_1, mesh_2));
	}
	for (auto& feature_pairs : feature_pairs_collection) {
		if (best_result == compute_in_pair_distances(feature_pairs, mesh_1, mesh_2))
			return feature_pairs;
	}
	// unreachable
	CGAL_assertion(false);
	return feature_pairs_0;
}

std::vector<std::pair<size_t, size_t>> get_spatial_simple_features(const Polyhedron& mesh_1, const Polyhedron& mesh_2) {
	VCI vertices_1 = mesh_1.vertices_begin();
	std::vector<std::pair<double, size_t>> features_1;
	for (size_t i = 0; i < mesh_1.size_of_vertices(); i++, vertices_1++) {
		double curv_1 = curvature_at_vert(mesh_1, vertices_1);
		features_1.push_back(std::make_pair(curv_1, i));
	}
	VCI vertices_2 = mesh_2.vertices_begin();
	std::vector<std::pair<double, size_t>> features_2;
	for (size_t i = 0; i < mesh_2.size_of_vertices(); i++, vertices_2++) {
		double curv_2 = curvature_at_vert(mesh_2, vertices_2);
		features_2.push_back(std::make_pair(curv_2, i));
	}
	std::sort(features_1.begin(), features_1.end());
	std::sort(features_2.begin(), features_2.end());

	// Could also consider N = min(mesh_1.size_of_vertices(), mesh_2.size_of_vertices()) / c
	const size_t N = 10;

	std::vector<size_t> candidate_features_1{};
	std::vector<size_t> candidate_features_2{};
	for (size_t i = features_1.size() - N; i < features_1.size(); i++) {
		candidate_features_1.push_back(features_1[i].second);
	}
	for (size_t i = features_2.size() - N; i < features_2.size(); i++) {
		candidate_features_2.push_back(features_2[i].second);
	}

	// Choose three high-curv. verts from each mesh such that the triangles have large circumference
	std::vector<size_t> unordered_features_1 = choose_three_to_maximize_distance(candidate_features_1, mesh_1);
	std::vector<size_t> unordered_features_2 = choose_three_to_maximize_distance(candidate_features_2, mesh_2);

	return order_feature_triples(unordered_features_1, unordered_features_2, mesh_1, mesh_2);
}

// Testing
void test_decimation() {
	std::string input_path = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\banana3.obj";
	std::string output_path = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_auto_decimated_banana3.obj";
	std::vector<Point> points;
	std::vector<std::vector<size_t>> faces;
	Polyhedron mesh;
	CGAL::Verbose_ostream vout(true);

	bool succ = input_mesh(points, faces, mesh, vout, input_path);
	if (!succ) {
		std::cerr << "Error reading input meshes." << std::endl;
		CGAL_assertion(false);
	}

	Polyhedron low_res = low_res_mesh(mesh, LOW_RES_NUMBER_OF_EDGES, vout);
	output_mesh(output_path, low_res, vout);
}