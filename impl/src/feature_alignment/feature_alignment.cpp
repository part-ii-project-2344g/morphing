#include "feature_alignment.h"
#include "../overlay/cpp_utils.h"
#include "../embedding/geometry_utils.h"
#include "../embedding/io.h"

#include "grid_minimization.h"

#define MAX_DISPLACEMENT_THRESHOLD 5e-16

double feature_distance(const Polyhedron& mesh_1, const Polyhedron& mesh_2, std::vector<std::pair<size_t, size_t>> features, double rx, double ry, double rz) {
	double res = 0.0;
	Transformation t = rotZ(rz)*rotY(ry)*rotX(rx);
	for (auto pair : features) {
		size_t vert1 = pair.first;
		size_t vert2 = pair.second;
		Vector v1 = get_vertex(mesh_1, vert1)->point() - CGAL::ORIGIN;
		Vector v2 = get_vertex(mesh_2, vert2)->point() - CGAL::ORIGIN;
		Vector v1r = t.transform(v1);
		double distSq = (v1r - v2).squared_length();
		res += distSq;
	}
	return res;
}

Polyhedron rotate_mesh_1(
	const Polyhedron& mesh_1,
	const Polyhedron& mesh_2,
	const std::vector<std::pair<size_t, size_t>>& features
) {
	// f = the function to minimize
	std::function<double(double, double, double)> f = [mesh_1, mesh_2, features](double x, double y, double z) { return feature_distance(mesh_1, mesh_2, features, x, y, z); };
	const size_t grid_size = 4;
	double grid_centre_x = 0.0;
	double grid_centre_y = 0.0;
	double grid_centre_z = 0.0;
	double grid_width_x = 360.0;
	double grid_width_y = 360.0;
	double grid_width_z = 360.0;

#define STEPS 2
	// Hardcoded for now.
	double grid_width_multipliers[STEPS]{ 0.25, 0.25/4.0 };

	double best_x = 0.0;
	double best_y = 0.0;
	double best_z = 0.0;
	double best_val = CGAL_IA_MAX_DOUBLE;

	std::vector<std::pair<Grid3D, size_t/*iteration_number*/>> grids_to_do{};

	grids_to_do.push_back(std::make_pair(Grid3D{ grid_size, grid_centre_x - grid_width_x / 2.0, grid_centre_x + grid_width_x / 2.0, grid_centre_y - grid_width_y / 2.0, grid_centre_y + grid_width_y / 2.0, grid_centre_z - grid_width_z / 2.0, grid_centre_z + grid_width_z / 2.0, f }, 0));

	// Do the grid minimization a few (STEPS) times in a row with progressively finer detail and find the global 'minimum', i.e. approximately the best rotation of mesh_1.	
	while (!grids_to_do.empty()) {
		auto pair = grids_to_do[grids_to_do.size() - 1];
		grids_to_do.pop_back();
		size_t iteration_number = pair.second;
		Grid3D grid = pair.first;

		std::vector<GridPoint3D> gps = grid.local_minima();
		for (GridPoint3D new_local_min : gps) {
			if (new_local_min.val < best_val) {
				best_val = new_local_min.val;
				best_x = new_local_min.x;
				best_y = new_local_min.y;
				best_z = new_local_min.z;
			}
			// If we haven't crossed the step limit, push new grids to the stack to investigate the newly found local minima closer.
			if (iteration_number < STEPS) {
				double this_grid_width_x = grid_width_x * grid_width_multipliers[iteration_number];
				double this_grid_width_y = grid_width_y * grid_width_multipliers[iteration_number];
				double this_grid_width_z = grid_width_z * grid_width_multipliers[iteration_number];
				grids_to_do.push_back(std::make_pair(
					Grid3D{ grid_size, new_local_min.x - this_grid_width_x / 2.0, new_local_min.x + this_grid_width_x / 2.0, new_local_min.y - this_grid_width_y / 2.0, new_local_min.y + this_grid_width_y / 2.0, new_local_min.z - this_grid_width_z / 2.0, new_local_min.z + this_grid_width_z / 2.0, f },
					iteration_number + 1
				));
			}
		}
	}
#undef STEPS
	// Debug.
	std::cout << "Optimal rotation is rotY(" << best_y << ")*rotX(" << best_x << ")." << std::endl;

	Transformation t = rotZ(best_z)*rotY(best_y)*rotX(best_x);
	Polyhedron mesh_1r{ mesh_1 };
	for (auto vertex_it = mesh_1r.vertices_begin(); vertex_it != mesh_1r.vertices_end(); vertex_it++) {
		(*vertex_it).point() = t.transform(vertex_it->point());
	}
	
	return mesh_1r;
}

// Apply the map to a single point of mesh_1.
Point map(Vector p, double c, double d, const Point& arg, Point feature_1) {
	double dist = sqrt((arg - feature_1).squared_length());
	if (dist < d) {
		return normalized(arg + c * p * (d - dist));
	}
	else { // dist >= d
		return arg;
	}
}

// Apply the map to all vertices of mesh_1 and normalize them.
Polyhedron apply_map_for_feature_i(Polyhedron mesh_1/*pass by copy*/, Polyhedron& mesh_2, 
	size_t v_ind_1, // index of feature vertex at mesh_1 
	size_t v_ind_2, // index of corresponding feature vertex at mesh_2
	double c, double d) {
	Point feature_1 = get_vertex(mesh_1, v_ind_1)->point();
	Point feature_2 = get_vertex(mesh_2, v_ind_2)->point();
	Vector p = feature_2 - feature_1;

	Vertex_iterator it1 = mesh_1.vertices_begin();
	for (; it1 != mesh_1.vertices_end(); it1++) {
		it1->point() = map(p, c, d, it1->point(), feature_1);
	}

	return mesh_1;
}

bool is_orientation_consistent(Polyhedron& mesh) {
	std::vector<Point> points;
	std::vector<std::vector<size_t>> faces;
	get_points(points, mesh);
	get_faces(faces, mesh);
	return is_orientation_consistent_on_sphere(points, faces);
}

// Perform step 2 (step 3 if mesh_1 and mesh_2 are passed swapped).
std::pair<Polyhedron, Polyhedron> step23(Polyhedron mesh_1/* pass by copy */, Polyhedron mesh_2/* pass by copy */, std::vector<std::pair<size_t, size_t>> features, bool meshes_swapped, double d, bool& out_success, const std::string& debug_output_path, CGAL::Verbose_ostream& vout) {
	for (auto p : features) {
		std::cout << "Trying next feature in step23()." << std::endl;
		size_t v_ind_1 = meshes_swapped ? p.second : p.first;
		size_t v_ind_2 = meshes_swapped ? p.first : p.second;
		double c = 5.0; // The first c attempted will be 0.5.
		
		bool orientation_consistent = false;
		while (!orientation_consistent) {
			c *= 0.1; // TODO: Tweak?
			std::cout << "Trying c=" << c << std::endl;

			if (c <= 0.0) {
				std::cout << "step23's c dropped to 0.0, failing." << std::endl;
				out_success = false;
				return std::make_pair(Polyhedron{}, Polyhedron{});
			}

			Polyhedron mesh_1_new = apply_map_for_feature_i(mesh_1, mesh_2, v_ind_1, v_ind_2, c, d);
			orientation_consistent = is_orientation_consistent(mesh_1_new);
			
			// Debug output.
			std::string path = debug_output_path;
			std::string debug_output_path_specified = path.substr(0, path.find_last_of("\\/") + 1) + "[m" + (meshes_swapped ? "2" : "1") + "c" + std::to_string(c) + ",d" + std::to_string(d) + "]" + path.substr(path.find_last_of("\\/") + 1);
			// output_mesh(debug_output_path_specified, mesh_1_new, vout);

			if (orientation_consistent) {
				// Update state and move on to the next feature. 
				// Otherwise we discard mesh_1_new and retry with smaller c.
				mesh_1 = mesh_1_new;
			}
		}
	}
	out_success = true;
	return std::make_pair(mesh_1, mesh_2);
}

std::pair<Polyhedron, Polyhedron> step2(Polyhedron& mesh_1, Polyhedron& mesh_2, std::vector<std::pair<size_t, size_t>> features, double d, bool& out_success, const std::string& debug_output_path, CGAL::Verbose_ostream& vout) {
	return step23(mesh_1, mesh_2, features, false, d, out_success, debug_output_path, vout);
}

std::pair<Polyhedron, Polyhedron> step3(Polyhedron& mesh_1, Polyhedron& mesh_2, std::vector<std::pair<size_t, size_t>> features, double d, bool& out_success, const std::string& debug_output_path, CGAL::Verbose_ostream& vout) {
	return step23(mesh_2, mesh_1, features, true, d, out_success, debug_output_path, vout);
}

double get_max_displacement(Polyhedron& mesh_1, Polyhedron& mesh_1_new) {
	CGAL_assertion(mesh_1.size_of_vertices() == mesh_1_new.size_of_vertices());

	double res_sq = 0.0;

	Vertex_iterator it1 = mesh_1.vertices_begin();
	Vertex_iterator it2 = mesh_1_new.vertices_begin();
	for (size_t i = 0; i < mesh_1.size_of_vertices(); i++, it1++, it2++) {
		double dist_sq = (it1->point() - it2->point()).squared_length();
		res_sq = std::max(res_sq, dist_sq);
		//std::cout << "Old: " << it1->point() << ", new: " << it2->point() << std::endl;
	}
	return std::sqrt(res_sq);
}

double get_max_displacement(Polyhedron& mesh_1, Polyhedron& mesh_2, Polyhedron& mesh_1_new, Polyhedron& mesh_2_new) {
	return std::max(get_max_displacement(mesh_1, mesh_1_new), get_max_displacement(mesh_2, mesh_2_new));
}

bool are_points_identical(const Point& p1, const Point& p2) {
	return sqrt((p1 - p2).squared_length()) < 1e-14; // TODO: Tweak?
}

bool do_features_coincide(Polyhedron& mesh_1, Polyhedron& mesh_2, std::vector<std::pair<size_t, size_t>> features) {
	for (auto p : features) {
		size_t v1 = p.first;
		size_t v2 = p.second;
		Point p1 = get_vertex(mesh_1, v1)->point();
		Point p2 = get_vertex(mesh_2, v2)->point();
		if (!are_points_identical(p1, p2)) return false;
	}
	return true;
}

std::pair<Polyhedron, Polyhedron> align_features(
	Polyhedron& mesh_1_raw, 
	Polyhedron mesh_2, // pass by copy
	std::vector<std::pair<size_t, size_t>> features,
	const std::string& debug_out_path,
	CGAL::Verbose_ostream& vout
) {
	Polyhedron mesh_1 = rotate_mesh_1(mesh_1_raw, mesh_2, features);

	double d = 4.0; // 2.0 will be the first used value and results in a fully global transformation.
	bool features_coincide = false;

	while (!features_coincide) {
		if (d == 0.0) {
			std::cerr << "d dropped to 0.0 :(" << std::endl;
			CGAL_assertion(false);
			return std::make_pair(Polyhedron{}, Polyhedron{});
		}
		d *= 0.5; // TODO: Tweak?
		std::cout << "Trying d=" << d << std::endl;

		Polyhedron mesh_1_copy{ mesh_1 };
		Polyhedron mesh_2_copy{ mesh_2 };

		double max_displacement = CGAL_IA_MAX_DOUBLE;
		bool success = false;
		while (max_displacement > MAX_DISPLACEMENT_THRESHOLD) {
			// Step 2.
			auto meshes_new = step2(mesh_1_copy, mesh_2_copy, features, d, success, debug_out_path, vout);
			if (!success) break;
			Polyhedron mesh_1_new = meshes_new.first;
			Polyhedron mesh_2_new = meshes_new.second;
			
			// Step 3.
			// If c drops to 0.0 inside of step3 (or step2), set success to false, and retry this outer while loop with a smaller d.
			meshes_new = step3(mesh_1_new, mesh_2_new, features, d, success, debug_out_path, vout); 
			if (!success) break;
			// The output of step3 is reversed!
			mesh_1_new = meshes_new.second;
			mesh_2_new = meshes_new.first;

			// Compute max_displacement of step 2 and 3 above.
			max_displacement = get_max_displacement(mesh_1_copy, mesh_2_copy, mesh_1_new, mesh_2_new);
			// TODO: try to get displacement only of feature verts!
			std::cout << "max_displacement = " << max_displacement << std::endl;
			//std::cout << "mesh_1_copy.vertices_begin()->point() = " << mesh_1_copy.vertices_begin()->point() << std::endl;
			//std::cout << "++mesh_1_copy.vertices_begin()->point() = " << (++mesh_1_copy.vertices_begin())->point() << std::endl;
			//std::cout << "++++mesh_1_copy.vertices_begin()->point() = " << (++(++mesh_1_copy.vertices_begin()))->point() << std::endl;
			//std::cout << "++++++mesh_1_copy.vertices_begin()->point() = " << (++(++(++mesh_1_copy.vertices_begin())))->point() << std::endl;

			// Update state.
			mesh_1_copy = mesh_1_new;
			mesh_2_copy = mesh_2_new;
		}

		if (!success) continue;

		features_coincide = do_features_coincide(mesh_1_copy, mesh_2_copy, features);
		if (features_coincide) {
			// We're done.
			return std::make_pair(mesh_1_copy, mesh_2_copy);
		}
	}

	CGAL_assertion(false); // Unreachable code.
	return std::make_pair(Polyhedron{}, Polyhedron{});
}