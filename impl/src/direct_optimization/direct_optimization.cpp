#include "direct_optimization.h"
#include "../feature_alignment/grid_minimization.h"
#include "../embedding/geometry_utils.h"

Polyhedron rotated_polyhedron(Polyhedron p, double x, double y, double z) {
	Transformation t = rotZ(z)*rotY(y)*rotX(x);
	Polyhedron res{ p };
	for (auto vertex_it = res.vertices_begin(); vertex_it != res.vertices_end(); vertex_it++) {
		(*vertex_it).point() = t.transform(vertex_it->point());
	}
	return res;
}

Polyhedron optimal_rotation(
	Polyhedron& mesh_1, // the embedding to be modified
	Polyhedron& mesh_2,
	Polyhedron& shape_1,
	Polyhedron& shape_2,
	std::function<double(Polyhedron&, Polyhedron&, Polyhedron&, Polyhedron&)> f_to_minimize
) {
	std::function<double(double, double, double)> f_grid = [mesh_1, &mesh_2, &shape_1, &shape_2, f_to_minimize](double x, double y, double z) {
		Polyhedron mesh_1r = rotated_polyhedron(mesh_1, x, y, z);
		return f_to_minimize(mesh_1r, mesh_2, shape_1, shape_2);
	};

	const size_t grid_size = 2;
	double grid_centre_x = 0.0;
	double grid_centre_y = 0.0;
	double grid_centre_z = 0.0;
	double grid_width_x = 360.0;
	double grid_width_y = 360.0;
	double grid_width_z = 360.0;

#define STEPS 0
	// Hardcoded for now.
	double grid_width_multipliers[((STEPS) == 0) ? 1 : (STEPS)]{ 1.0/((double)grid_size), };

	double best_x = 0.0;
	double best_y = 0.0;
	double best_z = 0.0;
	double best_val = CGAL_IA_MAX_DOUBLE;

	std::vector<std::pair<Grid3D, size_t/*iteration_number*/>> grids_to_do{};

	grids_to_do.push_back(std::make_pair(Grid3D{ grid_size, grid_centre_x - grid_width_x / 2.0, grid_centre_x + grid_width_x / 2.0, grid_centre_y - grid_width_y / 2.0, grid_centre_y + grid_width_y / 2.0, grid_centre_z - grid_width_z / 2.0, grid_centre_z + grid_width_z / 2.0, f_grid }, 0));

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
					Grid3D{ grid_size, new_local_min.x - this_grid_width_x / 2.0, new_local_min.x + this_grid_width_x / 2.0, new_local_min.y - this_grid_width_y / 2.0, new_local_min.y + this_grid_width_y / 2.0, new_local_min.z - this_grid_width_z / 2.0, new_local_min.z + this_grid_width_z / 2.0, f_grid },
					iteration_number + 1
				));
			}
		}
	}
#undef STEPS
	// Debug.
	std::cout << "Optimal rotation is rotZ(" << best_z << ")*rotY(" << best_y << ")*rotX(" << best_x << ")." << std::endl;

	return rotated_polyhedron(mesh_1, best_x, best_y, best_z);
}

std::pair<Polyhedron, Polyhedron> optimize_embedding(
	Polyhedron& mesh_1, // the embedding to be modified
	Polyhedron& mesh_2,
	Polyhedron& shape_1,
	Polyhedron& shape_2,
	std::function<double(Polyhedron&, Polyhedron&, Polyhedron&, Polyhedron&)> f_to_minimize
) {
	Polyhedron mesh_1_rotated = optimal_rotation(mesh_1, mesh_2, shape_1, shape_2, f_to_minimize);

	// TODO More interesting optimization.

	return std::make_pair(mesh_1_rotated, mesh_2);
}