#pragma once
#include "../common/typedefs.h"

void simple_features_to_file(const Polyhedron& mesh_1, const Polyhedron& mesh_2, std::string output_path);
std::vector<std::pair<size_t, size_t>> get_simple_features(const Polyhedron& mesh_1, const Polyhedron& mesh_2);

void spatial_simple_features_to_file(const Polyhedron& mesh_1, const Polyhedron& mesh_2, std::string output_path);
std::vector<std::pair<size_t, size_t>> get_spatial_simple_features(const Polyhedron& mesh_1, const Polyhedron& mesh_2);

void poisson_features_to_file(const Polyhedron& mesh_1, const Polyhedron& mesh_2, const Polyhedron& low_res_1, const Polyhedron& low_res_2, 
	double PERCENTAGE_OF_VERTS_TO_CHECK, size_t N_FEATURES, double POISSON_RADIUS, double DIST_WEIGHT, double AVERAGING_RADIUS, std::string output_path, 
	CGAL::Verbose_ostream& vout);

std::vector<std::pair<size_t, size_t>> get_poisson_features(const Polyhedron& mesh_1, const Polyhedron& mesh_2, const Polyhedron& low_res_1, const Polyhedron& low_res_2,
	double PERCENTAGE_OF_VERTS_TO_CHECK /* 0.1 */, size_t N_FEATURES /* 3 */, double POISSON_RADIUS /* 0.3 */, double DIST_WEIGHT /* 1.0 */, double AVERAGING_RADIUS /* 0.3 */,
	CGAL::Verbose_ostream& vout);


// Testing.
void test_decimation();
