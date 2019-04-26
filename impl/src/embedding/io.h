#pragma once
#include "../common/typedefs.h"

std::string get_path_with_underscore(const std::string& path);
std::string get_path_without_underscore(const std::string& path);
std::string get_path_with_underscore_toggled(const std::string& path);

bool input_mesh(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh, CGAL::Verbose_ostream& vout, std::string& input_filename);
bool output_mesh(std::string path, Polyhedron& mesh, CGAL::Verbose_ostream& vout, bool printing_on = true);
void parse_args(const int argc, char** argv, std::string& input_path, std::string& output_path, bool& verbose, double& RELAXATION_EPSILON, Sphere& sph);
void make_mesh_from_points_and_faces(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh);
void make_mesh_from_points_and_faces_skipping_orientation(std::vector<Point>& points, std::vector<std::vector<size_t>>& faces, Polyhedron& mesh);

void ask_before_termination();
void ask_before_termination(bool disable_user_interaction);

void get_faces(std::vector<std::vector<size_t>>& faces, const Polyhedron& mesh);
void get_points(std::vector<Point>& points, const Polyhedron& mesh);

bool custom_output_mesh_obj_from_points_and_faces(std::string path, const std::vector<Point>& points, const std::vector<std::vector<size_t>>& faces, CGAL::Verbose_ostream& vout, bool printing_on);




void test_path_underscore_functions();