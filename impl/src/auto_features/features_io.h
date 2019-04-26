#pragma once
#include "../common/typedefs.h"

bool output_features(const std::vector<std::pair<size_t, size_t>>& features, std::string path);

bool output_features(const std::vector<std::pair<size_t, size_t>>& features, std::string path, const Polyhedron& mesh_1, const Polyhedron& mesh_2);