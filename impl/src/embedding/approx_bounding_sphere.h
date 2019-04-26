#pragma once

#include "../common/typedefs.h"

Sphere compute_initial_sphere(std::vector<Point>& points);
Sphere compute_approx_bounding_sphere(std::vector<Point>& points);
void normalize_bounding_sphere(std::vector<Point>& points, Sphere sph);