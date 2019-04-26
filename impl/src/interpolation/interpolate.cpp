#include "interpolate.h"

#include "../embedding/io.h"
#include "../overlay/cpp_utils.h"

std::string generate_path(const std::string& base_output_path, size_t ind) {
	// "f" stands for frame, and we number them from 1.
	return base_output_path.substr(0, base_output_path.find_last_of(".")) + "[f" + pad_number(ind + 1, 3) + "]" + base_output_path.substr(base_output_path.find_last_of("."));
}

void interpolate_and_export(const std::vector<Point>& points_1, const std::vector<Point>& points_2, const std::vector<std::vector<size_t>>& faces, double t, std::string output_path_i, CGAL::Verbose_ostream& vout) {
	std::vector<Point> points_interpolated{};
	CGAL_assertion(points_1.size() == points_2.size());
	for (size_t i = 0; i < points_1.size(); i++) {
		points_interpolated.push_back(
			CGAL::ORIGIN + (
				(points_1[i] - CGAL::ORIGIN)*(1.0-t) + 
				(points_2[i] - CGAL::ORIGIN)*(t)
			)
		);
	}
	custom_output_mesh_obj_from_points_and_faces(output_path_i, points_interpolated, faces, vout, true);
}