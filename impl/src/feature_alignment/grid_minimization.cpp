#include "grid_minimization.h"
#include<unordered_map>

std::vector<GridPoint> Grid::local_minima(){
	std::unordered_map<size_t, std::unordered_map<size_t, double>> values{};
	// collect values
	for (size_t xi = 0; xi < size; xi++) {
		double x = x_min + (x_max - x_min)*((double(xi)) / ((double)size));
		for (size_t yi = 0; yi < size; yi++) {
			double y = y_min + (y_max - y_min)*((double(yi)) / ((double)size));
			double val = f(x, y);
			values[xi][yi] = val;
		}
	}
	std::vector<GridPoint> minima{};
	// find minima
	for (size_t xi = 0; xi < size; xi++) {
		double x = x_min + (x_max - x_min)*((double(xi)) / ((double)size));
		for (size_t yi = 0; yi < size; yi++) {
			double y = y_min + (y_max - y_min)*((double(yi)) / ((double)size));

			double val = values.at(xi).at(yi);
			bool is_a_minimum = true;
			for (int dx = -1; dx <= 1; dx++) {
				if (!is_a_minimum) break;
				for (int dy = -1; dy <= 1; dy++) {
					if (dx == 0 && dy == 0) continue;
					if (xi + dx < 0 || xi + dx >= size || yi + dy < 0 || yi + dy >= size) continue;
					if (values.at(xi + dx).at(yi + dy) < val) {
						is_a_minimum = false;
						break;
					}
				}
			}
			if (is_a_minimum) {
				minima.push_back(GridPoint{ x, y, val });
			}
		}
	}
	return minima;
}


std::vector<GridPoint3D> Grid3D::local_minima() {
	std::unordered_map<size_t, std::unordered_map<size_t, std::unordered_map<size_t, double>>> values{};
	// collect values
	for (size_t xi = 0; xi < size; xi++) {
		double x = x_min + (x_max - x_min)*((double(xi)) / ((double)size));
		for (size_t yi = 0; yi < size; yi++) {
			double y = y_min + (y_max - y_min)*((double(yi)) / ((double)size));
			for (size_t zi = 0; zi < size; zi++) {
				double z = z_min + (z_max - z_min)*((double(zi)) / ((double)size));
				double val = f(x, y, z);
				values[xi][yi][zi] = val;
			}
		}
	}
	std::vector<GridPoint3D> minima{};
	// find minima
	for (size_t xi = 0; xi < size; xi++) {
		double x = x_min + (x_max - x_min)*((double(xi)) / ((double)size));
		for (size_t yi = 0; yi < size; yi++) {
			double y = y_min + (y_max - y_min)*((double(yi)) / ((double)size));
			for (size_t zi = 0; zi < size; zi++) {
				double z = z_min + (z_max - z_min)*((double(zi)) / ((double)size));
				double val = values.at(xi).at(yi).at(zi);
				bool is_a_minimum = true;
				for (int dx = -1; dx <= 1; dx++) {
					if (!is_a_minimum) break;
					for (int dy = -1; dy <= 1; dy++) {
						for (int dz = -1; dz <= 1; dz++) {
							if (dx == 0 && dy == 0 && dz == 0) continue;
							if (xi + dx < 0 || xi + dx >= size || yi + dy < 0 || yi + dy >= size || zi + dz < 0 || zi + dz >= size) continue;
							if (values.at(xi + dx).at(yi + dy).at(zi + dz) < val) {
								is_a_minimum = false;
								break;
							}
						}
					}
				}
				if (is_a_minimum) {
					minima.push_back(GridPoint3D{ x, y, z, val });
				}
			}
			
		}
	}
	return minima;
}