#pragma once
#include"../common/typedefs.h"

struct GridPoint {
public:
	double x;
	double y;
	double val;
};

class Grid {
public:
	size_t size;
	double x_min;
	double x_max;
	double y_min;
	double y_max;

	std::function<double(double, double)> f;

	Grid(size_t _size, double _x_min, double _x_max, double _y_min, double _y_max, std::function<double(double, double)> _f) :
		size(_size), x_min(_x_min), x_max(_x_max), y_min(_y_min), y_max(_y_max), f(_f) {}

	std::vector<GridPoint> local_minima();
};


struct GridPoint3D {
public:
	double x;
	double y;
	double z;
	double val;
};

class Grid3D {
public:
	size_t size;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;

	std::function<double(double, double, double)> f;

	Grid3D(size_t _size, double _x_min, double _x_max, double _y_min, double _y_max, double _z_min, double _z_max, std::function<double(double, double, double)> _f) :
		size(_size), x_min(_x_min), x_max(_x_max), y_min(_y_min), y_max(_y_max), z_min(_z_min), z_max(_z_max), f(_f) {}

	std::vector<GridPoint3D> local_minima();
};
