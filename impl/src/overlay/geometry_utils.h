#pragma once

#include "../common/typedefs.h"

bool identify_vertices(Vector v1, Vector v2);

bool intersect(const Point& p1, const Point& p2, const Point& q1, const Point& q2, Point& output);
bool intersect(const HCI& h1, const HCI& h2, Point& output);

size_t containing_face(const Point& v, const Vector& offset, const Polyhedron& mesh, const FInvIndex& f_index);
size_t containing_face(const Vector& v, const Vector& offset, const Polyhedron& mesh, const FInvIndex& f_index);
size_t containing_face_arbitrary(const Point& v, const Polyhedron& mesh, const FInvIndex& f_index);
size_t containing_face_arbitrary(const Vector& v, const Polyhedron& mesh, const FInvIndex& f_index);

std::vector<double> barycentric_coords_planar(Vector p, Vector a, Vector b, Vector c);
std::vector<double> barycentric_coords_spherical(Vector p, Vector a, Vector b, Vector c);

bool vert_edge_intersection(const Point& vert, const Point& edge_start, const Point& edge_end);
bool vert_edge_intersection(const Point& vert, const HCI& edge);

Vector get_edge_direction_normalized(const HCI& edge);

Vector line_projection(const Vector& a, const Vector&b, const Vector& p);