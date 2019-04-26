#pragma once

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Random.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double>                      Kernel;
typedef Kernel::Point_3                                     Point;
typedef Kernel::Vector_3                                    Vector;
typedef CGAL::Polyhedron_3<Kernel>                          Polyhedron;
typedef Polyhedron::Vertex_iterator                         Vertex_iterator;
typedef CGAL::Random                                        Random;
typedef CGAL::Aff_transformation_3<Kernel>                  Transformation;
typedef Polyhedron::Vertex_const_iterator                   VCI;
typedef Polyhedron::Halfedge_const_iterator					HCI;
typedef Polyhedron::Facet_const_iterator                    FCI;
typedef Polyhedron::Vertex_const_handle                     Vertex;
typedef Polyhedron::Halfedge_const_handle                   Halfedge;
typedef Polyhedron::Facet_const_handle                      Facet;
typedef Polyhedron::Halfedge_around_facet_const_circulator  HFCC;
typedef Polyhedron::Halfedge_around_vertex_const_circulator HVCC;
typedef CGAL::Inverse_index<VCI>                            VInvIndex;
typedef CGAL::Inverse_index<HCI>                            HInvIndex;
typedef CGAL::Inverse_index<FCI>                            FInvIndex;

typedef struct _Sphere {
	double radius;
	Point centre;
} Sphere;