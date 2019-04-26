#pragma once

// How close two vertices from different meshes must be to be 'merged'.
#define VERTICES_IDENTIFICATION_THRESHOLD 1e-10

// How close a vertex must be to an edge for an intersection to be signalled.
#define VERTEX_EDGE_INTERSECTION_THRESHOLD 1e-10

// How long of a step we take after crossing a vertex to determine which face we've progressed into.
#define STEP_SIZE_ALONG_HEDGE 1e-15

#define DISABLE_USER_INTERACTION true