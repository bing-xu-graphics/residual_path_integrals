#pragma once

#include "vector.h"

/// Convienent class used for storing a point on a surface.
/// Sometimes we will also use it for storing an infinitely far points (environment maps).
/// In that case the normal is the direction of the infinitely far point pointing towards the origin,
/// and position is a point on the scene bounding sphere.
struct PointAndNormal {
	Vector3 position;
	Vector3 normal;
	
	///TODO change to other data structure later. right now add more info for a vertex builder
	//the shape_id is usually known when calling the sample_point, so no need to store that
	int primitive_id = -1;
	Vector2 st; // A 2D parametrization of the surface. Irrelavant to UV mapping.
                // for triangle this is the barycentric coordinates, which we use
                // for interpolating the uv map. 
};
