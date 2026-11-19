#pragma once

#include "lajolla.h"
#include "matrix.h"
#include "point_and_normal.h"
#include "shape.h"
#include "spectrum.h"
#include "texture.h"
#include "vector.h"
#include <variant>
struct Scene;

struct Dynamic
{
    int shape_id;
    Matrix4x4 trans2old; //TODO remove?
    Matrix4x4 inv_trans2old; //TODO remove
    //directly refer to a transformed old shape? for easy access or intersection use later
    Shape old_shape; 
    Shape base_shape;//TODO remove?

    /// for embree use
    RTCGeometry rtc_geometry;
    uint32_t dynamic_geomID;
    uint32_t ghost_geomID;
};



