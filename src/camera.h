#pragma once

#include "lajolla.h"
#include "filter.h"
#include "matrix.h"
#include "vector.h"
#include "ray.h"
#include "spectrum.h"
#include "intersection.h"

/// Currently we only support a pinhole perspective camera
struct Camera
{
    Camera() {}
    Camera(const Matrix4x4 &cam_to_world,
           Real fov, // in degree
           int width, int height,
           Vector3 cam_normal,
           const Filter &filter,
           int medium_id, Real focalDistance = 1.0);

    Matrix4x4 sample_to_cam, cam_to_sample;
    Matrix4x4 cam_to_world, world_to_cam;
    int width, height;
    Real focalDistance;
    Vector3 cam_normal;
    Real image_plane_distance;
    Filter filter;

    int medium_id; // for participating media rendering in homework 2
};

/// Given screen position in [0, 1] x [0, 1],
/// generate a camera ray.
Ray sample_primary(const Camera &camera,
                   const Vector2 &screen_pos,
                   const int time = 0);

Spectrum sample_Wi(const Camera &camera, const PathVertex &vertex,
                   Vector3 &sensorPos, Vector3 &dirToSensor,
                   Real &cosAtSensor, Vector2 &pRaster2,
                   Real &pdfSensorW,
                   const int time = 0);
Spectrum sample_Wi(const Camera &camera, const PathVertex &vertex,
                   Vector3 &sensorPos, Vector3 &dirToSensor,
                   Real &cosAtSensor, Vector2 &pRaster2,
                   Real &pdfSensorW,Real &pdfAreaVertex,
                   const int time = 0);
