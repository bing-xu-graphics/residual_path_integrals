#include "camera.h"
#include "lajolla.h"
#include "transform.h"

#include <cmath>

Camera::Camera(const Matrix4x4 &cam_to_world,
               Real fov,
               int width, int height,
               Vector3 cam_normal,
               const Filter &filter,
               int medium_id,
               Real focalDistance
               )
    : cam_to_world(cam_to_world),
      world_to_cam(inverse(cam_to_world)),
      width(width), height(height),
      filter(filter), medium_id(medium_id), focalDistance(focalDistance), cam_normal(cam_normal)
{
  Real aspect = (Real)width / (Real)height;
  cam_to_sample = scale(Vector3(-Real(0.5), -Real(0.5) * aspect, Real(1.0))) *
                  translate(Vector3(-Real(1.0), -Real(1.0) / aspect, Real(0.0))) *
                  perspective(fov);
  sample_to_cam = inverse(cam_to_sample);
  image_plane_distance = (Real)width / (2.0 * tan(fov * c_PI / 360.0));
}

Ray sample_primary(const Camera &camera,
                   const Vector2 &screen_pos,
                   const int time)
{
  // screen_pos' domain is [0, 1]^2
  Vector2 pixel_pos{screen_pos.x * camera.width, screen_pos.y * camera.height};

  // Importance sample from the pixel filter (see filter.h for more explanation).
  // We first extract the subpixel offset.
  Real dx = pixel_pos.x - floor(pixel_pos.x);
  Real dy = pixel_pos.y - floor(pixel_pos.y);
  // dx/dy are uniform variables in [0, 1],
  // so let's remap them using importance sampling.
  Vector2 offset = sample(camera.filter, Vector2{dx, dy});
  // Filters are placed at pixel centers.
  Vector2 remapped_pos{
      (floor(pixel_pos.x) + Real(0.5) + offset.x) / camera.width,
      (floor(pixel_pos.y) + Real(0.5) + offset.y) / camera.height};

  Vector3 pt = xform_point(camera.sample_to_cam,
                           Vector3(remapped_pos[0], remapped_pos[1], Real(0)));
  Vector3 dir = normalize(pt);
  return Ray{xform_point(camera.cam_to_world, Vector3{0, 0, 0}),
             // the last normalize might not be necessary
             normalize(xform_vector(camera.cam_to_world, dir)),
             Real(0), infinity<Real>(), time};
}

Spectrum sample_Wi(const Camera &camera, const PathVertex &vertex,
                   Vector3 &sensorPos, Vector3 &dirToSensor,
                   Real &cosAtSensor, Vector2 &pRaster2,
                   Real &pdfSensorW,
                   int time)
{
  Vector3 pLensWorld = xform_point(camera.cam_to_world, Vector3{0.0, 0.0, 0.0});
  sensorPos = pLensWorld;
  Vector3 camNormal = camera.cam_normal;
  dirToSensor = normalize(pLensWorld - vertex.position);
  if (dot(camNormal, -dirToSensor) <= 0.0)
    return make_zero_spectrum(); // check point is in front of the camera
  // check it projects onto the screen (and where)
  Vector3 refPoint = xform_point(camera.world_to_cam, vertex.position); // in local coordinate
  Vector3 screenPos = xform_point(camera.cam_to_sample, refPoint);
  pRaster2.x = screenPos.x;
  pRaster2.y = screenPos.y;
  if (pRaster2.x < 0 || pRaster2.x > 1.0 || pRaster2.y < 0 || pRaster2.y > 1.0)
  {
    pdfSensorW = 0.0;
    return make_zero_spectrum();
  }
  pRaster2.x *= camera.width;
  pRaster2.y *= camera.height;

  // compute distance and normalize direction to the camera
  Real dist = length(dirToSensor);
  dirToSensor /= dist;
  cosAtSensor = dot(camNormal, -dirToSensor);

  // compute pdf conversion factor from image plane area to surface area
  Real imagePointToSensorDist = camera.image_plane_distance / cosAtSensor;
  Real imageToSolidAngleFactor = imagePointToSensorDist * imagePointToSensorDist / cosAtSensor; // pdfSensorW
  Real cosAtVertex = dot(dirToSensor, vertex.geometric_normal);
  Real imageToSurfaceFactor = imageToSolidAngleFactor * fabs(cosAtVertex) / dist / dist; // pdfVertexA
  Real cameraPdfA = imageToSurfaceFactor;
  pdfSensorW = imageToSolidAngleFactor;
  // We divide the contribution by surfaceToImageFactor to convert the (already
  // divided) pdf from surface area to image plane area, w.r.t. which the
  // pixel integral is actually defined.
  Real G = cosAtSensor / dist / dist;
  return make_const_spectrum(imageToSurfaceFactor);
}
Spectrum sample_Wi(const Camera &camera, const PathVertex &vertex,
                   Vector3 &sensorPos, Vector3 &dirToSensor,
                   Real &cosAtSensor, Vector2 &pRaster2,
                   Real &pdfSensorW,Real &pdfAreaVertex,
                   const int time)
{
  Vector3 pLensWorld = xform_point(camera.cam_to_world, Vector3{0.0, 0.0, 0.0});
  Vector3 camNormal = camera.cam_normal;
  sensorPos = pLensWorld;
  dirToSensor = pLensWorld - vertex.position;
  if (dot(camNormal, -dirToSensor) <= 0.0)
    return make_zero_spectrum(); // check point is in front of the camera
  // check it projects onto the screen (and where)
  Vector3 refPoint = xform_point(camera.world_to_cam, vertex.position); // in local coordinate
  Vector3 screenPos = xform_point(camera.cam_to_sample, refPoint);
  pRaster2.x = screenPos.x;
  pRaster2.y = screenPos.y;
  if (pRaster2.x < 0 || pRaster2.x > 1.0 || pRaster2.y < 0 || pRaster2.y > 1.0)
  {
    pdfSensorW = 0.0;
    return make_zero_spectrum();
  }
  pRaster2.x *= camera.width;
  pRaster2.y *= camera.height;

  // compute distance and normalize direction to the camera
  Real dist = length(dirToSensor);
  dirToSensor /= dist;
  cosAtSensor = dot(camNormal, -dirToSensor);

  // compute pdf conversion factor from image plane area to surface area
  Real imagePointToSensorDist = camera.image_plane_distance / cosAtSensor;
  Real imageToSolidAngleFactor = imagePointToSensorDist * imagePointToSensorDist / cosAtSensor; // pdfSensorW
  Real cosAtVertex = dot(dirToSensor, vertex.geometric_normal);
  //TODO what if cosAtVertex <=0
  Real imageToSurfaceFactor = imageToSolidAngleFactor * fabs(cosAtVertex) / dist / dist; // pdfVertexA
  Real cameraPdfA = imageToSurfaceFactor;
  pdfSensorW = imageToSolidAngleFactor;
  pdfAreaVertex = cameraPdfA;
  // We divide the contribution by surfaceToImageFactor to convert the (already
  // divided) pdf from surface area to image plane area, w.r.t. which the
  // pixel integral is actually defined.
  Real G = cosAtSensor / dist / dist;
  return make_const_spectrum(imageToSurfaceFactor);
}
/**
 * refactored from mitsuba
 * \brief Compute the directional sensor response function
 * of the camera multiplied with the cosine foreshortening
 * factor associated with the image plane
 *
 * \param d
 *     A normalized direction vector from the aperture position to the
 *     reference point in question (all in local camera space)
 */

/**
 * @brief We() method takes a ray with origin p and direction w
 *  and evaluates the importance emitted from the point on the camera  in a direction w.
 *  pRaster2 stores the raster position associated with the ray on the film;
 */
