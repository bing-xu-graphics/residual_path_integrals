#include "intersection.h"
#include "material.h"
#include "ray.h"
#include "scene.h"

std::optional<PathVertex> intersect(const Scene &scene,
                                    Ray &ray,
                                    const RayDifferential &ray_diff)
{
    IntersectContext user_context; // there are rtc_context inside
    InitIntersectionContext(&user_context);
    RTCRayHit rtc_rayhit;
    RTCRay &rtc_ray = rtc_rayhit.ray;
    RTCHit &rtc_hit = rtc_rayhit.hit;
    rtc_ray = RTCRay{
        (float)ray.org.x, (float)ray.org.y, (float)ray.org.z,
        (float)ray.tnear,
        (float)ray.dir.x, (float)ray.dir.y, (float)ray.dir.z,
        (float)ray.time, // time
        (float)ray.tfar,
        (unsigned int)(-1), // mask
        0,                  // ray ID
        0                   // ray flags
    };
    rtc_hit = RTCHit{
        0, 0, 0, // Ng_x, Ng_y, Ng_z
        0,
        0,                        // u, v
        RTC_INVALID_GEOMETRY_ID,  // primitive ID
        RTC_INVALID_GEOMETRY_ID,  // geometry ID
        {RTC_INVALID_GEOMETRY_ID} // instance IDs
    };
    rtcIntersect1(scene.embree_scene, &user_context.rtc_context, &rtc_rayhit);
    if (rtc_hit.geomID == RTC_INVALID_GEOMETRY_ID)
    {
        { // don't contaminate the ghost record; will test ghost separately as needed
            const Payload *payload = (const Payload *)user_context.payload;
            assert(payload->dist_t.size() == payload->ghost_shape_ids.size());
            int n = payload->dist_t.size();
            for (int i = 0; i < n; ++i)
            {
                Real dist = payload->dist_t[i];
                if (dist <= rtc_ray.tfar && dist > rtc_ray.tnear)
                {
                    ray.num_ghost_hits += 1;
                    ray.ghost_shape_ids.push_back(payload->ghost_shape_ids[i]);
                    ray.ghost_hit_positions.push_back(payload->ghost_hit_positions[i]);
                    ray.ghost_hit_geom_normals.push_back(payload->ghost_hit_geom_normals[i]);
                }
            }
        }
        return {};
    };
    assert(rtc_hit.geomID < scene.shapes.size());

    PathVertex vertex;
    vertex.position = Vector3{ray.org.x, ray.org.y, ray.org.z} +
                      Vector3{ray.dir.x, ray.dir.y, ray.dir.z} * Real(rtc_ray.tfar);
    vertex.geometric_normal = normalize(Vector3{rtc_hit.Ng_x, rtc_hit.Ng_y, rtc_hit.Ng_z});
    vertex.shape_id = rtc_hit.geomID;
    vertex.primitive_id = rtc_hit.primID;
    const Shape &shape = scene.shapes[vertex.shape_id];
    vertex.material_id = get_material_id(shape);
    vertex.interior_medium_id = get_interior_medium_id(shape);
    vertex.exterior_medium_id = get_exterior_medium_id(shape);
    vertex.st = Vector2{rtc_hit.u, rtc_hit.v};

    ShadingInfo shading_info = compute_shading_info(scene.shapes[vertex.shape_id], vertex);
    vertex.shading_frame = shading_info.shading_frame;
    vertex.uv = shading_info.uv;
    vertex.mean_curvature = shading_info.mean_curvature;
    vertex.ray_radius = transfer(ray_diff, distance(ray.org, vertex.position));
    // vertex.ray_radius stores approximatedly dp/dx,
    // we get uv_screen_size (du/dx) using (dp/dx)/(dp/du)
    vertex.uv_screen_size = vertex.ray_radius / shading_info.inv_uv_size;
    vertex.time = ray.time;
    // Flip the geometry normal to the same direction as the shading normal
    if (dot(vertex.geometric_normal, vertex.shading_frame.n) < 0)
    {
        vertex.geometric_normal = -vertex.geometric_normal;
    }
    // get the user-defined payload
    {
        const Payload *payload = (const Payload *)user_context.payload;
        int n = payload->dist_t.size();
        for (int i = 0; i < n; ++i)
        {
            Real dist = payload->dist_t[i];
            if (dist <= rtc_ray.tfar && dist > rtc_ray.tnear) // change to include the equal case
            {
                ray.num_ghost_hits += 1;
                ray.ghost_shape_ids.push_back(payload->ghost_shape_ids[i]);
                ray.ghost_hit_positions.push_back(payload->ghost_hit_positions[i]);
                ray.ghost_hit_geom_normals.push_back(payload->ghost_hit_geom_normals[i]);
            }
        }
    }

    return vertex;
}

bool occluded(const Scene &scene, Ray &ray)
{
    IntersectContext user_context;
    InitIntersectionContext(&user_context);
    RTCRay rtc_ray;
    rtc_ray.org_x = (float)ray.org[0];
    rtc_ray.org_y = (float)ray.org[1];
    rtc_ray.org_z = (float)ray.org[2];
    rtc_ray.dir_x = (float)ray.dir[0];
    rtc_ray.dir_y = (float)ray.dir[1];
    rtc_ray.dir_z = (float)ray.dir[2];
    rtc_ray.tnear = (float)ray.tnear;
    rtc_ray.tfar = (float)ray.tfar;
    rtc_ray.mask = (unsigned int)(-1);
    rtc_ray.time = (float)ray.time; // time
    rtc_ray.flags = 0;
    // TODO: switch to rtcOccluded16
    rtcOccluded1(scene.embree_scene, &user_context.rtc_context, &rtc_ray);
    // get the user-defined payload

    {
        const Payload *payload = (const Payload *)user_context.payload;
        ray.num_ghost_hits = payload->num_ghost_hits;
    }
    return rtc_ray.tfar < 0;
}

Spectrum emission(const PathVertex &v,
                  const Vector3 &view_dir,
                  const Scene &scene)
{
    int light_id = get_area_light_id(scene.shapes[v.shape_id]);
    assert(light_id >= 0);
    const Light &light = scene.lights[light_id];
    return emission(light,
                    view_dir,
                    v.uv_screen_size,
                    PointAndNormal{v.position, v.geometric_normal},
                    scene);
}

/////////////for incremental rendering//////
std::optional<PathVertex> intersect_dynamic(const Scene &scene,
                                            Ray &ray,
                                            const RayDifferential &ray_diff)
{
    IntersectContext user_context; // there are rtc_context inside
    InitIntersectionContext(&user_context);
    RTCRayHit rtc_rayhit;
    RTCRay &rtc_ray = rtc_rayhit.ray;
    RTCHit &rtc_hit = rtc_rayhit.hit;
    rtc_ray = RTCRay{
        (float)ray.org.x, (float)ray.org.y, (float)ray.org.z,
        (float)ray.tnear,
        (float)ray.dir.x, (float)ray.dir.y, (float)ray.dir.z,
        (float)ray.time, // time
        (float)ray.tfar,
        (unsigned int)(-1), // mask
        0,                  // ray ID
        0                   // ray flags
    };
    rtc_hit = RTCHit{
        0, 0, 0, // Ng_x, Ng_y, Ng_z
        0,
        0,                        // u, v
        RTC_INVALID_GEOMETRY_ID,  // primitive ID
        RTC_INVALID_GEOMETRY_ID,  // geometry ID
        {RTC_INVALID_GEOMETRY_ID} // instance IDs
    };
    rtcIntersect1(scene.dynamic_scene, &user_context.rtc_context, &rtc_rayhit);
    if (rtc_hit.geomID == RTC_INVALID_GEOMETRY_ID)
    {
        return {};
    };
    assert(rtc_hit.geomID < scene.shapes.size());

    PathVertex vertex;
    vertex.position = Vector3{ray.org.x, ray.org.y, ray.org.z} +
                      Vector3{ray.dir.x, ray.dir.y, ray.dir.z} * Real(rtc_ray.tfar);
    vertex.geometric_normal = normalize(Vector3{rtc_hit.Ng_x, rtc_hit.Ng_y, rtc_hit.Ng_z});
    vertex.shape_id = rtc_hit.geomID;
    vertex.primitive_id = rtc_hit.primID;
    const Shape &shape = scene.shapes[vertex.shape_id];
    vertex.material_id = get_material_id(shape);
    vertex.interior_medium_id = get_interior_medium_id(shape);
    vertex.exterior_medium_id = get_exterior_medium_id(shape);
    vertex.st = Vector2{rtc_hit.u, rtc_hit.v};

    ShadingInfo shading_info = compute_shading_info(scene.shapes[vertex.shape_id], vertex);
    vertex.shading_frame = shading_info.shading_frame;
    vertex.uv = shading_info.uv;
    vertex.mean_curvature = shading_info.mean_curvature;
    vertex.ray_radius = transfer(ray_diff, distance(ray.org, vertex.position));
    // vertex.ray_radius stores approximatedly dp/dx,
    // we get uv_screen_size (du/dx) using (dp/dx)/(dp/du)
    vertex.uv_screen_size = vertex.ray_radius / shading_info.inv_uv_size;
    vertex.time = ray.time;
    // Flip the geometry normal to the same direction as the shading normal
    if (dot(vertex.geometric_normal, vertex.shading_frame.n) < 0)
    {
        vertex.geometric_normal = -vertex.geometric_normal;
    }

    return vertex;
}

bool occluded_dynamic(const Scene &scene, Ray &ray)
{
    IntersectContext user_context;
    InitIntersectionContext(&user_context);
    RTCRay rtc_ray;
    rtc_ray.org_x = (float)ray.org[0];
    rtc_ray.org_y = (float)ray.org[1];
    rtc_ray.org_z = (float)ray.org[2];
    rtc_ray.dir_x = (float)ray.dir[0];
    rtc_ray.dir_y = (float)ray.dir[1];
    rtc_ray.dir_z = (float)ray.dir[2];
    rtc_ray.tnear = (float)ray.tnear;
    rtc_ray.tfar = (float)ray.tfar;
    rtc_ray.mask = (unsigned int)(-1);
    rtc_ray.time = (float)ray.time; // time
    rtc_ray.flags = 0;
    // TODO: switch to rtcOccluded16
    rtcOccluded1(scene.dynamic_scene, &user_context.rtc_context, &rtc_ray);
    return rtc_ray.tfar < 0;
}

/**
 * @brief re-test the times of intersection with ghost objects
 * fill the ray with num_ghost_hits and return whether it ever hit a ghost
 * @param scene
 * @param ray
 * @return true
 * @return false
 */
bool intersect_ghost(const Scene &scene, Ray &ray)
{
    IntersectContext user_context; // there are rtc_context inside
    InitIntersectionContext(&user_context);
    RTCRayHit rtc_rayhit;
    RTCRay &rtc_ray = rtc_rayhit.ray;
    RTCHit &rtc_hit = rtc_rayhit.hit;
    rtc_ray = RTCRay{
        (float)ray.org.x, (float)ray.org.y, (float)ray.org.z,
        (float)ray.tnear,
        (float)ray.dir.x, (float)ray.dir.y, (float)ray.dir.z,
        (float)ray.time, // time
        (float)ray.tfar,
        (unsigned int)(-1), // mask
        0,                  // ray ID
        0                   // ray flags
    };
    rtc_hit = RTCHit{
        0, 0, 0, // Ng_x, Ng_y, Ng_z
        0,
        0,                        // u, v
        RTC_INVALID_GEOMETRY_ID,  // primitive ID
        RTC_INVALID_GEOMETRY_ID,  // geometry ID
        {RTC_INVALID_GEOMETRY_ID} // instance IDs
    };
    rtcIntersect1(scene.ghost_scene, &user_context.rtc_context, &rtc_rayhit);
    if (rtc_hit.geomID == RTC_INVALID_GEOMETRY_ID)
    {
        {
            const Payload *payload = (const Payload *)user_context.payload;
            int n = payload->dist_t.size();
            for (int i = 0; i < n; ++i)
            {
                Real dist = payload->dist_t[i];
                if (dist <= rtc_ray.tfar && dist > rtc_ray.tnear)
                {
                    ray.num_ghost_hits += 1;
                    ray.ghost_shape_ids.push_back(payload->ghost_shape_ids[i]);
                    ray.ghost_hit_positions.push_back(payload->ghost_hit_positions[i]);
                    ray.ghost_hit_geom_normals.push_back(payload->ghost_hit_geom_normals[i]);
                }
            }
        }
        return ray.num_ghost_hits > 0;
    };
    assert(rtc_hit.geomID < scene.shapes.size());
    {
        const Payload *payload = (const Payload *)user_context.payload;
        int n = payload->dist_t.size();
        for (int i = 0; i < n; ++i)
        {
            Real dist = payload->dist_t[i];
            if (dist <= rtc_ray.tfar && dist > rtc_ray.tnear)
            {
                ray.num_ghost_hits += 1;
                ray.ghost_shape_ids.push_back(payload->ghost_shape_ids[i]);
                ray.ghost_hit_positions.push_back(payload->ghost_hit_positions[i]);
                ray.ghost_hit_geom_normals.push_back(payload->ghost_hit_geom_normals[i]);
            }
        }
    }
    return ray.num_ghost_hits > 0;
}
