#include "scene.h"
#include "table_dist.h"
#include "transform.h"
Scene::Scene(const RTCDevice &embree_device,
             const Camera &camera,
             const std::vector<Material> &materials,
             const std::vector<Shape> &shapes,
             const std::vector<Light> &lights,
             std::vector<Dynamic> &in_dynamics, // not constant here
             const std::vector<Medium> &media,
             int envmap_light_id,
             const TexturePool &texture_pool,
             const RenderOptions &options,
             const std::string &output_filename) : embree_device(embree_device), camera(camera), materials(materials),
                                                   shapes(shapes), lights(lights), dynamics(in_dynamics), media(media),
                                                   envmap_light_id(envmap_light_id),
                                                   texture_pool(texture_pool), options(options),
                                                   output_filename(output_filename)
{
    // Register the geometry to Embree
    embree_scene = rtcNewScene(embree_device);
    dynamic_scene = rtcNewScene(embree_device);
    ghost_scene = rtcNewScene(embree_device);
    // We don't care about build time.
    rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
    rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_ROBUST);

    //////////// incremental rendering ///////////////////////
    if (dynamics.size() > 0)
    {
        rtcSetSceneBuildQuality(ghost_scene, RTC_BUILD_QUALITY_HIGH);
        rtcSetSceneFlags(ghost_scene, RTC_SCENE_FLAG_ROBUST);
        rtcSetSceneBuildQuality(dynamic_scene, RTC_BUILD_QUALITY_HIGH);
        rtcSetSceneFlags(dynamic_scene, RTC_SCENE_FLAG_ROBUST);
    }

    for (const Shape &shape : this->shapes)
    {
        if (is_dynamic(shape))
        {                                                       // if dynamic, register as instance geometry in embree
            Dynamic &dynamic = dynamics[get_dynamic_id(shape)]; // output
            uint32_t solid_dyn_geomID = register_embree_dynamic(shape, embree_device, embree_scene, 2, dynamic);
            dynamic.dynamic_geomID = solid_dyn_geomID;

            {
                init_sampling_dist(dynamic.base_shape);
                init_sampling_dist(dynamic.old_shape);
            }
            // register the dynamic shape again in the specific dynamic scene.
            uint32_t solid_dyn_geomID_scene1 = register_embree_dynamic(shape, embree_device, dynamic_scene, 2, dynamic);
        }
        else
        { // if static register as usual User Geometry
            uint32_t geomID = register_embree(shape, embree_device, embree_scene);
        }
    }
    for (const Shape &shape : this->shapes)
    {
        if (is_dynamic(shape))
        {
            // also add the ghost sibling
            Dynamic &dynamic = dynamics[get_dynamic_id(shape)];
            uint32_t ghost_dyn_geomID = register_embree_ghost_dynamic(shape, embree_device, embree_scene, ghost_scene, 2, dynamic);
            dynamic.ghost_geomID = ghost_dyn_geomID;
            ghost_shapes.push_back(shape);
        }
    }

    rtcCommitScene(dynamic_scene);
    rtcCommitScene(ghost_scene);
    rtcCommitScene(embree_scene);

    // TODO currently assume that the BBox of the scene is not expanded by the dynamic objects
    //  Get scene bounding box from Embree
    RTCBounds embree_bounds;
    rtcGetSceneBounds(embree_scene, &embree_bounds);
    Vector3 lb{embree_bounds.lower_x, embree_bounds.lower_y, embree_bounds.lower_z};
    Vector3 ub{embree_bounds.upper_x, embree_bounds.upper_y, embree_bounds.upper_z};
    bounds = BSphere{distance(ub, lb) / 2, (lb + ub) / Real(2)};

    // build shape & light sampling distributions if necessary
    // TODO: const_cast is a bit ugly...
    std::vector<Shape> &mod_shapes = const_cast<std::vector<Shape> &>(this->shapes);
    for (Shape &shape : mod_shapes)
    {
        init_sampling_dist(shape);
    }
    std::vector<Light> &mod_lights = const_cast<std::vector<Light> &>(this->lights);
    for (Light &light : mod_lights)
    {
        init_sampling_dist(light, *this);
    }

    // build a sampling distributino for all the lights
    std::vector<Real> power(this->lights.size());
    for (int i = 0; i < (int)this->lights.size(); i++)
    {
        power[i] = light_power(this->lights[i], *this);
    }
    light_dist = make_table_dist_1d(power);
}

Scene::Scene(const Scene &t) : embree_device(t.embree_device), embree_scene(t.embree_scene),
                               dynamic_scene(t.dynamic_scene), ghost_scene(t.ghost_scene),
                               camera(t.camera), materials(t.materials),
                               shapes(t.shapes), ghost_shapes(t.ghost_shapes), dynamics(t.dynamics),
                               lights(t.lights), media(t.media),
                               envmap_light_id(t.envmap_light_id),
                               texture_pool(t.texture_pool), bounds(t.bounds), options(t.options),
                               output_filename(t.output_filename), light_dist(t.light_dist)
{
    // embree_scene is managed by Embree using reference counting.
    // it's crucial for us to also maintain the reference count when copying scenes.
    rtcRetainScene(t.embree_scene);
    rtcRetainScene(t.dynamic_scene);
    rtcRetainScene(t.ghost_scene);
}

Scene::~Scene()
{
    // This decreses the reference count of embree_scene in Embree,
    // if it reaches zero, Embree will deallocate the scene.
    rtcReleaseScene(embree_scene);
    rtcReleaseScene(dynamic_scene);
    rtcReleaseScene(ghost_scene);
}

int sample_light(const Scene &scene, Real u)
{
    return sample(scene.light_dist, u);
}

Real light_pmf(const Scene &scene, int light_id)
{
    return pmf(scene.light_dist, light_id);
}

// adapted from Iliyan's VCM code

Spectrum connectToSensor(
    const PathVertex &vertex, // intput
    Vector3 &dirToSensor,     // output
    Real &cosAtSensor,        // output
    Vector2 &pRaster2,        // output
    Real &pdfSensorW,         // output
    const Scene &scene,
    int time)
{
    const Camera &camera = scene.camera;
    Vector3 sensorPos{0., 0., 0.};
    Spectrum value = sample_Wi(camera, vertex, sensorPos, dirToSensor, cosAtSensor, pRaster2, pdfSensorW);
    // test whether the connection to sensor gets occluded by other objects
    if (pdfSensorW == 0.0)
        return make_zero_spectrum();
    Real dist = distance(sensorPos, vertex.position);
    Ray shadow_ray{vertex.position, dirToSensor, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * dist, time};
    if (occluded(scene, shadow_ray))
        return make_zero_spectrum();
    return value;
}
// the only difference is also to pass out a vertex
Spectrum connectToSensorVertex(
    const PathVertex &vertex, // intput
    PathVertex &sampled,      // output
    Vector3 &dirToSensor,     // output
    Vector2 &pRaster2,        // output
    const Scene &scene,
    int time,
    bool testGhost) // whether to test the number of ghost hits for the shadow rays
{
    const Camera &camera = scene.camera;
    Vector3 sensorPos{0., 0., 0.};
    Real cosAtSensor=0., pdfSensorW=0., pdfAreaVertex=0.;
    Spectrum value = sample_Wi(camera, vertex, sensorPos, dirToSensor, cosAtSensor, pRaster2, pdfSensorW, pdfAreaVertex, time);
    if (pdfSensorW == 0.0)
        return make_zero_spectrum();

    Real dist = distance(sensorPos, vertex.position);
    {
        sampled.type = VertexType::Camera;
        sampled.time = time;
        sampled.position = xform_point(scene.camera.cam_to_world, Vector3{0.0, 0.0, 0.0});
        sampled.geometric_normal = normalize(xform_normal(scene.camera.cam_to_world, Vector3{0.0, 0.0, 1.0}));
        sampled.transport.contrib = value; // TODO confirm or ?make_const_spectrum(1.0);
        sampled.transport.pathLength = 1;
        sampled.transport.pdfFwd = pdfAreaVertex;
        sampled.transport.pdfDir = pdfSensorW;
        sampled.transport.pdfOriginA = pdfAreaVertex;//pdfA of vertex
    }

    Ray shadow_ray{vertex.position, dirToSensor, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * dist, time};
    if (!testGhost)
    {
        if (occluded(scene, shadow_ray))
            return make_zero_spectrum();
    }
    else
    { // use intersection rays for the test
        std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, RayDifferential{});
        if (vertex_)
            return make_zero_spectrum(); // hit something
        sampled.hit_ghost_times = shadow_ray.num_ghost_hits;
        sampled.ghost_ids = shadow_ray.ghost_shape_ids;
        sampled.ghost_hit_positions = shadow_ray.ghost_hit_positions;
        sampled.ghost_hit_geom_normals = shadow_ray.ghost_hit_geom_normals;
    }

    return value;
}

int sample_dynamic_shape(const Scene &scene, Real u)
{
    return int(u*scene.dynamics.size());
}
/// The probability mass function of the sampling procedure above.
Real dynamic_pmf(const Scene &scene, int dynamic_id)
{
    return Real(1.0) / scene.dynamics.size();
}
