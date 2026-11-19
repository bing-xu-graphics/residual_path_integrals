#pragma once

#include "scene.h"
#include "pcg.h"
#include "../utils.h"
#include "portal_utils.h"
#include "portal_utils_ghost.h"
#include "portal_utils_mix.h"
#include "incremental_basepath.h"

#define BETTER_MIS_GHOST
/**
 * @brief generate two subpaths passing a ghost dynamic object.
 * Make the starting point a ghost dymmy
 * 1) sample a ghost point on the ghost bbox;
 * 2) sample another point on the ghost bbox and this produces a dir; test the intersection with the actual ghost
 *     without loss of generality, make this dir2Light
 * 3) dir2Sensor = -dir2Light
 *
 * @param lightPath
 * @param sensorPath
 * @param rng
 * @param scene
 * @param time
 * @return true
 * @return false
 */
bool genGhostSubpaths(std::vector<PathVertex> &lightPath,
                      std::vector<PathVertex> &sensorPath,
                      pcg32_state &rng,
                      const Scene &scene,
                      int time = 0)
{
    int maxDepth = scene.options.max_depth;
    PathVertex vGhost;
    Vector3 dir2Light, dir2Sensor;
    bool success = sampleGhost(vGhost, dir2Light, rng, scene, time);
    if (!success)
    {
        return false;
    }
    Real pdfDir2Light = vGhost.transport.pdfDir, pdfDir2Sensor = 1.0;
    // the other dir is already defined
    dir2Sensor = -dir2Light;
    lightPath.push_back(vGhost);

    RayDifferential ray_diff;
    int k = maxDepth - 1;
    randomWalk(lightPath, make_const_spectrum(1.0), 1.0, dir2Light, ray_diff, k, 1, TransportDirection::TO_LIGHT, rng, scene, time);
    const PathVertex &vertex1 = lightPath[1]; // the vertex at S end of the ghost edge
    lightPath[1].hit_ghost_times = 0;         // TODO check
    vGhost.transport.pdfDir = pdfDir2Sensor;  // set to 1.0
    vGhost.transport.pdfFwd = pdfDir2Sensor;
    sensorPath.push_back(vertex1); // if change this; change the lightPath[1].hit_ghost_times = 0; together
    randomWalk(sensorPath, make_const_spectrum(1.0), 1.0, dir2Sensor, ray_diff, k, 1, TransportDirection::TO_VIEW, rng, scene, time);
    return true;
}

/**
 * @brief for general case where the two intersections from a ghost edge are already found
 * s>=2 && t>=2
 * @param lightPath
 * @param sensorPath
 * @param vDynamicS
 * @param vDynamicT
 * @param s
 * @param t
 * @param screenPos
 * @param rng
 * @param scene
 * @param numSamples
 * @param timee
 * @return Spectrum
 */
Spectrum render_general_ghost_path(const std::vector<PathVertex> &lightPath,
                                   const std::vector<PathVertex> &sensorPath,
                                   const PathVertex &vGhostS,
                                   const PathVertex &lastHitLight,
                                   int s, int t,
                                   const DSubPathState& lightPathState,
                                   const DSubPathState& sensorPathState,
                                   Vector2 &screenPos,
                                   pcg32_state &rng,
                                   const Scene &scene,
                                   const int &numSamples,
                                   int time)
{ // the general case where s >=2 and t>=2
    assert(s >= 2 && t >= 2);
    const PathVertex &vLight = lightPath[s - 1];
    const PathVertex &vSensor = sensorPath[t - 1];
    // got the two intersections from this ghost edge
    const PathVertex &prevLight = lightPath[s - 2];
    const PathVertex &prevSensor = sensorPath[t - 2];
    assert(vSensor.type != VertexType::Light);
    
    Spectrum Llight = make_zero_spectrum(), Wsensor = make_zero_spectrum();
    PathVertex sampledVLight;
    int nLight = lightPath.size();

    { // connect to light at s side
        Vector3 dir2Light;
        Llight = directIllumination(vLight, dir2Light, sampledVLight, rng, scene, time, true); // Le * G / pdfLight
        if (average(Llight) == 0.0)
            return make_zero_spectrum();
        const Material &mat = scene.materials[vLight.material_id];
        Vector3 dir2View = normalize(prevLight.position - vLight.position);
        Spectrum f = eval(mat, dir2View, dir2Light, vLight, scene.texture_pool) / fabs(dot(dir2Light, vLight.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
        Llight *= vLight.transport.contrib * f;
    }

    Vector2 pRaster2{0.0, 0.0};
    PathVertex sampledVSensor;
    { // connect to sensor at t side
        Vector3 dir2Sensor;
        Spectrum We = connectToSensorVertex(vSensor, sampledVSensor, dir2Sensor, pRaster2, scene, time, true);
        if (average(We) == 0.0)
            return make_zero_spectrum();
        Wsensor = We / (Real)numSamples;
        Vector3 dir2Light = normalize(prevSensor.position - vSensor.position);
        assert(vSensor.material_id >= 0);
        Spectrum f = eval(scene.materials[vSensor.material_id], dir2Light, dir2Sensor, vSensor, scene.texture_pool) / fabs(dot(dir2Sensor, vSensor.geometric_normal));
        Wsensor *= vSensor.transport.contrib * f;
    }
    { // handle the ghost edge in the middle
        // this is like one NULL node
        Spectrum radiance = Llight * Wsensor / vGhostS.transport.pdfDir * vGhostS.transport.contrib;
        {
            Real cosTheta1 = fabs(dot(normalize(lightPath[1].position - vGhostS.position), vGhostS.geometric_normal)); /// distance_squared(lightPath[1].position, vGhostS.position);
            radiance *= cosTheta1;
        }

        // only need G on one side of the ghost
        screenPos = pRaster2;

#if defined(BETTER_MIS_GHOST)
        Real base_pdf;
        Real weight = computeMixMISWeight_better_ghost(vGhostS, vLight, vSensor, sampledVLight, sampledVSensor, lightPath, sensorPath,
                s, t, lightPathState, sensorPathState, SamplingTechnique::GHOST_TWO_ENDS, scene, time, base_pdf);
#elif defined(NAIVE_MIS_MIX_GHOST)
        Real weight = compute_mis_ghost_part(sampledVLight, vLight, vSensor, sampledVSensor, SamplingTechnique::GHOST_TWO_ENDS);
#else
        Real weight = computeGhostPortalMISWeight(&sampledVLight, &vLight, &vSensor, &sampledVSensor);
#endif
        radiance *= weight;
        assert(average(radiance) == average(radiance));
        assert(!isinf(average(radiance)));
        return radiance;
    }
    return make_zero_spectrum();
}

void render_ghost_path_1_from_light(const PathVertex &vGhostS,
                                    const Spectrum &Llight,
                                    const std::vector<PathVertex> &lightPath,
                                    const std::vector<PathVertex> &sensorPath,
                                    pcg32_state &rng,
                                    const Scene &scene,
                                    const int &numSamples,
                                    const int time,
                                    Image3 &img)
{
    int nSensor = sensorPath.size();
    if (nSensor < 2)
        return;
    // prepare the accumulated pdfFwd and pdfRevs
    DSubPathState pathState(nSensor);
    cumprod_pdfs_except_origin(sensorPath, pathState);
    for (int t = 2; t <= nSensor; ++t)
    {
        const PathVertex &vertex = sensorPath[t - 1];
        Vector3 dirToSensor{0., 0., 0.};
        Vector2 pRaster2{0.0, 0.0};
        PathVertex sampledVSensor;
        Spectrum We = connectToSensorVertex(vertex, sampledVSensor, dirToSensor, pRaster2, scene, time, true);
        if (average(We) == 0.0)
            continue;
        ////////////// connect each path vertex to the sensor
        const PathVertex &prev = sensorPath[t - 2];
        Vector3 dirLight = normalize(prev.position - vertex.position);
        Spectrum f = eval(scene.materials[vertex.material_id], dirLight, dirToSensor, vertex, scene.texture_pool) / fabs(dot(dirToSensor, vertex.geometric_normal));
        Spectrum radiance = Llight * vGhostS.transport.contrib / (Real)numSamples * We * vertex.transport.contrib * f;
#if defined(BETTER_MIS_GHOST)
        Real base_pdf;
        Real weight = computeMixMISWeight_better_ghost(vGhostS, vertex, vertex, sensorPath[0], sampledVSensor, lightPath,sensorPath, 0, t, 
        NULL, pathState, SamplingTechnique::GHOST_FROM_LIGHT, scene, time, base_pdf);
#elif defined(NAIVE_MIS_MIX_GHOST)
        Real weight = compute_mis_ghost_part({}, {}, vertex, sampledVSensor, SamplingTechnique::GHOST_FROM_LIGHT);
#else
        Real weight = computeGhostPortalMISWeight(nullptr, nullptr, &vertex, &sampledVSensor);
#endif
        radiance *= weight;
        assert(average(radiance) == average(radiance));
        assert(!isinf(average(radiance)));
        splatSample(img, pRaster2, radiance);
    }
}

void render_ghost_path2_from_sensor(const PathVertex &vGhostS,
                                    const Spectrum &We,
                                    const Vector2 &pRaster2, // passed from the base
                                    std::vector<PathVertex> &lightPath,
                                    const std::vector<PathVertex> &sensorPath,
                                    pcg32_state &rng,
                                    const Scene &scene,
                                    const int &numSamples,
                                    const int time,
                                    Image3 &img)
{
    PathVertex lastHitLight;
    if(lightPath.back().type == VertexType::Light){
        lastHitLight = lightPath.back();
        //then remove it from the list 
        lightPath.pop_back();
        assert(lastHitLight.type == VertexType::Light);
    }

    int nLight = lightPath.size();
    if (nLight < 2) // TODO check this!!
        return;
    // prepare the accumulated pdfFwd and pdfRevs
    DSubPathState pathState(nLight);
    cumprod_pdfs_except_origin(lightPath, pathState);

    for (int s = 2; s <= nLight; ++s)
    {
        const PathVertex &vertex = lightPath[s - 1];
        // do NEE or directly hit light here
        PathVertex sampledVLight;
        Spectrum Llight = make_zero_spectrum();

        {
            Vector3 dir2Light;
            Llight = directIllumination(vertex, dir2Light, sampledVLight, rng, scene, time, true); // L*G/pdfLight
            if (average(Llight) == 0.0)
                continue;
            assert(vertex.material_id >= 0);
            const PathVertex &prev = lightPath[s - 2];
            Vector3 dir2Sensor = normalize(prev.position - vertex.position);
            Spectrum f = eval(scene.materials[vertex.material_id], dir2Light, dir2Sensor, vertex, scene.texture_pool);
            f /= fabs(dot(dir2Sensor, vertex.geometric_normal));
            Llight *= vertex.transport.contrib * f;
        }

        Spectrum radiance = We * Llight * vGhostS.transport.contrib / (Real)numSamples;
        if (average(radiance) == 0.0)
            continue;
#if defined(BETTER_MIS_GHOST)
        Real base_pdf;
        Real weight = computeMixMISWeight_better_ghost(vGhostS, vertex, vertex, sampledVLight, lightPath[0], lightPath, sensorPath, s, 0, 
            pathState, NULL, SamplingTechnique::GHOST_FROM_SENSOR, scene, time, base_pdf);
#elif defined(NAIVE_MIS_MIX_GHOST)
        Real weight = compute_mis_ghost_part(sampledVLight, vertex, {}, {}, SamplingTechnique::GHOST_FROM_SENSOR);
#else
        Real weight = computeGhostPortalMISWeight(&sampledVLight, &vertex, nullptr, nullptr);
#endif
        radiance *= weight;
        assert(average(radiance) == average(radiance));
        assert(!isinf(average(radiance)));
        splatSample(img, pRaster2, radiance);
    }
}

void render_ghost_path_3_two_ends_noredundant(const PathVertex &vGhostS,
                                  const Vector3 &ghost_dir2Light,
                                  std::vector<PathVertex> &lightPath,
                                  const std::vector<PathVertex> &sensorPath,
                                  pcg32_state &rng,
                                  const Scene &scene,
                                  const int &numSamples,
                                  const int baseTimeStamp,
                                  Image3 &img,
                                  bool negative = false)
{
    int nLight = lightPath.size(), nSensor = sensorPath.size();
    if (nLight < 2 || nSensor < 2)
        return;

    DSubPathState lightPathState(nLight), sensorPathState(nSensor);
    cumprod_pdfs_except_origin(lightPath, lightPathState);
    cumprod_pdfs_except_origin(sensorPath, sensorPathState);

    // to accelarate the performance, unfold the nested loop
    std::vector<Spectrum> Llight_cache;       // cache the radiance for the light subpath, at index of 0, it holds information from s==2
    std::vector<Spectrum> Wsensor_cache;      // cache the radiance for the sensor subpath
    std::vector<Vector2> base_pRaster2_cache; // cache the raster position for the sensor subpath
    // lets also record the sampledVLight and sampledVsensor for the base path for later use, we can only record the hit_times, but let's just settle for now
    std::vector<PathVertex> sampledVLight_cache;
    std::vector<PathVertex> sampledVSensor_cache;

    // 2) to fill in the cache for both side of the subpath seperately
    for (int s = 2; s <= nLight; ++s)
    {
        // 1) render the base path
        const PathVertex &vLight = lightPath[s - 1];
        const PathVertex &prevLight = lightPath[s - 2];
        Spectrum Llight = make_zero_spectrum();
        PathVertex base_sampledVLight;

        { // connect to light at s side
            Vector3 dir2Light;
            Llight = directIllumination(vLight, dir2Light, base_sampledVLight, rng, scene, baseTimeStamp, true); // Le * G / pdfLight
            const Material &mat = scene.materials[vLight.material_id];
            Vector3 dir2View = normalize(prevLight.position - vLight.position);
            Spectrum f = eval(mat, dir2View, dir2Light, vLight, scene.texture_pool) / fabs(dot(dir2Light, vLight.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
            Llight *= vLight.transport.contrib * f;
        }
        Llight_cache.push_back(Llight);
        sampledVLight_cache.push_back(base_sampledVLight);
    }

    for (int t = 2; t <= nSensor; ++t)
    {
       // 1) render the base sensor subpath
        const PathVertex &vSensor = sensorPath[t - 1];
        const PathVertex &prevSensor = sensorPath[t - 2];

        Spectrum Wsensor = make_zero_spectrum();
        Vector3 base_dir2Sensor;
        Vector2 base_pRaster2{0.0, 0.0};
        PathVertex base_sampledVSensor;
        if (vSensor.type == VertexType::Light)
        {
            Wsensor_cache.push_back(Wsensor);
            base_pRaster2_cache.push_back(base_pRaster2);        // put in dummy to hold the place
            sampledVSensor_cache.push_back(base_sampledVSensor); // put in dummy to hold the place
            continue;                                         // TODO do this later
        }
        Spectrum We = connectToSensorVertex(vSensor, base_sampledVSensor, base_dir2Sensor, base_pRaster2, scene, baseTimeStamp, true);
        if (average(We) == 0.0)
        {
            Wsensor_cache.push_back(Wsensor);
            base_pRaster2_cache.push_back(base_pRaster2);        // put in dummy to hold the place
            sampledVSensor_cache.push_back(base_sampledVSensor); // put in dummy to hold the place
            continue;                                         // TODO do this later
        }

        Vector3 dir2Light = normalize(prevSensor.position - vSensor.position);
        Spectrum f_base_T = eval(scene.materials[vSensor.material_id], dir2Light, base_dir2Sensor, vSensor, scene.texture_pool) / fabs(dot(base_dir2Sensor, vSensor.geometric_normal));
        Wsensor = We / (Real)numSamples * vSensor.transport.contrib * f_base_T;

        Wsensor_cache.push_back(Wsensor);
        base_pRaster2_cache.push_back(base_pRaster2);
        sampledVSensor_cache.push_back(base_sampledVSensor);

    }

    Real cosTheta_atGhost = fabs(dot(normalize(lightPath[1].position - vGhostS.position), vGhostS.geometric_normal)); /// distance_squared(lightPath[1].position, vGhostS.position)
    
    // 3) go over the actually combinations (s, t) of the path and compose the radiance there.
    for (int s = 2; s <= nLight; ++s) /// num of vertices now;
    {
        for (int t = 2; t <= nSensor; ++t)
        {
            int depth = s + t - 1; // get rid of one depth
            if (depth < 1 || depth > scene.options.max_depth)
                continue;
            // render the base path
            
            const PathVertex &vLight = lightPath[s - 1];
            const PathVertex &vSensor = sensorPath[t - 1];
            
            int s_cache = s-2;
            int t_cache = t-2;
            // render the base path radiance
            const Spectrum &base_Llight = Llight_cache[s_cache];
            const Spectrum &base_Wsensor = Wsensor_cache[t_cache];
            const PathVertex &base_sampledVLight = sampledVLight_cache[s_cache];
            const PathVertex &base_sampledVSensor = sampledVSensor_cache[t_cache];
            if (average(base_Llight) == 0.0 || average(base_Wsensor) == 0.0)
                continue;

            Spectrum baseRadiance = base_Llight * base_Wsensor / vGhostS.transport.pdfDir * vGhostS.transport.contrib;
            if (average(baseRadiance) < 1e-10)
                continue;
            Vector2 base_pRaster2 = base_pRaster2_cache[t_cache];
            baseRadiance *= cosTheta_atGhost;
            // calculate the MIS weight
            Real base_pdf;
            Real base_weight = computeMixMISWeight_better_ghost(vGhostS, vLight, vSensor, 
                        base_sampledVLight, base_sampledVSensor, lightPath, sensorPath,
                        s, t, lightPathState, sensorPathState, SamplingTechnique::GHOST_TWO_ENDS, scene, baseTimeStamp, base_pdf);

            baseRadiance *= base_weight;
           
            splatSample(img, base_pRaster2, baseRadiance);
        }
    }
}

void render_ghost_path_3_two_ends(const PathVertex &vGhostS,
                                  const Vector3 &ghost_dir2Light,
                                  std::vector<PathVertex> &lightPath,
                                  const std::vector<PathVertex> &sensorPath,
                                  pcg32_state &rng,
                                  const Scene &scene,
                                  const int &numSamples,
                                  const int time,
                                  Image3 &img,
                                  bool negative = false)
{
    PathVertex lastHitLight;
    if(lightPath.back().type == VertexType::Light){
        lastHitLight = lightPath.back();
        //then remove it from the list 
        lightPath.pop_back();
        assert(lastHitLight.type == VertexType::Light);
    }

    int nLight = lightPath.size(), nSensor = sensorPath.size();
    if (nLight < 2 || nSensor < 2)
        return;

    DSubPathState lightPathState(nLight), sensorPathState(nSensor);
    cumprod_pdfs_except_origin(lightPath, lightPathState);
    cumprod_pdfs_except_origin(sensorPath, sensorPathState);

    for (int s = 2; s <= nLight; ++s) /// num of vertices now;
    {
        for (int t = 2; t <= nSensor; ++t)
        {
            int depth = s + t - 1; // get rid of one depth
            if (depth < 1 || depth > scene.options.max_depth)
                continue;

            Vector2 screenPos{0., 0.}; // the splatting screen position and the radiance for the path
            Spectrum radiance = render_general_ghost_path(lightPath, sensorPath, vGhostS, lastHitLight, s, t, 
                                                            lightPathState, sensorPathState,
                                                          screenPos, rng, scene, numSamples, time);
            
            if(negative) splatSample(img, screenPos, -radiance);
            else splatSample(img, screenPos, radiance);
        }
    }
}

Spectrum render_pt_frame_selector_ghost(const Scene &scene,
                                  int x, int y,              /* pixel coordinates */
                                  RAY_VISIBLE_MASK ray_mask, // only render the rays intersected with dynamic or pass the ghost
                                  pcg32_state &rng,
                                  int time)
{
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos, time);
    RayDifferential ray_diff = init_ray_differential(w, h);
    PathVertex sampledVSensor;
    {
        sampledVSensor.type = VertexType::Camera;
        sampledVSensor.time = time;
        sampledVSensor.position = xform_point(scene.camera.cam_to_world, Vector3{0.0, 0.0, 0.0});
        sampledVSensor.geometric_normal = normalize(xform_normal(scene.camera.cam_to_world, Vector3{0.0, 0.0, 1.0}));
        sampledVSensor.transport.contrib = make_const_spectrum(1.0); // TODO confirm or ?make_const_spectrum(1.0);
        sampledVSensor.transport.pathLength = 0;
        sampledVSensor.transport.pdfOriginA = 1.0;

        const Real cosAtSensor = fabs(dot(sampledVSensor.geometric_normal, ray.dir));
        // Compute pdf conversion factor from area on image plane to solid angle on ray
        Real imagePointToSensorDist = scene.camera.image_plane_distance / cosAtSensor;
        Real imageToSolidAdngleFactor = imagePointToSensorDist * imagePointToSensorDist / cosAtSensor;
        // We put the virtual image plane at such a distance from the camera origin
        // that the pixel area is one and thus the image plane sampling pdf is 1.
        // The solid angle ray pdf is then equal to the conversion factor from
        // image plane area density to ray solid angle density
        const Real cameraPdfW = imageToSolidAdngleFactor;
        sampledVSensor.transport.pdfDir = cameraPdfW;
        sampledVSensor.transport.pdfFwd = sampledVSensor.transport.pdfRev = 1.0;
    }

    std::vector<PathVertex> lightPath;//, sensorPath; // record the vertices for misWeight evaluation;
    lightPath.push_back(sampledVSensor);
    randomWalk_easy(lightPath, make_const_spectrum(1.0), 1.0, ray.dir, ray_diff, scene.options.max_depth, 1, TransportDirection::TO_LIGHT, rng, scene, time);

    if(lightPath.size()<2) return make_zero_spectrum();
    lightPath.erase(lightPath.begin()); //get rid of sampledVSensor
    
    PathVertex lastHitLight;
    if(lightPath.back().type == VertexType::Light){
        lastHitLight = lightPath.back();
        //then remove it from the list 
        lightPath.pop_back();
        assert(lastHitLight.type == VertexType::Light);
    }

    int nLight = lightPath.size();
    // prepare the accumulated pdfFwd and pdfRevs
    DSubPathState pathState(nLight);
    lightPath[0].transport.pdfFwd = 1.0;
    cumprod_pdfs(lightPath, pathState);
    {
        Real pdfAPrimaryHit = sampledVSensor.transport.pdfDir 
        * fabs(dot(normalize(sampledVSensor.position - lightPath[0].position), lightPath[0].geometric_normal))
        / distance_squared(sampledVSensor.position, lightPath[0].position); // pdfA of sampling the first vertex
        sampledVSensor.transport.pdfFwd = sampledVSensor.transport.pdfOriginA = pdfAPrimaryHit;
    }
    //render the path//////////////// only show the contribution selected//////////
    bool isSelected = false; // to mark whether the current ray's contribution is selected i.e. whether satisfy the mask selector
    Spectrum radiance_selected = make_zero_spectrum();
    for (int s = 1; s <= nLight; ++s)
    {
        const PathVertex &vertex = lightPath[s - 1];
        if((ray_mask & RAY_VISIBLE_DYNAMIC) && is_vertex_dynamic(vertex, scene))
            isSelected = true;
        if((ray_mask & RAY_VISIBLE_GHOST) && vertex.hit_ghost_times > 0)
            isSelected = true;

        // next event estimation TODO test whether this edge hit ghost objects
        PathVertex sampledVLight;
        Spectrum Llight = make_zero_spectrum();

        {
            Vector3 dir_light;
            Llight = directIllumination(vertex, dir_light, sampledVLight, rng, scene, time, true);
            
            if (average(Llight) == 0.0)
                continue;
            const Material &mat = scene.materials[vertex.material_id];
            Vector3 dir_view = (s==1)? -ray.dir : normalize(lightPath[s-2].position - vertex.position);
            Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool) / fabs(dot(dir_light, vertex.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
            Llight *= vertex.transport.contrib * f; 
        }

        if((!isSelected) && !(sampledVLight.hit_ghost_times>0 && (ray_mask & RAY_VISIBLE_GHOST))) continue;
        
#if defined(BETTER_MIS_GHOST)
        Real base_pdf;
        Real weight = computeMixMISWeight_better_ghost({}, vertex, sampledVSensor, sampledVLight, sampledVSensor, lightPath, {}, s, 0, 
            pathState, NULL, SamplingTechnique::UNI_PT, scene, time, base_pdf);
#else 
        Real weight = 1.0 /(vertex.hit_ghost_times + sampledVLight.hit_ghost_times + 1.0);
#endif
        radiance_selected += weight * Llight;
    }

    return radiance_selected;
}
