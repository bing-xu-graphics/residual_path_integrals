

#pragma once

#include "scene.h"
#include "pcg.h"
#include "parallel.h"
#include "progress_reporter.h"
#include "transform.h"
#include "utils.h"

/**
 * @brief base path for dynamic part 
 */

Spectrum prepare_dynamic_path_1_from_light(PathVertex &vDynamic,
                                           std::vector<PathVertex> &lightPath,
                                           std::vector<PathVertex> &sensorPath,
                                           pcg32_state &rng,
                                           const Scene &scene,
                                           const int &numSamples,
                                           int time) // output
{
    PathVertex sampledVLight;
    Vector3 dynamic_dir2Light;
    Spectrum Llight = directIllumination(vDynamic, dynamic_dir2Light, sampledVLight, rng, scene, time, true); // Le/pdfLight * G_to_ghost??  theG_to_ghost is not right inside
    if (average(Llight) == 0.0)
        return Llight;
    // sample the other direction and then do light tracing
    Vector3 dynamic_dir2Sensor;
    Real pdfDir2Sensor;
    bool isOtherDir = sampleScattering(vDynamic, dynamic_dir2Light, dynamic_dir2Sensor, pdfDir2Sensor, rng, scene);
    if (!isOtherDir)
    {
        return make_zero_spectrum();
    }
    lightPath.push_back(vDynamic);
    lightPath.push_back(sampledVLight);
    sensorPath.push_back(vDynamic);
    sensorPath[0].transport.pdfDir = pdfDir2Sensor;
    sensorPath[0].transport.pdfFwd = 1.0; // TODO temporarilly
    sensorPath[0].transport.pdfRev = 1.0;
    
    RayDifferential ray_diff;

    randomWalk_easy(sensorPath, make_const_spectrum(1.0), pdfDir2Sensor, dynamic_dir2Sensor,
               ray_diff, scene.options.max_depth-1, 1, TransportDirection::TO_VIEW, rng, scene, time);
    return Llight;
}

Spectrum prepare_dynamic_path_2_from_sensor(PathVertex &vDynamic,
                                            std::vector<PathVertex> &lightPath,
                                            std::vector<PathVertex> &sensorPath,
                                            pcg32_state &rng,
                                            const Scene &scene,
                                            const int &numSamples,
                                            int time, Vector2 &o_pRaster2)
{
    // connect to camera first
    Vector3 dynamic_dir2Sensor;
    PathVertex sampledVSensor;
    Spectrum We = connectToSensorVertex(vDynamic, sampledVSensor, dynamic_dir2Sensor, o_pRaster2, scene, time, true);
    if (average(We) == 0.0)
        return We;
    // sample the other direction and then do path tracing
    Vector3 dynamic_dir2Light;
    Real pdfDir2Light;
    bool isOtherDir = sampleScattering(vDynamic, dynamic_dir2Sensor, dynamic_dir2Light, pdfDir2Light, rng, scene);
    if (!isOtherDir)
    {
        return make_zero_spectrum();
    }
    sensorPath.push_back(vDynamic);
    sensorPath.push_back(sampledVSensor);
    lightPath.push_back(vDynamic);
    lightPath[0].transport.pdfDir = pdfDir2Light;
    lightPath[0].transport.pdfFwd = 1.0;
    lightPath[0].transport.pdfRev = 1.0;
    RayDifferential ray_diff;
    randomWalk_easy(lightPath, make_const_spectrum(1.0), pdfDir2Light, dynamic_dir2Light,
               ray_diff, scene.options.max_depth - 1, 1, TransportDirection::TO_LIGHT, rng, scene, time);
    
    return We;
}

void prepare_dynamic_path_3_two_ends(
    PathVertex &vDynamic,
    const Vector3 &dir2Light,
    std::vector<PathVertex> &lightPath,
    std::vector<PathVertex> &sensorPath,
    pcg32_state &rng,
    const Scene &scene,
    const int &numSamples,
    int time)
{
    int k = scene.options.max_depth;

    Vector3 dir2Sensor;
    Real pdfDir2Sensor;
    bool isOtherDir = sampleScattering(vDynamic, dir2Light, dir2Sensor, pdfDir2Sensor, rng, scene);
    if (!isOtherDir)
    {
        return;
    }
    lightPath.push_back(vDynamic);
    sensorPath.push_back(vDynamic);
    sensorPath[0].transport.pdfDir = pdfDir2Sensor;
    sensorPath[0].transport.pdfFwd = 1.0;
    RayDifferential ray_diff;
    pcg32_state rng_cache = rng;
    randomWalk_easy(lightPath, make_const_spectrum(1.0), vDynamic.transport.pdfDir, dir2Light, ray_diff, k-1, 1, TransportDirection::TO_LIGHT, rng, scene, time);
    rng = rng_cache;
    randomWalk_easy(sensorPath, make_const_spectrum(1.0), pdfDir2Sensor, dir2Sensor, ray_diff, k-1, 1, TransportDirection::TO_VIEW, rng, scene, time);
}

/////////////////////////////////////////////////////////////////////////////////
//base path for the ghost part

// fill in only the sensorPath
Spectrum prepare_ghost_path_1_from_light(
    const PathVertex &vGhostS,
    std::vector<PathVertex> &lightPath,
    std::vector<PathVertex> &sensorPath,
    pcg32_state &rng,
    const Scene &scene,
    const int &numSamples,
    int time)
{

    int maxDepth = scene.options.max_depth;
    // re-sample the light source
    PathVertex sampledVLight;
    Vector3 ghotst_dir2Light;
    Spectrum Llight = directIllumination(vGhostS, ghotst_dir2Light, sampledVLight, rng, scene, time, false); // Le/pdfLight * G_to_ghost??  theG_to_ghost is not right inside
    if (average(Llight) == 0.0)
        return Llight;
    Vector3 ghost_dir2Sensor = -ghotst_dir2Light;
    // do light tracing.
    sensorPath.push_back(sampledVLight);
    RayDifferential ray_diff; // dummy_variable

    randomWalk_easy(sensorPath, make_const_spectrum(1.0), 1.0, ghost_dir2Sensor,
            ray_diff, maxDepth, 1, TransportDirection::TO_VIEW, rng, scene, time);
    if(sensorPath.size()>1) sensorPath[1].transport.pdfFwd = 1.0;
    sensorPath[0].transport.pdfFwd = 1.0;
    sensorPath[0].transport.pdfRev = 1.0;
    assert(lightPath.empty());
    return Llight;
}

// fill in only the light path
Spectrum prepare_ghost_path_2_from_sensor(
    const PathVertex &vGhostS,
    std::vector<PathVertex> &lightPath,
    std::vector<PathVertex> &sensorPath,
    pcg32_state &rng,
    const Scene &scene,
    const int &numSamples,
    int time,
    Vector2 &o_pRaster2)
{
    int maxDepth = scene.options.max_depth;
    Vector3 ghost_dir2Sensor;
    PathVertex sampledVSensor;
    Spectrum We = connectToSensorVertex(vGhostS, sampledVSensor, ghost_dir2Sensor, o_pRaster2, scene, time, false);
    if (average(We) == 0.0)
        return We;
    Vector3 ghost_dir2Light = -ghost_dir2Sensor;
    lightPath.push_back(sampledVSensor);
    RayDifferential ray_diff; // dummy_variable
    randomWalk_easy(lightPath, make_const_spectrum(1.0), 1.0, ghost_dir2Light,
            ray_diff, maxDepth, 1, TransportDirection::TO_LIGHT, rng, scene, time);
    if(lightPath.size()>1) lightPath[1].transport.pdfFwd = 1.0;
    lightPath[0].transport.pdfFwd = 1.0;//TODO put the pdfFwd for sensor somewhere else
    lightPath[0].transport.pdfRev = 1.0;
    assert(sensorPath.empty());
    return We;
}
/*
    lightPath = [vGhostS, v1, v2, ...,]
    sensorPath = [v1, v2, ...] //the number of hits with ghost of vector(v1, v2) is stored in v2.
*/
void prepare_ghost_path_3_two_ends(
    const PathVertex &vGhostS,
    const Vector3 &dir2Light,
    std::vector<PathVertex> &lightPath,
    std::vector<PathVertex> &sensorPath,
    pcg32_state &rng,
    const Scene &scene,
    const int &numSamples,
    int time)
{
    int k = scene.options.max_depth;
    lightPath.push_back(vGhostS);
    RayDifferential ray_diff;
    randomWalk_easy(lightPath, make_const_spectrum(1.0), 1.0, dir2Light, ray_diff, k-1, 1, TransportDirection::TO_LIGHT, rng, scene, time);
    
    if(lightPath.size()<2) return;
    lightPath[1].hit_ghost_times = 0;
    const PathVertex &vertex1 = lightPath[1];
    sensorPath.push_back(vertex1);
    const Vector3 &dir2Sensor = -dir2Light;
    randomWalk_easy(sensorPath, make_const_spectrum(1.0), 1.0, dir2Sensor, ray_diff, k-1, 1, TransportDirection::TO_VIEW, rng, scene, time);

    if(sensorPath.size()<2) return;
    if(lightPath.size()>1) lightPath[1].transport.pdfFwd = 1.0;
    if(sensorPath.size()>1) sensorPath[1].transport.pdfFwd = 1.0;
    lightPath[0].transport.pdfFwd = lightPath[0].transport.pdfRev = 1.0;
    sensorPath[0].transport.pdfFwd = sensorPath[0].transport.pdfRev = 1.0;
}

