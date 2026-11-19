#pragma once

#include "scene.h"
#include "pcg.h"
#include "../utils.h"
#include "portal_utils.h"
#include "portal_utils_mix.h"
#include "incremental_basepath.h"

#define BETTER_MIS

void render_dynamic_path_1_from_light(std::vector<PathVertex> &lightPath,
                                      std::vector<PathVertex> &sensorPath,
                                      pcg32_state &rng,
                                      const Scene &scene,
                                      const int &numSamples,
                                      int time,
                                      const Spectrum &Llight,
                                      Image3 &img)
{
    assert(lightPath.size() == 2);
    const PathVertex &vDynamicS = lightPath[0], &vDynamicT = sensorPath[0];
    int nSensor = sensorPath.size();
    const PathVertex &vLight = lightPath[0]; // s-1 here s==1
    // prepare the accumulated pdfFwd and pdfRevs
    DSubPathState pathState(nSensor);
    cumprod_pdfs(sensorPath, pathState);

    for (int t = 2; t <= nSensor; ++t)
    {
        const PathVertex &vSensor = sensorPath[t - 1];
        Vector3 dirToSensor;
        Vector2 pRaster2{0.0, 0.0};
        PathVertex sampledVSensor, sampledVLight;
        Spectrum We = connectToSensorVertex(vSensor, sampledVSensor, dirToSensor, pRaster2, scene, time, true);
        if (average(We) == 0.0)
            continue;
        Spectrum Wsensor = We;
        if (t > 1)
        {
            const PathVertex &prevSensor = sensorPath[t - 2];
            Vector3 dirLight = normalize(prevSensor.position - vSensor.position); // also just dir
            Spectrum f = eval(scene.materials[vSensor.material_id], dirLight, dirToSensor, vSensor, scene.texture_pool);
            f /= fabs(dot(dirToSensor, vSensor.geometric_normal)); // TODO remove the costheta from the bsdf_eval
            Wsensor *= vSensor.transport.contrib * f;              // if s==1 and t==1 at the same time , then f is already included in the directlighting
        }

        Spectrum radiance = Llight * Wsensor * vDynamicS.transport.contrib / (Real)numSamples;
        {
            const PathVertex &vDynamicSnext = lightPath[1];
            const PathVertex &vDynamicTnext = t >= 2 ? sensorPath[1] : sampledVSensor;
            const Material &matDynamic = scene.materials[get_material_id(scene.shapes[vDynamicS.shape_id])];
            const Vector3 dir2Light = normalize(vDynamicSnext.position - vDynamicS.position);
            const Vector3 dir2Sensor = normalize(vDynamicTnext.position - vDynamicS.position);
            Spectrum fDynamic2Light = eval(matDynamic, dir2Sensor, dir2Light, vDynamicS, scene.texture_pool) / fabs(dot(dir2Light, vDynamicS.geometric_normal));
            radiance *= fDynamic2Light;
            if (t > 1)
            {
                radiance *= fabs(dot(dir2Sensor, vDynamicS.geometric_normal)); // if directly connecting to camera, no cosine on the sensor side
                radiance /= vDynamicT.transport.pdfDir;                        // this is conditional on S.pdfDir;
            }
        }
#if defined(BETTER_MIS)
        Real base_pdf;
        Real weight = computeMixMISWeight_better_dynamic(vLight, vSensor, lightPath[1], sampledVSensor,
                                                    lightPath, sensorPath, 1, t, NULL,
                                                    pathState, SamplingTechnique::DYNAMIC_FROM_LIGHT, scene, time, base_pdf);
#elif defined(NAIVE_MIS_MIX)
        Real weight = compute_mis_solid_part(lightPath[1], vLight, vSensor, sampledVSensor, SamplingTechnique::DYNAMIC_FROM_LIGHT);
#else 
        Real weight = computePortalMISWeight(vLight, vSensor); 
#endif
        radiance *= weight;
        assert(average(radiance) == average(radiance));
        assert(!isinf(radiance.x) && !isinf(radiance.y) && !isinf(radiance.z));
        splatSample(img, pRaster2, radiance);
    }
} 

void render_dynamic_path_2_from_sensor(std::vector<PathVertex> &lightPath,
                                       std::vector<PathVertex> &sensorPath,
                                       pcg32_state &rng,
                                       const Scene &scene,
                                       const int &numSamples,
                                       int time,
                                       const Spectrum &Wsensor,
                                       const Vector2 &pRaster2,
                                       Image3 &img)
{
    assert(sensorPath.size() == 2);
    const PathVertex &vDynamicS = lightPath[0], &vDynamicT = sensorPath[0];
    const PathVertex &vSensor = sensorPath[0]; // t-1 here t==1

    PathVertex lastHitLight;
    if(lightPath.back().type == VertexType::Light){
        lastHitLight = lightPath.back();
        //then remove it from the list 
        lightPath.pop_back();
        assert(lastHitLight.type == VertexType::Light);
    }
    // prepare the accumulated pdfFwd and pdfRevs
    int nLight = lightPath.size();
    DSubPathState pathState(nLight);
    cumprod_pdfs(lightPath, pathState);

    for (int s = 1; s <= nLight; ++s)
    {   
        const PathVertex &vLight = lightPath[s - 1];
        assert(vLight.material_id >= 0);
        PathVertex sampledVLight;
        Spectrum Llight = make_zero_spectrum();
                
        {
            Vector3 dir_light;
            Llight = directIllumination(vLight, dir_light, sampledVLight, rng, scene, time, true);
            if (average(Llight) == 0.0)
                continue;
            if (s > 1) //its gauranteed.
            {
                const PathVertex &prevLight = lightPath[s - 2];
                const Material &mat = scene.materials[vLight.material_id];
                Vector3 dir_view = normalize(prevLight.position - vLight.position);
                Spectrum f = eval(mat, dir_view, dir_light, vLight, scene.texture_pool) / fabs(dot(dir_light, vLight.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
                Llight *= vLight.transport.contrib * f;
            }
        }
        
        Spectrum radiance = Llight * Wsensor * vDynamicS.transport.contrib / (Real)numSamples;
        {
            const PathVertex &vDynamicSnext = s >= 2 ? lightPath[1] : sampledVLight;
            const PathVertex &vDynamicTnext = sensorPath[1];
            const Material &matDynamic = scene.materials[get_material_id(scene.shapes[vDynamicS.shape_id])];
            const Vector3 dir2Light = normalize(vDynamicSnext.position - vDynamicS.position);
            const Vector3 dir2Sensor = normalize(vDynamicTnext.position - vDynamicT.position);
            Spectrum fDynamic2Light = eval(matDynamic, dir2Sensor, dir2Light, vDynamicS, scene.texture_pool) / fabs(dot(dir2Light, vDynamicS.geometric_normal));
            radiance *= fDynamic2Light;
            if (s > 1){
                radiance *= fabs(dot(dir2Light, vDynamicS.geometric_normal));
                radiance /= vDynamicS.transport.pdfDir;
            }
            
        }
#if defined(BETTER_MIS)
        Real base_pdf;
        Real weight = computeMixMISWeight_better_dynamic(vLight, vSensor, sampledVLight, sensorPath[1],
                                                    lightPath, sensorPath, s, 1, pathState, NULL,
                                                    SamplingTechnique::DYNAMIC_FROM_SENSOR, scene, time, base_pdf);
#elif defined(NAIVE_MIS_MIX)
        Real weight = compute_mis_solid_part(sampledVLight, vLight, vSensor, sensorPath[1], SamplingTechnique::DYNAMIC_FROM_SENSOR);
#else 
        Real weight = computePortalMISWeight(vLight, vSensor);
#endif        
        radiance *= weight;
        splatSample(img, pRaster2, radiance);
    }
}

void render_dynamic_path_3_two_ends(std::vector<PathVertex> &lightPath,
                                    std::vector<PathVertex> &sensorPath,
                                    pcg32_state &rng,
                                    const Scene &scene,
                                    const int &numSamples,
                                    int time,
                                    Image3 &img,
                                    bool negative = false)
{
    if (lightPath.size() < 2 || sensorPath.size() < 2)
        return;

    PathVertex lastHitLight;
    if(lightPath.back().type == VertexType::Light){
        lastHitLight = lightPath.back();
        //then remove it from the list 
        lightPath.pop_back();
        assert(lastHitLight.type == VertexType::Light);
    }
    // /////////////////////////////////
    // // 0 stands for directly hit the dynamic instead of sampling it explicitly,
    // // which is not possible to hit the exactly same point on the surface, pdf=0
    int nLight = lightPath.size(), nSensor = sensorPath.size();
    const PathVertex &vDynamicS = lightPath[0], &vDynamicT = sensorPath[0];

    // prepare the accumulated pdfFwd and pdfRevs for both subpaths
    DSubPathState lightPathState(nLight), sensorPathState(nSensor);
    cumprod_pdfs(lightPath, lightPathState);
    cumprod_pdfs(sensorPath, sensorPathState);

    // the bsdf in the middle
    const PathVertex &vDynamicSnext = lightPath[1];
    const PathVertex &vDynamicTnext = sensorPath[1];
    const Material &matDynamic = scene.materials[get_material_id(scene.shapes[vDynamicS.shape_id])];
    const Vector3 dynamic_dir2Light = normalize(vDynamicSnext.position - vDynamicS.position);
    const Vector3 dynamic_dir2Sensor = normalize(vDynamicTnext.position - vDynamicT.position);
    Spectrum fDynamic2Light = eval(matDynamic, dynamic_dir2Sensor, dynamic_dir2Light, vDynamicS, scene.texture_pool) / fabs(dot(dynamic_dir2Light, vDynamicS.geometric_normal));

    for (int s = 2; s <= nLight; ++s) /// num of vertices now;
    {
        for (int t = 2; t <= nSensor; ++t)
        {
            // general case handeled special cases inside.
            int depth = s + t;
            if (depth < 1 || depth > scene.options.max_depth)
                continue;
            const PathVertex &vLight = lightPath[s - 1];
            const PathVertex &vSensor = sensorPath[t - 1];
            // connect to light source on the one end
            
            assert(vLight.material_id >= 0);
            PathVertex sampledVLight;
            Spectrum Llight = make_zero_spectrum();

            {
                Vector3 dir_light;    
                Llight = directIllumination(vLight, dir_light, sampledVLight, rng, scene, time, true);
                if (average(Llight) == 0.0)
                    continue;
                
                const PathVertex &prevLight = lightPath[s - 2];
                const Material &mat = scene.materials[vLight.material_id];
                Vector3 dir_view = normalize(prevLight.position - vLight.position);
                Spectrum f = eval(mat, dir_view, dir_light, vLight, scene.texture_pool)/ fabs(dot(dir_light, vLight.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
                Llight *= vLight.transport.contrib * f;

            }

            Vector3 dirToSensor;
            Vector2 pRaster2{0.0, 0.0};
            PathVertex sampledVSensor;
            Spectrum We = connectToSensorVertex(vSensor, sampledVSensor, dirToSensor, pRaster2, scene, time, true);
            if (average(We) == 0.0)
                continue;
            Spectrum Wsensor = We;
            { 
                const PathVertex &prevSensor = sensorPath[t - 2];
                Vector3 dirLight = normalize(prevSensor.position - vSensor.position); // also just dir
                Spectrum f = eval(scene.materials[vSensor.material_id], dirLight, dirToSensor, vSensor, scene.texture_pool);
                f /= fabs(dot(dirToSensor, vSensor.geometric_normal)); // TODO remove the costheta from the bsdf_eval
                Wsensor *= vSensor.transport.contrib * f;              // if s==1 and t==1 at the same time , then f is already included in the directlighting
            }

            Spectrum radiance = Llight * Wsensor * vDynamicS.transport.contrib / (Real)numSamples;

            radiance *= fDynamic2Light;
            // s>1
            radiance *= fabs(dot(dynamic_dir2Light, vDynamicS.geometric_normal));
            radiance /= vDynamicS.transport.pdfDir;
            // t>1
            radiance *= fabs(dot(dynamic_dir2Sensor, vDynamicT.geometric_normal));
            radiance /= vDynamicT.transport.pdfDir;

#ifdef BETTER_MIS
            Real base_pdf;
            Real weight = computeMixMISWeight_better_dynamic(vLight, vSensor, sampledVLight, sampledVSensor,
                                                        lightPath, sensorPath, s, t,
                                                        lightPathState, sensorPathState, SamplingTechnique::DYNAMIC_TWO_ENDS, scene, time, base_pdf);
#elif defined(NAIVE_MIS_MIX)
        Real weight = compute_mis_solid_part(sampledVLight, vLight, vSensor, sampledVSensor, SamplingTechnique::DYNAMIC_TWO_ENDS);
#else 
            Real weight = computePortalMISWeight(vLight, vSensor); // TODO the mis weight
#endif
            radiance *= weight;
            if(negative) splatSample(img, pRaster2, -radiance);
            else splatSample(img, pRaster2, radiance);
        }
    }
}

void render_dynamic_path_3_two_ends_noredundant(std::vector<PathVertex> &lightPath,
                                    std::vector<PathVertex> &sensorPath,
                                    pcg32_state &rng,
                                    const Scene &scene,
                                    const int &numSamples,
                                    int baseTimeStamp,
                                    Image3 &img)
{
    if (lightPath.size() < 2 || sensorPath.size() < 2)
        return;

    // /////////////////////////////////
    // // 0 stands for directly hit the dynamic instead of sampling it explicitly,
    // // which is not possible to hit the exactly same point on the surface, pdf=0
    int nLight = lightPath.size(), nSensor = sensorPath.size();
    const PathVertex &vDynamicS = lightPath[0], &vDynamicT = sensorPath[0];

    // prepare the accumulated pdfFwd and pdfRevs for both subpaths
    DSubPathState lightPathState(nLight), sensorPathState(nSensor);
    cumprod_pdfs(lightPath, lightPathState);
    cumprod_pdfs(sensorPath, sensorPathState);

    // the bsdf in the middle
    const PathVertex &vDynamicSnext = lightPath[1];
    const PathVertex &vDynamicTnext = sensorPath[1];
    const Material &matDynamic = scene.materials[get_material_id(scene.shapes[vDynamicS.shape_id])];
    const Vector3 dynamic_dir2Light = normalize(vDynamicSnext.position - vDynamicS.position);
    const Vector3 dynamic_dir2Sensor = normalize(vDynamicTnext.position - vDynamicT.position);
    Spectrum fDynamic2Light = eval(matDynamic, dynamic_dir2Sensor, dynamic_dir2Light, vDynamicS, scene.texture_pool) / fabs(dot(dynamic_dir2Light, vDynamicS.geometric_normal));

    std::vector<Spectrum> base_Llight_cache;
    std::vector<Spectrum> base_Wsensor_cache;
    std::vector<PathVertex> base_sampledVLight_cache;
    std::vector<PathVertex> base_sampledVSensor_cache;
    std::vector<Vector2> base_pRaster2_cache;
    Spectrum base_f_dynamic = eval(matDynamic, dynamic_dir2Sensor, dynamic_dir2Light, vDynamicS, scene.texture_pool) / fabs(dot(dynamic_dir2Light, vDynamicS.geometric_normal));

    for (int s = 2; s <= nLight; ++s)
    {
        // 1) render the base light subpath
        const PathVertex &vLight = lightPath[s - 1];
        const PathVertex &prevLight = lightPath[s - 2];

        assert(vLight.material_id >= 0);
        Vector3 dir_light;
        PathVertex sampledVLight;
        Spectrum Llight = directIllumination(vLight, dir_light, sampledVLight, rng, scene, baseTimeStamp, true);
        {//f next to the light source
            const Material &mat = scene.materials[vLight.material_id];
            Vector3 dir_view = normalize(prevLight.position - vLight.position);
            Spectrum f = eval(mat, dir_view, dir_light, vLight, scene.texture_pool) / fabs(dot(dir_light, vLight.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
            Llight *= vLight.transport.contrib * f;
        }
        base_Llight_cache.push_back(Llight); // even if this is zero spectrum
        base_sampledVLight_cache.push_back(sampledVLight);

    }

    for (int t = 2; t <= nSensor; ++t)
    {
        // 1) render the base sensor subpath
        const PathVertex &vSensor = sensorPath[t - 1];
        const PathVertex &prevSensor = sensorPath[t - 2];

        Vector3 dirToSensor;
        Vector2 base_pRaster2{0.0, 0.0};
        PathVertex base_sampledVSensor;
        Spectrum We = connectToSensorVertex(vSensor, base_sampledVSensor, dirToSensor, base_pRaster2, scene, baseTimeStamp, true);
        
        Spectrum Wsensor = We;
        { 
            Spectrum f = eval(scene.materials[vSensor.material_id], normalize(prevSensor.position - vSensor.position), dirToSensor, vSensor, scene.texture_pool) / fabs(dot(dirToSensor, vSensor.geometric_normal));
            Wsensor *= vSensor.transport.contrib * f; // if s==1 and t==1 at the same time , then f is already included in the directlighting
        }

        base_Wsensor_cache.push_back(Wsensor);
        base_pRaster2_cache.push_back(base_pRaster2);
        base_sampledVSensor_cache.push_back(base_sampledVSensor);

    }

    for (int s = 2; s <= nLight; ++s) /// num of vertices now;
    {
        for (int t = 2; t <= nSensor; ++t)
        {
            int depth = s + t;
            if (depth < 1 || depth > scene.options.max_depth)
                continue;
            // render the base Path
            const PathVertex &vLight = lightPath[s - 1];
            const PathVertex &vSensor = sensorPath[t - 1];
            int s_cache = s - 2;
            int t_cache = t - 2;

            // render the base Path
            const Spectrum &base_Llight = base_Llight_cache[s_cache];
            const Spectrum &base_Wsensor = base_Wsensor_cache[t_cache];
            const Vector2 &pRaster2 = base_pRaster2_cache[t_cache];

            if (average(base_Llight) == 0.0 || average(base_Wsensor) == 0.0)
                continue;

            Spectrum baseRadiance = base_Llight * base_Wsensor * vDynamicS.transport.contrib / (Real)numSamples;
            baseRadiance *= base_f_dynamic;
            // s>1
            baseRadiance *= fabs(dot(dynamic_dir2Light, vDynamicS.geometric_normal));
            baseRadiance /= vDynamicS.transport.pdfDir;
            // t>1
            baseRadiance *= fabs(dot(dynamic_dir2Sensor, vDynamicT.geometric_normal));
            baseRadiance /= vDynamicT.transport.pdfDir;
            const PathVertex& sampledVLight = base_sampledVLight_cache[s_cache];
            const PathVertex& sampledVSensor = base_sampledVSensor_cache[t_cache];
            Real base_pdf;
            Real base_weight = computeMixMISWeight_better_dynamic(vLight, vSensor, sampledVLight, sampledVSensor,
                                                                    lightPath, sensorPath, s, t,
                                                                    lightPathState, sensorPathState, 
                                                                    SamplingTechnique::DYNAMIC_TWO_ENDS, scene, baseTimeStamp, base_pdf);

            splatSample(img, pRaster2, baseRadiance);
        }
    }
}

/**
 * @brief sample a point on a dynamic object first;
 * then do unidirectional path tracing in both directions and nee/camera projection
 * 3 sampling strategies: (only need to cover s==1&t==1 case once)
 * 1) connect vDynamic to light first, then do unidirectional PT with NEE towards sensor, s==2, t==1..k
 * 2) connect vDynamic to sensor first, then do unidirectional PT with NEE towards light, t==1, s==2..k
 * 3) do unidirectional PT with NEE in both ends, s==2..k, t==2..k
 * @param scene
 * @param img
 * @param spps
 * @param rng
 * @param numSamples
 */
