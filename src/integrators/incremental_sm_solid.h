
#pragma once

#include "scene.h"
#include "pcg.h"
#include "parallel.h"
#include "progress_reporter.h"
#include "transform.h"
#include "utils.h"
#include "incremental_basepath.h"
#include "incremental_shiftmapping_utils.h"

/////////some utility functions////////////////
void sampleDynamicDir(PathVertex &vertex,
                   Vector3 &diro,
                   pcg32_state &rng,
                   const Scene &scene,
                   int time = 0)
{

    //////////////////
    // sample the outgoing direction
    const Dynamic &dyn_agent = scene.dynamics[vertex.dynamic_id];
    int dynamic_shape_id = dyn_agent.shape_id;
    const Shape &dynamic = (time == NEW_TIME) ? dyn_agent.base_shape : dyn_agent.old_shape;

    Real pdfDynamicDir;
    Vector2 dynamic_dir_rdn{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Vector3 dirFromDynamic = sample_dir_from_shape(dynamic, vertex, dynamic_dir_rdn, pdfDynamicDir, time); // sample the hemisphere

    // output
    diro = dirFromDynamic;
    vertex.transport.pdfDir = pdfDynamicDir;
}

void sampleDynamicPos(PathVertex &vertex,
                   pcg32_state &rng,
                   const Scene &scene,
                   int time = 0)
{
    // sample a dynamic object first using inversion according to the area of the objects
    // sample a 4d, a point on the dynamic surface first and then sample a 2d direction
    Real dynamic_w = next_pcg32_real<Real>(rng);
    Vector2 dynamic_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real shape_w = next_pcg32_real<Real>(rng);

    int dynamic_id = sample_dynamic_shape(scene, dynamic_w); // TODO temporarily set it to 0 //
    const Dynamic &dyn_agent = scene.dynamics[dynamic_id];
    int dynamic_shape_id = dyn_agent.shape_id;
    const Shape &dynamic = (time == NEW_TIME) ? dyn_agent.base_shape : dyn_agent.old_shape;
    PointAndNormal point_on_dynamic =
        sample_point_on_shape(dynamic, Vector3{}, dynamic_uv, shape_w, time);
    Real pdfDynamic = dynamic_pmf(scene, dynamic_id) *
                      pdf_point_on_shape(dynamic, point_on_dynamic, Vector3{}, time);

    vertex.type = VertexType::Surface;
    vertex.time = time;
    vertex.dynamic_id = dynamic_id;
    vertex.shape_id = dynamic_shape_id;
    vertex.material_id = get_material_id(scene.shapes[vertex.shape_id]);
    vertex.position = point_on_dynamic.position;
    vertex.geometric_normal = point_on_dynamic.normal;
    //////////////////////////////////////check below////////
    vertex.primitive_id = point_on_dynamic.primitive_id;
    vertex.st = point_on_dynamic.st;
    ShadingInfo shadingInfo = compute_shading_info(dynamic, vertex);
    vertex.shading_frame = shadingInfo.shading_frame;
    vertex.uv = shadingInfo.uv;

    vertex.hit_dynamic_times = 1;
    vertex.transport.pdfOriginA = pdfDynamic;
    vertex.transport.pdfFwd = 1.0; // TODO change this; this should be the pdf to generate this dynamic vertex along random walk direction
    vertex.transport.pathLength = 1;
    vertex.transport.contrib = make_const_spectrum(1.0 / pdfDynamic); // FIXME attention
}
///////////////////////////////////////////////////
/**
 * @brief this file is used to render the difference between 
 *  1) the paths interacted with the solid dynamic object and 
 *  2) its corresponding path created using shift mapping
 *  Hybrid shift mapping strategy: using the same random seed replay until hitting a diffuse surface. 
 *  
 *  step1: We generate the base path
 *  step2: we generate the shfited path using shift mapping
 *  step3: we calculate the radiance and calculate the difference radiance of the pair.  
 */
void compute_shifted_dynamic_path_1_from_light(const PathVertex &vDynamicS,
                                               const PathVertex &shiftedVDynamic,
                                               const std::vector<PathVertex> &lightPath,
                                               const std::vector<PathVertex> &sensorPath,
                                               pcg32_state &base_rng,
                                               pcg32_state &shifted_rng,  //copy of the previous base rng after sampling the dynamic position
                                               const Scene &scene,
                                               const int &numSamples,
                                               int baseTimeStamp,
                                               int shiftedTimeStamp,
                                               const Spectrum &base_Llight,
                                               Image3 &img)
{
    assert(lightPath.size() == 2);
    const PathVertex &base_sampledVLight = lightPath[1];
    int nSensor = sensorPath.size();
    const PathVertex &vLight = lightPath[0]; // s-1 here s==1
    /////////////////////step2: compose the shifted path//////////////////////////
    std::vector<PathVertex> shifted_lightPath, shifted_sensorPath;
    std::vector<SHIFT_MAP_STATUS> shifted_status;
    shifted_lightPath.push_back(shiftedVDynamic);
    shifted_sensorPath.push_back(shiftedVDynamic);
    shifted_status.push_back(SHIFT_MAP_STATUS::DIVERGED);
    
    SHIFT_MAP_STATUS status = SHIFT_MAP_STATUS::DIVERGED;
    PathVertex shifted_sampledVLight;
    Vector3 dynamic_dir2Light;
    Spectrum shifted_Llight = directIllumination(shiftedVDynamic, dynamic_dir2Light, shifted_sampledVLight, shifted_rng, scene, shiftedTimeStamp, true); // Le/pdfLight * G_to_ghost??  theG_to_ghost is not right inside

    if (average(shifted_Llight) == 0.0)
        status = SHIFT_MAP_STATUS::FAILED;
    else
        shifted_lightPath.push_back(shifted_sampledVLight);

    for(int t = 2; t <= nSensor;++t){
        if (status == SHIFT_MAP_STATUS::FAILED) break;
        const PathVertex &vSensor = sensorPath[t - 1]; //the base path; to find the (t-1)th vertex of the shifted path
        const PathVertex &prevSensor = sensorPath[t - 2];
        PathVertex &shifted_prevSensor = shifted_sensorPath.back();//this is also t-2
        
        if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED || status == SHIFT_MAP_STATUS::CONNECTED){ //the main difference is when calculating the radiance, not here
            if(is_vertex_dynamic(vSensor, scene)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }            
             // test whether the new segment is blocked by the dynamic
            
            Ray shadow_ray{vSensor.position, normalize(prevSensor.position - vSensor.position), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(vSensor.position, prevSensor.position), shiftedTimeStamp};
            if (occluded_dynamic(scene, shadow_ray))
            {
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            //otherwise, udpate the ghost intersection, the pdfFwd at shifted_vSensor is actually affected
            PathVertex shifted_vSensor(vSensor); 
            shifted_vSensor.time = shiftedTimeStamp;
            assert(t>=3);
            updatePdfsOfSegment(shifted_sensorPath[t-3], shifted_prevSensor, shifted_vSensor, t<nSensor? sensorPath[t]:PathVertex{}, scene, shiftedTimeStamp);
            updateGhostIsectInfo(shifted_vSensor, shifted_prevSensor, scene, shiftedTimeStamp); 
            shifted_sensorPath.push_back(shifted_vSensor);
            if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED)
                status = SHIFT_MAP_STATUS::CONNECTED;
            shifted_status.push_back(status);

        }else if(connectability_condition(vSensor, prevSensor, shifted_prevSensor, scene)){  //getVertexRoughnessType(vSensor, scene) == MATType::VERTEX_TYPE_DIFFUSE){
            
            assert(status == SHIFT_MAP_STATUS::DIVERGED);
            if(is_vertex_dynamic(vSensor, scene)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }

            //do the path reconnection here from shifted_prevSensor to vSensor
            bool isConnected = connectVertices(shifted_prevSensor, vSensor, scene, shiftedTimeStamp);
            if(!isConnected){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            //if the reconnection succeed,
            // 1) need to retest the ghost intersections for shifted_vSensor and shifted_sampledVSensor
            // 2) need to recompute the pdfFwd and pdfRev //TODO check what others are changed??
            PathVertex shifted_vSensor(vSensor); 
            shifted_vSensor.time = shiftedTimeStamp; 
            updatePdfsOfSegment(t==2? shifted_sampledVLight: shifted_sensorPath[t-3], shifted_prevSensor, shifted_vSensor, t<nSensor? sensorPath[t]:PathVertex{}, scene, shiftedTimeStamp);
            updateGhostIsectInfo(shifted_vSensor, shifted_prevSensor, scene, shiftedTimeStamp); 
            shifted_sensorPath.push_back(shifted_vSensor);
            status = SHIFT_MAP_STATUS::RECENTLY_CONNECTED;
            shifted_status.push_back(status);

        }else{ //otherwise random walk using the same random seed

            assert(status == SHIFT_MAP_STATUS::DIVERGED);
            Vector3 shifted_nextDir;
            Real pdfNextDir;
            bool isOtherDir = sampleScattering(shifted_prevSensor, normalize((t==2? shifted_sampledVLight.position: shifted_sensorPath[t-3].position) - shifted_prevSensor.position), 
                    shifted_nextDir, pdfNextDir, shifted_rng, scene);
            if(!isOtherDir){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            if(t==2) shifted_sensorPath[t - 2].transport.pdfDir = pdfNextDir;///filling in the pdfDir used later.
            //just one step forward
            RayDifferential ray_diff;
            randomWalk(shifted_sensorPath, 
                t==2? make_const_spectrum(1.0) : shifted_prevSensor.transport.contrib, 
                pdfNextDir, shifted_nextDir, ray_diff, t+1, t-1, 
                TransportDirection::TO_VIEW, shifted_rng, scene, shiftedTimeStamp);

            if(shifted_sensorPath.size() != t){ //random walk fails
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            shifted_status.push_back(status);
        }
        assert(shifted_sensorPath.size() == shifted_status.size());
    }    

    /////////////////////step3: calculate the radiance and calculate the difference radiance of the pair //////////////////////////
    int nSensor_shifted = shifted_sensorPath.size();
    //prepare the accumulated pdfFwd and pdfRevs used for MIS
    DSubPathState pathState(nSensor);
    cumprod_pdfs(sensorPath, pathState);
    DSubPathState shifted_pathState(nSensor_shifted);
    cumprod_pdfs(shifted_sensorPath, shifted_pathState);

    Spectrum ratio_T = make_const_spectrum(1.0);
    for (int t = 2; t <= nSensor; ++t) //put s==1 and t==1 into DYNAMIC_FROM_SENSOR
    {
        // 1) render the base path
        const PathVertex &vSensor = sensorPath[t - 1];
        Vector3 dirToSensor;
        Vector2 base_pRaster2{0.0, 0.0};
        PathVertex base_sampledVSensor;
        Spectrum base_We = connectToSensorVertex(vSensor, base_sampledVSensor, dirToSensor, base_pRaster2, scene, baseTimeStamp, true);
        Spectrum Wsensor = base_We;
        //record part of the ratiance for ratio calculation later
        Spectrum cache_for_ratio = base_Llight;
        if (t > 1)
        { //f at the last vertex next to sensor
            const PathVertex &prevSensor = sensorPath[t - 2];
            Vector3 dirLight = normalize(prevSensor.position - vSensor.position); // also just dir
            Spectrum f = eval(scene.materials[vSensor.material_id], dirLight, dirToSensor, vSensor, scene.texture_pool) / fabs(dot(dirToSensor, vSensor.geometric_normal)); 
            Wsensor *= f;              // if s==1 and t==1 at the same time , then f is already included in the directlighting
            cache_for_ratio *= vSensor.transport.contrib; 
        }

        Spectrum base_fDynamic2Light;
        { //this part is repeated
            const PathVertex &vDynamicSnext = lightPath[1];
            const PathVertex &vDynamicTnext = t >= 2 ? sensorPath[1] : base_sampledVSensor;
            const Material &matDynamic = scene.materials[get_material_id(scene.shapes[vDynamicS.shape_id])];
            const Vector3 dir2Light = normalize(vDynamicSnext.position - vDynamicS.position);
            const Vector3 dir2Sensor = normalize(vDynamicTnext.position - vDynamicS.position);
            base_fDynamic2Light = eval(matDynamic, dir2Sensor, dir2Light, vDynamicS, scene.texture_pool) / fabs(dot(dir2Light, vDynamicS.geometric_normal));
            cache_for_ratio *= base_fDynamic2Light;
            if (t > 1)
            {
                cache_for_ratio *= fabs(dot(dir2Sensor, vDynamicS.geometric_normal)); // if directly connecting to camera, no cosine on the sensor side
                cache_for_ratio /= sensorPath[0].transport.pdfDir;                        // this is conditional on S.pdfDir;
            }
        }

        Spectrum baseRadiance = cache_for_ratio * Wsensor * vDynamicS.transport.contrib / (Real)numSamples;
        Spectrum unweighted_baseRadiance = baseRadiance;
        Real base_pdf;
        Real base_weight = computeMixMISWeight_better_dynamic(vLight, vSensor, base_sampledVLight, base_sampledVSensor,
                                                    lightPath, sensorPath, 1, t, NULL,
                                                    pathState, SamplingTechnique::DYNAMIC_FROM_LIGHT, scene, baseTimeStamp, base_pdf);
        baseRadiance *= base_weight;
        if(base_pdf==0.0)continue;
        if(average(baseRadiance)==0.0)
        {
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
            continue;
        }

        // 2) render the shifted path //////////////////////////////////////////
        if(t > nSensor_shifted){ //no more shifted vertex; just calculate the base path and continue
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
            continue;
        }

        const PathVertex &shifted_vSensor = shifted_sensorPath[t - 1];
        const SHIFT_MAP_STATUS &status = shifted_status[t-1];
        
        if(status == SHIFT_MAP_STATUS::DIVERGED){
            //record the ratio until the last divergence
            ratio_T = shifted_vSensor.transport.contrib / vSensor.transport.contrib;

            //connect to the sampled sensor
            PathVertex shifted_sampledVSensor;
            Vector3 dirToSensor;
            Vector2 shifted_pRaster2{0.0, 0.0};
            Spectrum shifted_We = connectToSensorVertex(shifted_vSensor, shifted_sampledVSensor, dirToSensor, shifted_pRaster2, scene, shiftedTimeStamp, true);
            if (average(shifted_We) == 0.0)//SM failed
            {
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                continue; // doesn't affect the status for the whole path
            }

            if(t==1){
                Spectrum shiftedRadiance = shifted_Llight * shifted_We * vDynamicS.transport.contrib / (Real)numSamples;
                Vector3 shifted_wi = normalize(shifted_sampledVLight.position - shiftedVDynamic.position);
                Vector3 wo = normalize(shifted_sampledVSensor.position - shiftedVDynamic.position);
                const Material &matVsensor = scene.materials[get_material_id(scene.shapes[shiftedVDynamic.shape_id])];
                Spectrum f_shifted = eval(matVsensor, shifted_wi, wo, shiftedVDynamic, scene.texture_pool) / fabs(dot(wo, shiftedVDynamic.geometric_normal));
                shiftedRadiance *= f_shifted;
                Real shifted_pdf;
                Real shifted_weight = computeMixMISWeight_better_dynamic(vLight, shifted_vSensor, shifted_sampledVLight, shifted_sampledVSensor,
                                                        shifted_lightPath, shifted_sensorPath, 1, t, NULL,
                                                        shifted_pathState, SamplingTechnique::DYNAMIC_FROM_LIGHT, scene, shiftedTimeStamp, shifted_pdf);
                shiftedRadiance *= shifted_weight;
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
                splatSample(img, shifted_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));
                continue;
            }

            Spectrum shifted_WSensor(shifted_We);

            {//f next to the sensor
                const PathVertex &shifted_prevSensor = shifted_sensorPath[t - 2];
                Vector3 shifted_dirLight = normalize(shifted_prevSensor.position - shifted_vSensor.position); // also just dir
                Vector3 shifted_dirToSensor = normalize(shifted_sampledVSensor.position - shifted_vSensor.position);
                Spectrum f = eval(scene.materials[shifted_vSensor.material_id], shifted_dirLight, shifted_dirToSensor, shifted_vSensor, scene.texture_pool) / fabs(dot(shifted_dirToSensor, shifted_vSensor.geometric_normal)); 
                shifted_WSensor *= f;              // if s==1 and t==1 at the same time , then f is already included in the directlighting                
            }
            Spectrum shiftedRadiance = shifted_Llight * shifted_WSensor * vDynamicS.transport.contrib / (Real)numSamples;

            shiftedRadiance *= shifted_vSensor.transport.contrib;
            {//f at dynamic //TODO this is reduntantly calculated
                Vector3 shifted_wi = normalize(shifted_sampledVLight.position - shiftedVDynamic.position);
                Vector3 wo = normalize(t>=2? shifted_sensorPath[1].position:shifted_sampledVSensor.position - shiftedVDynamic.position);
                const Material &matShiftedVDynamic = scene.materials[get_material_id(scene.shapes[shiftedVDynamic.shape_id])];
                Spectrum shifted_fDynamic2Light = eval(matShiftedVDynamic, shifted_wi, wo, shiftedVDynamic, scene.texture_pool) / fabs(dot(wo, shiftedVDynamic.geometric_normal));
                shiftedRadiance *= shifted_fDynamic2Light;
                shiftedRadiance *= fabs(dot(wo, shiftedVDynamic.geometric_normal)); // if directly connecting to camera, no cosine on the sensor side
                shiftedRadiance /= shifted_sensorPath[0].transport.pdfDir;                        // this is conditional on S.pdfDir;
            }
            //calculate the mis weight for the shifted path
            Real shifted_pdf;
            Real shifted_weight = computeMixMISWeight_better_dynamic(shifted_lightPath[1], shifted_vSensor, shifted_sampledVLight, shifted_sampledVSensor,
                                                        shifted_lightPath, shifted_sensorPath, 1, t, NULL,
                                                        shifted_pathState, SamplingTechnique::DYNAMIC_FROM_LIGHT, scene, shiftedTimeStamp, shifted_pdf);
            shiftedRadiance *= shifted_weight;
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
            splatSample(img, shifted_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));

        }else if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED){
            ratio_T *= shifted_Llight / division_protect(base_Llight);

            const PathVertex &prevSensor = sensorPath[t - 2];
            const PathVertex &shifted_prevSensor = shifted_sensorPath[t - 2];

            {//the change of f at dynamic
                const Material &matShiftedVDynamic = scene.materials[shiftedVDynamic.material_id];
                Vector3 shifted_wi = normalize(shifted_sampledVLight.position - shiftedVDynamic.position);
                Vector3 shifted_wo = normalize(shifted_sensorPath[1].position - shiftedVDynamic.position);
                Spectrum shifted_f_dynamic = eval(matShiftedVDynamic, shifted_wi, shifted_wo, shiftedVDynamic, scene.texture_pool) / fabs(dot(shifted_wo, shiftedVDynamic.geometric_normal));
                // change of the bsdf at the vDynamic
                ratio_T *= shifted_f_dynamic / division_protect(base_fDynamic2Light);
                Vector3 base_wo = normalize(shifted_sensorPath[1].position - vDynamicS.position);
                ratio_T *= fabs(dot(shifted_wo, shiftedVDynamic.geometric_normal) / dot(base_wo, vDynamicS.geometric_normal)); // cosine at vDynamic
            }
            { // account for the jacobian due to reconnection
                Real base_cosTheta = fabs(dot(normalize(vSensor.position - prevSensor.position), vSensor.geometric_normal));
                Real shifted_cosTheta = fabs(dot(normalize(vSensor.position - shifted_prevSensor.position), vSensor.geometric_normal));
                Real base_distance2 = distance_squared(vSensor.position, prevSensor.position);
                Real shifted_distance2 = distance_squared(vSensor.position, shifted_prevSensor.position);
                // the jacobian of the pdf
                ratio_T *= base_distance2 * shifted_cosTheta / base_cosTheta / shifted_distance2;
            }
            if(average(ratio_T)>10.)
            {
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                continue;
            }

            Spectrum shiftedRadiance = unweighted_baseRadiance * ratio_T;
            if(t < nSensor && t < nSensor_shifted){ 
                const Material &matVsensor = scene.materials[vSensor.material_id];
                Vector3 wi = normalize(prevSensor.position - vSensor.position);
                Vector3 wo = normalize(sensorPath[t].position - vSensor.position);
                Spectrum f_base = eval(matVsensor, wi, wo, vSensor, scene.texture_pool) / fabs(dot(wo, vSensor.geometric_normal));
                wi = normalize(shifted_prevSensor.position - vSensor.position);
                Spectrum f_shifted = eval(matVsensor, wi, wo, vSensor, scene.texture_pool) / fabs(dot(wo, vSensor.geometric_normal));
                ratio_T *= f_shifted / division_protect(f_base);
            }

            { // test whether the segment connected to sensor is blocked by the dynamic
                Ray shadow_ray{vSensor.position, dirToSensor, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(vSensor.position, base_sampledVSensor.position), shiftedTimeStamp};
                if (occluded_dynamic(scene, shadow_ray))
                {
                    splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                    continue;
                }
            }
            //update the ghost intersection of the base_sampledVSensor
            PathVertex shifted_sampledVSensor(base_sampledVSensor);
            shifted_sampledVSensor.time = shiftedTimeStamp;
            updateGhostIsectInfo(shifted_sampledVSensor, shifted_vSensor, scene, shiftedTimeStamp);
            Real shifted_pdf;
            Real shifted_weight = computeMixMISWeight_better_dynamic(shifted_lightPath[1], shifted_vSensor, shifted_sampledVLight, shifted_sampledVSensor,
                                                        shifted_lightPath, shifted_sensorPath, 1, t, NULL,
                                                        shifted_pathState, SamplingTechnique::DYNAMIC_FROM_LIGHT, scene, shiftedTimeStamp, shifted_pdf);
            shiftedRadiance *= shifted_weight;
            // the bsdf value at (vSensor) is changed due to the changed incoming direction
            {
                Vector3 base_wi = normalize(prevSensor.position - vSensor.position);
                Vector3 shifted_wi = normalize(shifted_prevSensor.position - vSensor.position);
                Vector3 wo = normalize(base_sampledVSensor.position - vSensor.position);
                const Material &matVsensor = scene.materials[vSensor.material_id];
                Spectrum f_base = eval(matVsensor, base_wi, wo, vSensor, scene.texture_pool) / fabs(dot(wo, vSensor.geometric_normal));
                Spectrum f_shifted = eval(matVsensor, shifted_wi, wo, vSensor, scene.texture_pool) / fabs(dot(wo, vSensor.geometric_normal));
                shiftedRadiance *= f_shifted / division_protect(f_base);
            }
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
            splatSample(img, base_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));

        }else if(status == SHIFT_MAP_STATUS::CONNECTED)
        {// use the cached ratio

            if(average(ratio_T)>10.)
            {
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                continue;
            }
            
            { // test whether the segment connected to sensor is blocked by the dynamic
                Ray shadow_ray{vSensor.position, dirToSensor, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(vSensor.position, base_sampledVSensor.position), shiftedTimeStamp};
                if (occluded_dynamic(scene, shadow_ray))
                {
                    splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                    continue;
                }
            }
            //update the ghost intersection of the base_sampledVSensor
            PathVertex shifted_sampledVSensor(base_sampledVSensor);
            shifted_sampledVSensor.time = shiftedTimeStamp;
            updateGhostIsectInfo(shifted_sampledVSensor, shifted_vSensor, scene, shiftedTimeStamp);
            Real shifted_pdf;
            Real shifted_weight = computeMixMISWeight_better_dynamic(shifted_lightPath[1], shifted_vSensor, shifted_sampledVLight, shifted_sampledVSensor,
                                                        shifted_lightPath, shifted_sensorPath, 1, t, NULL,
                                                        shifted_pathState, SamplingTechnique::DYNAMIC_FROM_LIGHT, scene, shiftedTimeStamp, shifted_pdf);
            
            Spectrum shiftedRadiance = shifted_weight * unweighted_baseRadiance * ratio_T;
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
            splatSample(img, base_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));

        }else{
            assert(false);
        }
    }
}

void compute_shifted_dynamic_path_2_from_sensor(const PathVertex &vDynamic,
                                                const PathVertex &shiftedVDynamic,
                                                const std::vector<PathVertex> &lightPath,
                                                const std::vector<PathVertex> &sensorPath,
                                                pcg32_state &rng,
                                                pcg32_state &shifted_rng,  //copy of the previous base rng after sampling the dynamic position
                                                const Scene &scene,
                                                const int &numSamples,
                                                int baseTimeStamp,
                                                int shiftedTimeStamp,
                                                const Spectrum &base_We,
                                                const Vector2 base_pRaster2,
                                                Image3 &img)
{
    assert(sensorPath.size() == 2);
    const PathVertex &vDynamicS = lightPath[0], &vDynamicT = sensorPath[0];
    int nLight = lightPath.size();
    const PathVertex &vSensor = sensorPath[0]; // t-1 here t==1; basically dummy vertex
    const PathVertex &base_sampledVSensor = sensorPath[1];
    /////////////////////step2: compose the shifted path//////////////////////////
    std::vector<PathVertex> shifted_lightPath, shifted_sensorPath;
    std::vector<SHIFT_MAP_STATUS> shifted_status;
    shifted_lightPath.push_back(shiftedVDynamic);
    shifted_sensorPath.push_back(shiftedVDynamic);
    shifted_status.push_back(SHIFT_MAP_STATUS::DIVERGED);

    SHIFT_MAP_STATUS status = SHIFT_MAP_STATUS::DIVERGED;
    PathVertex shifted_sampledVSensor;
    Vector3 dirToSensor;
    Vector2 shifted_pRaster2{0.0, 0.0};
    Spectrum shifted_We = connectToSensorVertex(shiftedVDynamic, shifted_sampledVSensor, dirToSensor, shifted_pRaster2, scene, shiftedTimeStamp, true);
    if (average(shifted_We) == 0.0)
        status = SHIFT_MAP_STATUS::FAILED;
    else 
        shifted_sensorPath.push_back(shifted_sampledVSensor);
    
    for (int s = 2; s <= nLight; ++s)
    {
        if (status == SHIFT_MAP_STATUS::FAILED) break;

        const PathVertex &vLight = lightPath[s - 1];//the base path; now let's find the (s-1)th shifted vertex
        PathVertex &shifted_prevLight = shifted_lightPath.back();
        
        if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED || status == SHIFT_MAP_STATUS::CONNECTED){ //the main difference is when calculating the radiance, not here
            if(is_vertex_dynamic(vLight, scene)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            const PathVertex &prevLight = lightPath[s - 2];
            Ray shadow_ray{vLight.position, normalize(prevLight.position - vLight.position), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(vLight.position, prevLight.position), shiftedTimeStamp};
            if (occluded_dynamic(scene, shadow_ray))
            {
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            //otherwise, update the ghost intersection information at shiftedVSensor for mis weight, the pdfs 
            PathVertex shifted_vLight(vLight); shifted_vLight.time = shiftedTimeStamp;
            assert(s>=3);
            updatePdfsOfSegment(shifted_lightPath[s-3], shifted_prevLight, shifted_vLight, s<nLight? lightPath[s]:PathVertex{}, scene, shiftedTimeStamp);
            updateGhostIsectInfo(shifted_vLight, shifted_prevLight, scene, shiftedTimeStamp);
            shifted_lightPath.push_back(shifted_vLight);
            if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED)
                status = SHIFT_MAP_STATUS::CONNECTED;
            shifted_status.push_back(status);

        }else if(getVertexRoughnessType(vLight, scene) == MATType::VERTEX_TYPE_DIFFUSE){
            
            assert(status == SHIFT_MAP_STATUS::DIVERGED);
            //do the path reconnection here from the shifted_prevLight to the vLight
            bool isConnected = connectVertices(shifted_prevLight, vLight, scene, shiftedTimeStamp);
            if (isConnected)
            {
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            //if the reconnection succeed,
            // 1) need to retest the ghost intersections for shifted_v and shifted_sampledVEnd
            // 2) need to recompute the pdfFwd and pdfRev //TODO check what others are changed??
            PathVertex shifted_vLight(vLight);
            shifted_vLight.time = shiftedTimeStamp;
            updatePdfsOfSegment(s==2?shifted_sampledVSensor:shifted_lightPath[s-3], shifted_prevLight, shifted_vLight, s<nLight? lightPath[s]:PathVertex{}, scene, shiftedTimeStamp);
            updateGhostIsectInfo(shifted_vLight, shifted_prevLight, scene, shiftedTimeStamp);
            shifted_lightPath.push_back(shifted_vLight);
            status = SHIFT_MAP_STATUS::RECENTLY_CONNECTED;
            shifted_status.push_back(status);

        }else{//if not reconnectable; just do random walk using the same random seed

            assert(status == SHIFT_MAP_STATUS::DIVERGED);
            Vector3 shifted_nextDir;
            Real pdfNextDir;
            bool isOtherDir = sampleScattering(shifted_prevLight, normalize((s==2? shifted_sampledVSensor.position:shifted_lightPath[s-3].position) - shifted_prevLight.position),
                shifted_nextDir, pdfNextDir, shifted_rng, scene); 
            if(!isOtherDir){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            if(s==2) shifted_lightPath[s-2].transport.pdfDir = pdfNextDir;
            //just one step forward
            RayDifferential ray_diff;
            randomWalk(shifted_lightPath, s==2?make_const_spectrum(1.0):shifted_prevLight.transport.contrib, 
                pdfNextDir, shifted_nextDir, ray_diff, s+1, s-1, TransportDirection::TO_LIGHT, 
                shifted_rng, scene, shiftedTimeStamp);
            {//TODO something to do with it.
                PathVertex lastHitLight;
                if(shifted_lightPath.back().type == VertexType::Light){
                    lastHitLight = shifted_lightPath.back();
                    //then remove it from the list 
                    shifted_lightPath.pop_back();
                    assert(lastHitLight.type == VertexType::Light);
                }
            
            }
            if(shifted_lightPath.size() != s){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            shifted_status.push_back(status);
        }

        assert(shifted_lightPath.size() == shifted_status.size());
    }

    /////////////////////step3: calculate the radiance and calculate the difference radiance of the pair //////////////////////////
    int nLight_shifted = shifted_lightPath.size();
    //prepare the accumulated pdfFwd and pdfRevs used for MIS
    DSubPathState pathState(nLight);
    cumprod_pdfs(lightPath, pathState);
    DSubPathState shifted_pathState(nLight_shifted);
    cumprod_pdfs(shifted_lightPath, shifted_pathState);

    Spectrum ratio_S = make_const_spectrum(1.0);
    for (int s = 1; s <= nLight; ++s)
    { // s==1 and t==1 is handled in case 1; so here start from 2
        const PathVertex &vLight = lightPath[s - 1];
        Vector3 dir_light;
        assert(vLight.material_id >= 0);
        PathVertex base_sampledVLight;
        Spectrum base_Llight = directIllumination(vLight, dir_light, base_sampledVLight, rng, scene, baseTimeStamp, true);

        if (s > 1){//f at the last vertex next to light
            const PathVertex &prevLight = lightPath[s - 2];
            const Material &mat = scene.materials[vLight.material_id];
            Vector3 dir_view = normalize(prevLight.position - vLight.position);
            Spectrum f = eval(mat, dir_view, dir_light, vLight, scene.texture_pool) / fabs(dot(dir_light, vLight.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
            base_Llight *= vLight.transport.contrib * f;
        }
        Spectrum baseRadiance = base_Llight * base_We * vDynamicS.transport.contrib / (Real)numSamples;

        Spectrum base_f_dynamic = make_const_spectrum(1.0);
        {//f at the dynamic vertex
            const PathVertex &vDynamicSnext = s >= 2 ? lightPath[1] : base_sampledVLight;
            const PathVertex &vDynamicTnext = sensorPath[1];
            const Material &matDynamic = scene.materials[vDynamicS.material_id];
            const Vector3 dir2Light = normalize(vDynamicSnext.position - vDynamicS.position);
            const Vector3 dir2Sensor = normalize(vDynamicTnext.position - vDynamicS.position);
            base_f_dynamic = eval(matDynamic, dir2Sensor, dir2Light, vDynamicS, scene.texture_pool) / fabs(dot(dir2Light, vDynamicS.geometric_normal));
            baseRadiance *= base_f_dynamic;
            if(s>1){
                baseRadiance *= fabs(dot(dir2Light, vDynamicS.geometric_normal));
                baseRadiance /= vDynamicS.transport.pdfDir;
            }

        }

        Spectrum unweighted_baseRadiance = baseRadiance;
        Real base_pdf;
        Real base_weight = computeMixMISWeight_better_dynamic(vLight, vSensor, base_sampledVLight, sensorPath[1],
                                                    lightPath, sensorPath, s, 1, pathState, NULL,
                                                    SamplingTechnique::DYNAMIC_FROM_SENSOR, scene, baseTimeStamp, base_pdf);

        baseRadiance *= base_weight;
        
        if(s > nLight_shifted){//no more shifted vertex; just calculate the base path and continue
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
            continue;
        }

        const PathVertex &shifted_vLight = shifted_lightPath[s - 1];
        const SHIFT_MAP_STATUS &status = shifted_status[s - 1];

        if(status == SHIFT_MAP_STATUS::DIVERGED){
            
            ratio_S = shifted_vLight.transport.contrib / vLight.transport.contrib;
            //connect to the sampled light source  
            PathVertex shifted_sampledVLight(base_sampledVLight);
            //TOOD connect to retest the ghost intersection and calculate the shifted_Llight
            Spectrum shifted_Llight = directIlluminationGivenLight(shifted_vLight, shifted_sampledVLight, scene, shiftedTimeStamp, true, true); //TODO
            if (average(shifted_Llight) == 0.0)//SM failed
            {
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                continue; // doesn't affect the status for the whole path
            }
            if(s == 1){
                Spectrum shiftedRadiance = shifted_Llight *shifted_We * vDynamicS.transport.contrib / (Real)numSamples;
                const Material &matShiftedVDynamic = scene.materials[shiftedVDynamic.material_id];
                const Vector3 shifted_wi = normalize(shifted_sampledVSensor.position - shiftedVDynamic.position);
                const Vector3 shifted_wo = normalize(shifted_sampledVLight.position - shiftedVDynamic.position);
                Spectrum shifted_f_dynamic = eval(matShiftedVDynamic, shifted_wi, shifted_wo, shiftedVDynamic, scene.texture_pool) / fabs(dot(shifted_wo, shiftedVDynamic.geometric_normal));
                shiftedRadiance *= shifted_f_dynamic;
                Real shifted_pdf;
                Real shifted_weight = computeMixMISWeight_better_dynamic(shifted_vLight, shifted_sensorPath[1], shifted_sampledVLight, shifted_sampledVSensor,
                                                    shifted_lightPath, shifted_sensorPath, s, 1, shifted_pathState,
                                                    NULL, SamplingTechnique::DYNAMIC_FROM_SENSOR, scene, shiftedTimeStamp, shifted_pdf);
                shiftedRadiance *= shifted_weight;
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
                splatSample(img, shifted_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));
                continue;
            }
            {//f next to light source
                const PathVertex &shifted_prevLight = shifted_lightPath[s - 2];
                const Material &mat = scene.materials[shifted_vLight.material_id];
                Vector3 shifted_dirToSensor = normalize(shifted_prevLight.position - shifted_vLight.position);
                Vector3 shifted_dirLight = normalize(base_sampledVLight.position - shifted_vLight.position);
                Spectrum f = eval(mat, shifted_dirToSensor, shifted_dirLight, shifted_vLight, scene.texture_pool) / fabs(dot(shifted_dirLight, shifted_vLight.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
                shifted_Llight *=  f;
            }

             Spectrum shiftedRadiance = shifted_Llight * shifted_We * vLight.transport.contrib * vDynamicS.transport.contrib / (Real)numSamples;
            {//TODO this is redundant
                //acount for the change for f at dynamic
                const Material &matShiftedVDynamic = scene.materials[shiftedVDynamic.material_id];
                const Vector3 shifted_wi = normalize(base_sampledVSensor.position - shiftedVDynamic.position);
                const Vector3 shifted_wo = normalize(vLight.position - shiftedVDynamic.position);
                Spectrum shifted_f_dynamic = eval(matShiftedVDynamic, shifted_wi, shifted_wo, shiftedVDynamic, scene.texture_pool) / fabs(dot(shifted_wo, shiftedVDynamic.geometric_normal));
                shiftedRadiance *= shifted_f_dynamic;
                shiftedRadiance *= fabs(dot(shifted_wo, shiftedVDynamic.geometric_normal));
                shiftedRadiance /= shifted_lightPath[0].transport.pdfDir;                        // this is conditional on S.pdfDir;
            }

            //calculate the mis weight for the shifted path
            Real shifted_pdf;
            Real shifted_weight = computeMixMISWeight_better_dynamic(shifted_vLight, shifted_sensorPath[1], shifted_sampledVLight, shifted_sampledVSensor,
                                                    shifted_lightPath, shifted_sensorPath, s, 1, shifted_pathState,
                                                    NULL, SamplingTechnique::DYNAMIC_FROM_SENSOR, scene, shiftedTimeStamp, shifted_pdf);
            shiftedRadiance *= shifted_weight;
            assert(average(baseRadiance) == average(baseRadiance));
            assert(average(shiftedRadiance) == average(shiftedRadiance));
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
            splatSample(img, shifted_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));
        
        }else if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED){

            ratio_S *= shifted_We / base_We;

            const PathVertex &prevLight = lightPath[s - 2];
            const PathVertex &shifted_prevLight = shifted_lightPath[s - 2];
            {//acount for the change for f at dynamic
                const Material &matShiftedVDynamic = scene.materials[shiftedVDynamic.material_id];
                const Vector3 base_wo = normalize(lightPath[1].position - vDynamicS.position);
                const Vector3 shifted_wi = normalize(shifted_sampledVSensor.position - shiftedVDynamic.position);
                const Vector3 shifted_wo = normalize(shifted_lightPath[1].position - shiftedVDynamic.position);
                Spectrum shifted_f_dynamic = eval(matShiftedVDynamic, shifted_wi, shifted_wo, shiftedVDynamic, scene.texture_pool) / fabs(dot(shifted_wo, shiftedVDynamic.geometric_normal));
                ratio_S *= shifted_f_dynamic / division_protect(base_f_dynamic);
                ratio_S *= fabs(dot(shifted_wo, shiftedVDynamic.geometric_normal) / dot(base_wo, vDynamicS.geometric_normal));
            }
            { // due to the path reconnection
                Real base_cosTheta = fabs(dot(normalize(prevLight.position - vLight.position), vLight.geometric_normal));
                Real shifted_cosTheta = fabs(dot(normalize(shifted_prevLight.position - vLight.position), vLight.geometric_normal));
                Real base_distance2 = distance_squared(prevLight.position, vLight.position);
                Real shifted_distance2 = distance_squared(shifted_prevLight.position, vLight.position);
                Real jacobian = shifted_cosTheta * base_distance2 / base_cosTheta / shifted_distance2;
                ratio_S *= jacobian;
            }
            if(average(ratio_S)>10.)
            {
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                continue;
            }
            
            Spectrum shiftedRadiance = unweighted_baseRadiance * ratio_S;

            // // account for the ratio of f at the reconnected vertex
            if(s < nLight && s < nLight_shifted){
                const Material &matVlight = scene.materials[vLight.material_id];
                Vector3 base_wi = normalize(prevLight.position - vLight.position);
                Vector3 wo = normalize(lightPath[s].position - vLight.position);
                Vector3 shifted_wi = normalize(shifted_vLight.position - vLight.position);
                Spectrum f_base = eval(matVlight, base_wi, wo, vLight, scene.texture_pool) / fabs(dot(wo, vLight.geometric_normal));
                Spectrum f_shifted = eval(matVlight, shifted_wi, wo, vLight, scene.texture_pool) / fabs(dot(wo, vLight.geometric_normal));
                ratio_S *= f_shifted / division_protect(f_base);
            }
            if(average(baseRadiance)==0.0){
                continue;
            }
            if(average(shiftedRadiance)==0.0){
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                    continue;
            }

            { // test whether the segment connect to light is blocked by the dynamic objects
                Ray shadow_ray{vLight.position, dir_light, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(vLight.position, base_sampledVLight.position), shiftedTimeStamp};
                if (occluded_dynamic(scene, shadow_ray))
                {
                    splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                    continue;
                }
            }

            //update the ghost intersection of the shifted_sampledVLight
            PathVertex shifted_sampledVLight(base_sampledVLight);
            shifted_sampledVLight.time = shiftedTimeStamp;
            updateGhostIsectInfo(shifted_sampledVLight, shifted_vLight, scene, shiftedTimeStamp);
            Real shifted_pdf;
            Real shifted_weight = computeMixMISWeight_better_dynamic(shifted_vLight, shifted_sensorPath[1], shifted_sampledVLight, shifted_sampledVSensor,
                                                    shifted_lightPath, shifted_sensorPath, s, 1, shifted_pathState,
                                                    NULL, SamplingTechnique::DYNAMIC_FROM_SENSOR, scene, shiftedTimeStamp, shifted_pdf);
            shiftedRadiance *= shifted_weight;
            { // the bsdf value at the last vertex; don't record this in the ratio_S
                Vector3 base_wi = normalize(prevLight.position - vLight.position);
                Vector3 shifted_wi = normalize(shifted_vLight.position - vLight.position);
                Vector3 wo = normalize(base_sampledVLight.position - vLight.position);

                const Material &matVlight = scene.materials[vLight.material_id];
                Spectrum f_base = eval(matVlight, base_wi, wo, vLight, scene.texture_pool) / fabs(dot(wo, vLight.geometric_normal));
                Spectrum f_shifted = eval(matVlight, shifted_wi, wo, vLight, scene.texture_pool) / fabs(dot(wo, vLight.geometric_normal));
                shiftedRadiance *= f_shifted / division_protect(f_base);
            }
            assert(average(baseRadiance) == average(baseRadiance));
            assert(average(shiftedRadiance) == average(shiftedRadiance));
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
            splatSample(img, shifted_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));

        }else if(status == SHIFT_MAP_STATUS::CONNECTED){//use the cached ratio
            if(average(ratio_S)>10.)
            {
                splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                continue;
            }
            { // test whether the segment connect to light is blocked by the dynamic objects
                Ray shadow_ray{vLight.position, dir_light, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(vLight.position, base_sampledVLight.position), shiftedTimeStamp};
                if (occluded_dynamic(scene, shadow_ray))
                {
                    splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                    continue;
                }
            }
            //update the ghost intersection of the shifted_sampledVLight
            PathVertex shifted_sampledVLight(base_sampledVLight);
            shifted_sampledVLight.time = shiftedTimeStamp;
            updateGhostIsectInfo(shifted_sampledVLight, shifted_vLight, scene, shiftedTimeStamp);
            Real shifted_pdf;
            Real shifted_weight = computeMixMISWeight_better_dynamic(shifted_vLight, shifted_sensorPath[1], shifted_sampledVLight, shifted_sampledVSensor,
                                                    shifted_lightPath, shifted_sensorPath, s, 1, shifted_pathState,
                                                    NULL, SamplingTechnique::DYNAMIC_FROM_SENSOR, scene, shiftedTimeStamp, shifted_pdf);
            Spectrum shiftedRadiance = shifted_weight * unweighted_baseRadiance * ratio_S;
            assert(average(baseRadiance) == average(baseRadiance));
            assert(average(shiftedRadiance) == average(shiftedRadiance));
            splatSample(img, base_pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
            splatSample(img, shifted_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));
        
        }else{
            assert(false);
        }        
    }
}

// a faster DP version to get rid of any redundant computation
void compute_shifted_dynamic_path_3_two_ends_dp(const PathVertex &vDynamic,
                                                const PathVertex &shiftedVDynamic,
                                                std::vector<PathVertex> &lightPath,
                                                std::vector<PathVertex> &sensorPath,
                                                pcg32_state &rng,
                                                pcg32_state &shifted_rng, //copy of the previous base rng after sampling the dynamic position and dir to the light
                                                const Scene &scene,
                                                const int &numSamples,
                                                int baseTimeStamp,
                                                int shiftedTimeStamp,
                                                Image3 &img)
{
    PathVertex lastHitLight;
    // /////////////////////////////////
    // // 0 stands for directly hit the dynamic instead of sampling it explicitly,
    // // which is not possible to hit the exactly same point on the surface, pdf=0
    
    if (lightPath.size() < 2 || sensorPath.size() < 2)
        return;
    
    if(lightPath.back().type == VertexType::Light){
        lastHitLight = lightPath.back();
        //then remove it from the list 
        lightPath.pop_back();
        assert(lastHitLight.type == VertexType::Light);
    }
    int nLight = lightPath.size(), nSensor = sensorPath.size();
    const PathVertex &vDynamicS = lightPath[0], &vDynamicT = sensorPath[0];

    ////////////// step 2 let's compose the shfited paths first
    std::vector<PathVertex> shifted_lightPath, shifted_sensorPath;
    std::vector<SHIFT_MAP_STATUS> shifted_status_light, shifted_status_sensor;
    
    //prepare the shiftedVDynamicS and shiftedVDynamicT
    PathVertex shiftedVDynamicS(shiftedVDynamic), shiftedVDynamicT(shiftedVDynamic);
    Vector3 shifted_dynamic_dir2Light, shifted_dynamic_dir2Sensor;
    {
        sampleDynamicDir(shiftedVDynamicS, shifted_dynamic_dir2Light, shifted_rng, scene, shiftedTimeStamp);
        shifted_lightPath.push_back(shiftedVDynamicS);
        shifted_status_light.push_back(SHIFT_MAP_STATUS::DIVERGED);

        Real pdfDir2Sensor;
        bool isOtherDir = sampleScattering(shiftedVDynamicT, shifted_dynamic_dir2Light, shifted_dynamic_dir2Sensor, pdfDir2Sensor, shifted_rng, scene);
        if (isOtherDir)
        {
            shiftedVDynamicT.transport.pdfDir = pdfDir2Sensor;
            shifted_sensorPath.push_back(shiftedVDynamicT);
            shifted_status_sensor.push_back(SHIFT_MAP_STATUS::DIVERGED);
        }else{
            assert(false);
        }
    }

    ///////////////////////////////////////
    
    pcg32_state shifted_rng_cache = shifted_rng;
    SHIFT_MAP_STATUS status = SHIFT_MAP_STATUS::DIVERGED;
    for (int s = 2; s <= nLight; ++s){//the shifted light path 
        if (status == SHIFT_MAP_STATUS::FAILED) break;
        const PathVertex &vLight = lightPath[s - 1]; //the base vertex
        PathVertex &shifted_prevLight = shifted_lightPath[s - 2]; //the shifted vertex
        const PathVertex& prevLight = lightPath[s-2]; 
        if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED || status == SHIFT_MAP_STATUS::CONNECTED){ //the main difference is when calculating the radiance, not here
            if(is_vertex_dynamic(vLight, scene)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }

            //here already after the reconnection, using the same subpath from now on
            // retest the visiblity of the new segment
            Ray shadow_ray{vLight.position, normalize(prevLight.position - vLight.position), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(vLight.position, prevLight.position), shiftedTimeStamp};
            if (occluded_dynamic(scene, shadow_ray))
            {
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            //otherwise, update hte ghost intersection information and pdfFwd of shifted_vLight is actually different because of shiftedPath[s-3]
            PathVertex shifted_vLight(vLight); 
            shifted_vLight.time = shiftedTimeStamp;
            assert(s >= 3); //after the reconnection, so s is at least 3
            if(!updatePdfsOfSegment(shifted_lightPath[s-3], shifted_prevLight, shifted_vLight, s<nLight? lightPath[s]:PathVertex{}, scene, shiftedTimeStamp)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            updateGhostIsectInfo(shifted_vLight, shifted_prevLight, scene, shiftedTimeStamp);
            shifted_lightPath.push_back(shifted_vLight);
            if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED) 
                status = SHIFT_MAP_STATUS::CONNECTED;
            shifted_status_light.push_back(status);

        }else if(connectability_condition(vLight, prevLight, shifted_prevLight, scene)){ //getVertexRoughnessType(vLight, scene) == MATType::VERTEX_TYPE_DIFFUSE){
            
            assert(status == SHIFT_MAP_STATUS::DIVERGED);
            if(is_vertex_dynamic(vLight, scene)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            bool isConnected = connectVertices(shifted_prevLight, vLight, scene, shiftedTimeStamp);
            if (!isConnected){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            //if the reconnection succeed,
            // 1) need to retest the ghost intersections for the reconnected segment
            // 2) need to recompute the pdfFwd and pdfRev //TODO check what others are changed??
            PathVertex shifted_vLight(vLight);
            shifted_vLight.time = shiftedTimeStamp;
            if(s==2){
                if(!updatePdfsOfSegment(shiftedVDynamicS, shifted_vLight, s<nLight? lightPath[s]:PathVertex{}, scene, shiftedTimeStamp)){
                    status = SHIFT_MAP_STATUS::FAILED;
                    break;
                }
                assert(shifted_vLight.transport.pdfFwd == shifted_vLight.transport.pdfFwd);
            } 
            else{
                if(!updatePdfsOfSegment(shifted_lightPath[s-3], shifted_prevLight, shifted_vLight, s<nLight? lightPath[s]:PathVertex{}, scene, shiftedTimeStamp)){
                    status = SHIFT_MAP_STATUS::FAILED;
                    break;
                }
            }//TODO if s==2 need to recompute the pdfs; light side should be fine; but sensor side is different
            updateGhostIsectInfo(shifted_vLight, shifted_prevLight, scene, shiftedTimeStamp);
            shifted_lightPath.push_back(shifted_vLight);
            status = SHIFT_MAP_STATUS::RECENTLY_CONNECTED;
            shifted_status_light.push_back(status);
            assert(shifted_lightPath.back().transport.pdfFwd==shifted_lightPath.back().transport.pdfFwd);
            
        }else{//if not reconnectable; just do random walk using the same random seed
            
            assert(status == SHIFT_MAP_STATUS::DIVERGED);
            Vector3 shifted_nextDir;
            Real pdfNextDir;
            if(s==2){
                pdfNextDir = shiftedVDynamicS.transport.pdfDir;
                assert(pdfNextDir>0);
                shifted_nextDir = shifted_dynamic_dir2Light;
                RayDifferential ray_diff;
                randomWalk(shifted_lightPath, make_const_spectrum(1.0), pdfNextDir, shifted_nextDir, ray_diff, s+1, s-1, TransportDirection::TO_LIGHT, shifted_rng, scene, shiftedTimeStamp);
            }else{
                bool isOtherDir = sampleScattering(shifted_prevLight, normalize(shifted_lightPath[s-3].position - shifted_prevLight.position),
                    shifted_nextDir, pdfNextDir, shifted_rng, scene); 
                if(!isOtherDir){
                    status = SHIFT_MAP_STATUS::FAILED;
                    break;
                }
                RayDifferential ray_diff; //during randomWalk the new vertex is inserted and its with the correct ghost intersection
                randomWalk(shifted_lightPath, shifted_prevLight.transport.contrib, pdfNextDir, shifted_nextDir, ray_diff, s+1, s-1, TransportDirection::TO_LIGHT, shifted_rng, scene, shiftedTimeStamp);
            }
            {//TODO something to do with it.
                PathVertex lastHitLight;
                if(shifted_lightPath.back().type == VertexType::Light){
                    lastHitLight = shifted_lightPath.back();
                    //then remove it from the list 
                    shifted_lightPath.pop_back();
                    assert(lastHitLight.type == VertexType::Light);
                }
            }
            if(shifted_lightPath.size()< s){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            shifted_status_light.push_back(status);
            assert(shifted_lightPath.back().transport.pdfFwd==shifted_lightPath.back().transport.pdfFwd);
        }

        assert(shifted_lightPath.size() == shifted_status_light.size());
    }

    status = SHIFT_MAP_STATUS::DIVERGED; //reset the status
    shifted_rng = shifted_rng_cache;
    // the sensor side
    for (int t = 2; t <= nSensor; ++t)
    {
        if (status == SHIFT_MAP_STATUS::FAILED) break;
        const PathVertex &vSensor = sensorPath[t - 1];
        PathVertex &shifted_prevSensor = shifted_sensorPath[t - 2];
        const PathVertex& prevSensor = sensorPath[t-2];

        if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED || status == SHIFT_MAP_STATUS::CONNECTED){
            if(is_vertex_dynamic(vSensor, scene)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }

            //retest the visibility of the new segment
            Ray shadow_ray{vSensor.position, normalize(prevSensor.position - vSensor.position), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(vSensor.position, prevSensor.position), shiftedTimeStamp};
            if (occluded_dynamic(scene, shadow_ray))
            {
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            PathVertex shifted_vSensor(vSensor);
            shifted_vSensor.time = shiftedTimeStamp;
            assert(t>=3); //after the reconnection, so t is at least 3
            if(!updatePdfsOfSegment(shifted_sensorPath[t-3], shifted_prevSensor, shifted_vSensor, t<nSensor? sensorPath[t]:PathVertex{}, scene, shiftedTimeStamp)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            updateGhostIsectInfo(shifted_vSensor, shifted_prevSensor, scene, shiftedTimeStamp);
            shifted_sensorPath.push_back(shifted_vSensor);
            if(status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED)
                status = SHIFT_MAP_STATUS::CONNECTED;
            shifted_status_sensor.push_back(status);

        }else if(connectability_condition(vSensor, prevSensor, shifted_prevSensor, scene)){ //getVertexRoughnessType(vSensor, scene) == MATType::VERTEX_TYPE_DIFFUSE){
            
            assert(status == SHIFT_MAP_STATUS::DIVERGED);
            if(is_vertex_dynamic(vSensor, scene)){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            bool isConnected = connectVertices(shifted_prevSensor, vSensor, scene, shiftedTimeStamp);
            if (!isConnected){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            PathVertex shifted_vSensor(vSensor);
            shifted_vSensor.time = shiftedTimeStamp;
            if(t==2){ 
                if(!updatePdfsOfSegment(shiftedVDynamicT, shifted_vSensor, t<nSensor? sensorPath[t]:PathVertex{}, scene, shiftedTimeStamp)){
                    status = SHIFT_MAP_STATUS::FAILED;
                    break;
                }
            }else{
                if(!updatePdfsOfSegment(shifted_sensorPath[t-3], shifted_prevSensor, shifted_vSensor, t<nSensor? sensorPath[t]:PathVertex{}, scene, shiftedTimeStamp)){
                    status = SHIFT_MAP_STATUS::FAILED;
                    break;
                }
            }
            updateGhostIsectInfo(shifted_vSensor, shifted_prevSensor, scene, shiftedTimeStamp);
            shifted_sensorPath.push_back(shifted_vSensor);
            status = SHIFT_MAP_STATUS::RECENTLY_CONNECTED;
            shifted_status_sensor.push_back(status);

        }else{//if not reconnectable; just do random walk using the same random seed
            
            assert(status == SHIFT_MAP_STATUS::DIVERGED);
            Vector3 shifted_nextDir;
            Real pdfNextDir;
            if(t==2){
                pdfNextDir = shiftedVDynamicT.transport.pdfDir;
                shifted_nextDir = shifted_dynamic_dir2Sensor;
                RayDifferential ray_diff;
                randomWalk(shifted_sensorPath, make_const_spectrum(1.0), pdfNextDir, shifted_nextDir, ray_diff, t+1, t-1, TransportDirection::TO_VIEW, shifted_rng, scene, shiftedTimeStamp);
            }else{
                bool isOtherDir = sampleScattering(shifted_prevSensor, normalize(shifted_sensorPath[t-3].position - shifted_prevSensor.position), 
                    shifted_nextDir, pdfNextDir, shifted_rng, scene);
                if(!isOtherDir){
                    status = SHIFT_MAP_STATUS::FAILED;
                    break;
                }
                RayDifferential ray_diff;
                randomWalk(shifted_sensorPath, shifted_prevSensor.transport.contrib, pdfNextDir, shifted_nextDir, ray_diff, t+1, t-1, TransportDirection::TO_VIEW, shifted_rng, scene, shiftedTimeStamp);
            }
            assert(shifted_sensorPath.back().type != VertexType::Light);
            if(shifted_sensorPath.size()<t){
                status = SHIFT_MAP_STATUS::FAILED;
                break;
            }
            shifted_status_sensor.push_back(status);
        }

        assert(shifted_sensorPath.size() == shifted_status_sensor.size());
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////// step3 actually calculate the difference radiance ///////////////////////
    
    int nLight_shifted = shifted_lightPath.size(), nSensor_shifted = shifted_sensorPath.size();

    const Material &matDynamic = scene.materials[vDynamicS.material_id];
    const Vector3 dynamic_dir2Light = normalize(lightPath[1].position - vDynamicS.position);
    const Vector3 dynamic_dir2Sensor = normalize(sensorPath[1].position - vDynamicT.position);
    Spectrum base_f_dynamic = eval(matDynamic, dynamic_dir2Sensor, dynamic_dir2Light, vDynamicS, scene.texture_pool) / fabs(dot(dynamic_dir2Light, vDynamicS.geometric_normal));

    ///////////////// prepare all the ratios needed
    // some cached ratio 
    Spectrum ratio_S = make_const_spectrum(1.0);
    Spectrum ratio_T = make_const_spectrum(1.0);
    Spectrum ratio_middle = make_const_spectrum(1.0); //only the middle f

    ////////////// the middle part: ratio_middle until hit both the reconnected vertices
    status = SHIFT_MAP_STATUS::DIVERGED; //reset the status
    if(nLight_shifted < 2 || nSensor_shifted < 2)
        status = SHIFT_MAP_STATUS::FAILED;

    if (status != SHIFT_MAP_STATUS::FAILED)
    {
        const Vector3 shifted_dynamic_dir2Light = normalize(shifted_lightPath[1].position - shiftedVDynamic.position);
        const Vector3 shifted_dynamic_dir2Sensor = normalize(shifted_sensorPath[1].position - shiftedVDynamic.position);
        Spectrum shifted_f_dynamic = eval(matDynamic, shifted_dynamic_dir2Sensor, shifted_dynamic_dir2Light, shiftedVDynamic, scene.texture_pool) / fabs(dot(shifted_dynamic_dir2Light, shiftedVDynamic.geometric_normal));
        ratio_middle *= shifted_f_dynamic;// / base_f_dynamic;

        const Vector3 dynamic_dir2Light = normalize(lightPath[1].position - vDynamic.position);
        const Vector3 dynamic_dir2Sensor = normalize(sensorPath[1].position - vDynamic.position);
        ratio_middle *= fabs(dot(shifted_dynamic_dir2Light, shiftedVDynamic.geometric_normal) );
        ratio_middle *= fabs(dot(shifted_dynamic_dir2Sensor, shiftedVDynamic.geometric_normal));
    }

    std::vector<Spectrum> base_Llight_cache;
    std::vector<Spectrum> base_Wsensor_cache;
    std::vector<PathVertex> base_sampledVLight_cache;
    std::vector<PathVertex> base_sampledVSensor_cache;
    std::vector<Vector2> base_pRaster2_cache;

    std::vector<Spectrum> shifted_Llight_cache;
    std::vector<Spectrum> shifted_Wsensor_cache;
    std::vector<PathVertex> shifted_sampledVLight_cache;
    std::vector<PathVertex> shifted_sampledVSensor_cache;
    std::vector<Vector2> shifted_pRaster2_cache;

    std::vector<SHIFT_MAP_STATUS> light_sm_status_cache;
    std::vector<SHIFT_MAP_STATUS> sensor_sm_status_cache;

    SHIFT_MAP_STATUS light_sub_status = status; //wether shifted_dynamicS and shifted_dynamicT are successfully generated?
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
        
        // 2) render the shifted light subpath
        if (light_sub_status == SHIFT_MAP_STATUS::FAILED)
        {
            shifted_Llight_cache.push_back(make_zero_spectrum());
            light_sm_status_cache.push_back(light_sub_status);
            shifted_sampledVLight_cache.push_back({});
            continue;
        }
        
        if(s > nLight_shifted) light_sub_status = SHIFT_MAP_STATUS::FAILED;
        else light_sub_status = shifted_status_light[s-1];

        if (light_sub_status == SHIFT_MAP_STATUS::FAILED)
        {
            shifted_Llight_cache.push_back(make_zero_spectrum());
            light_sm_status_cache.push_back(light_sub_status);
            shifted_sampledVLight_cache.push_back({});
            continue;
        }

        const PathVertex &shifted_vLight = shifted_lightPath[s - 1];
        const PathVertex &shifted_prevLight = shifted_lightPath[s - 2];
        Spectrum shifted_Llight = make_zero_spectrum();

        if( light_sub_status == SHIFT_MAP_STATUS::DIVERGED){
            //connect to the sampled light source
            PathVertex shifted_sampledVLight(sampledVLight);
            //record the ratio here for the diverged part
            ratio_S = shifted_vLight.transport.contrib / vLight.transport.contrib;

            shifted_Llight = directIlluminationGivenLight(shifted_vLight, shifted_sampledVLight, scene, shiftedTimeStamp, true, true); 
            {//f at the last vertex connected to the light source
                const PathVertex &shifted_prevLight = shifted_lightPath[s - 2];
                const Material &mat = scene.materials[shifted_vLight.material_id];
                Vector3 dir_in = normalize(shifted_prevLight.position - shifted_vLight.position);
                Vector3 dir_out = normalize(shifted_sampledVLight.position - shifted_vLight.position);
                Spectrum f = eval(mat, dir_in, dir_out, shifted_vLight, scene.texture_pool) / fabs(dot(dir_out, shifted_vLight.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
                shifted_Llight *= shifted_vLight.transport.contrib * f;
            }
            shifted_Llight_cache.push_back(shifted_Llight);
            shifted_sampledVLight_cache.push_back(shifted_sampledVLight);
            light_sm_status_cache.push_back(light_sub_status);

        }else if(light_sub_status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED){
            
            //update the ratio_S due to the reconnection
            ////G at the reconnection point 
            {           
                const Vector3 base_dir_reconnected = normalize(vLight.position - prevLight.position);
                const Vector3 shifted_dir_reconnected = normalize(vLight.position - shifted_prevLight.position);
                Real base_G_S = fabs(dot(base_dir_reconnected, vLight.geometric_normal)) / distance_squared(prevLight.position, vLight.position);
                Real shifted_G_S = fabs(dot(shifted_dir_reconnected, vLight.geometric_normal)) 
                    / distance_squared(shifted_prevLight.position, vLight.position);
                ratio_S *= shifted_G_S / base_G_S;
            }

            if(average(ratio_S)>10.)
            {
                shifted_Llight_cache.push_back(make_zero_spectrum());
                light_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
                shifted_sampledVLight_cache.push_back({});
                continue;
            }

            if(s>=3){ // ratio of f at the prev vertex before the reconnected vertex
                const Material &mat = scene.materials[prevLight.material_id];
                Vector3 wi = normalize(lightPath[s - 3].position - prevLight.position);
                Vector3 wo = normalize(vLight.position - prevLight.position);
                Spectrum base_prev_f = eval(mat, wi, wo, prevLight, scene.texture_pool) / fabs(dot(wo, prevLight.geometric_normal));
                const Material &shifted_mat = scene.materials[shifted_prevLight.material_id];
                wi = normalize(shifted_lightPath[s-3].position - shifted_prevLight.position);
                wo = normalize(vLight.position - shifted_prevLight.position);
                Spectrum shifted_prev_f = eval(shifted_mat, wi, wo, shifted_prevLight, scene.texture_pool) / fabs(dot(wo, shifted_prevLight.geometric_normal));
                ratio_S *= shifted_prev_f /division_protect(base_prev_f);
            }

            { // retest the visibility of this NEE segment
                const PathVertex &v1 = vLight;
                const PathVertex &v2 = sampledVLight;
                Ray shadow_ray{v1.position, normalize(v2.position - v1.position), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(v2.position, v1.position), shiftedTimeStamp};
                if (occluded_dynamic(scene, shadow_ray))
                {
                    shifted_Llight_cache.push_back(make_zero_spectrum());
                    light_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
                    shifted_sampledVLight_cache.push_back({});
                    continue;
                }
            }

            {//difference of the f of the last vertex next to light source
                const Material &matS = scene.materials[vLight.material_id];
                const Vector3 shifted_wi = normalize(shifted_prevLight.position - vLight.position);
                const Vector3 base_wi = normalize(prevLight.position - vLight.position);
                const Vector3 wo = normalize(sampledVLight.position - vLight.position);
                Spectrum shifted_f = eval(matS, shifted_wi, wo, vLight, scene.texture_pool) / fabs(dot(wo, vLight.geometric_normal));
                Spectrum base_f = eval(matS, base_wi, wo, vLight, scene.texture_pool) / fabs(dot(wo, vLight.geometric_normal));
                shifted_Llight = Llight * ratio_S * shifted_f / division_protect(base_f) ; //TODO only record L in light then later to assembly; then only need to multiply with shifted_f here. 
            }

            //retest the ghost intersection of the shifted_sampledVLight
            PathVertex shifted_sampledVLight(sampledVLight);
            shifted_sampledVLight.time = shiftedTimeStamp;
            updateGhostIsectInfo(shifted_sampledVLight, shifted_vLight, scene, shiftedTimeStamp);
            
            shifted_Llight_cache.push_back(shifted_Llight);
            shifted_sampledVLight_cache.push_back(shifted_sampledVLight);
            light_sm_status_cache.push_back(light_sub_status);

             // the ratio of f at the reconnected vertex changed       
            if(s<nLight && s < nLight_shifted){
                const Material &matS = scene.materials[vLight.material_id];
                const Vector3 shifted_wi = normalize(shifted_prevLight.position - vLight.position);
                const Vector3 base_wi = normalize(prevLight.position - vLight.position);
                const Vector3 wo = normalize(lightPath[s].position - vLight.position);
                Spectrum shifted_f = eval(matS, shifted_wi, wo, vLight, scene.texture_pool) / fabs(dot(wo, vLight.geometric_normal));
                Spectrum base_f = eval(matS, base_wi, wo, vLight, scene.texture_pool) / fabs(dot(wo, vLight.geometric_normal));
                ratio_S *= shifted_f /division_protect(base_f);
            }//for later use

        }else if(light_sub_status == SHIFT_MAP_STATUS::CONNECTED){
            if(average(ratio_S)>10.)
            {
                shifted_Llight_cache.push_back(make_zero_spectrum());
                light_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
                shifted_sampledVLight_cache.push_back({});
                continue;
            }

            { // retest the visiblity of the NEE segment
                const PathVertex &v1 = vLight;
                const PathVertex &v2 = sampledVLight;
                Ray shadow_ray{v1.position, normalize(v2.position - v1.position), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(v2.position, v1.position), shiftedTimeStamp};
                if (occluded_dynamic(scene, shadow_ray))
                {
                    shifted_Llight_cache.push_back(make_zero_spectrum());
                    light_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
                    shifted_sampledVLight_cache.push_back({});
                    continue;
                }
            }
            shifted_Llight = Llight * ratio_S;
            //retest the ghost intersection of the shifted_sampledVLight
            PathVertex shifted_sampledVLight(sampledVLight);
            shifted_sampledVLight.time = shiftedTimeStamp;
            updateGhostIsectInfo(shifted_sampledVLight, shifted_vLight, scene, shiftedTimeStamp);
            shifted_Llight_cache.push_back(shifted_Llight);
            shifted_sampledVLight_cache.push_back(shifted_sampledVLight);
            light_sm_status_cache.push_back(light_sub_status);

        }else{
            assert(false);
        }
        
    }

    ////////////////////////////// for the sensor side
    SHIFT_MAP_STATUS sensor_sub_status = status; // reset

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

         // 2) render the shfited sensor subpath
        if(t > nSensor_shifted) sensor_sub_status = SHIFT_MAP_STATUS::FAILED;
        else sensor_sub_status = shifted_status_sensor[t-1];

        if (sensor_sub_status == SHIFT_MAP_STATUS::FAILED)
        {
            shifted_Wsensor_cache.push_back(make_zero_spectrum());
            sensor_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
            shifted_sampledVSensor_cache.push_back({});
            shifted_pRaster2_cache.push_back({});
            continue;
        }
        const PathVertex &shifted_vSensor = shifted_sensorPath[t - 1];
        const PathVertex &shifted_prevSensor = shifted_sensorPath[t - 2];

        Spectrum shifted_Wsensor = make_zero_spectrum();

        if( sensor_sub_status == SHIFT_MAP_STATUS::DIVERGED){
           
           PathVertex shfted_sampledVSensor;
            Vector3 dirToSensor;
            Vector2 shifted_pRaster2{0.0, 0.0};
            Spectrum shifted_Wsensor = connectToSensorVertex(shifted_vSensor, shfted_sampledVSensor, dirToSensor, shifted_pRaster2, scene, shiftedTimeStamp, true);
            //f at the last vertex connected to the sensor
            const Material &mat = scene.materials[shifted_vSensor.material_id];
                Vector3 dir_in = normalize(shifted_prevSensor.position - shifted_vSensor.position);
            Spectrum f = eval(mat, dir_in, dirToSensor, shifted_vSensor, scene.texture_pool) / fabs(dot(dirToSensor, shifted_vSensor.geometric_normal)); // cosine included; because the G inside has no cosineAtVertex
            shifted_Wsensor *= shifted_vSensor.transport.contrib * f;
            assert(average(shifted_Wsensor)==average(shifted_Wsensor));
            shifted_Wsensor_cache.push_back(shifted_Wsensor);
            sensor_sm_status_cache.push_back(sensor_sub_status);
            shifted_pRaster2_cache.push_back(shifted_pRaster2);
            shifted_sampledVSensor_cache.push_back(shfted_sampledVSensor); //TODO right into this info package

            ratio_T = shifted_vSensor.transport.contrib / vSensor.transport.contrib;
        
        }else if(sensor_sub_status == SHIFT_MAP_STATUS::RECENTLY_CONNECTED){
            //update the ratio_T due to the reconnection
            //G at the reconnection point
            {           
                const Vector3 base_dir_reconnected = normalize(vSensor.position - prevSensor.position);
                const Vector3 shifted_dir_reconnected = normalize(vSensor.position - shifted_prevSensor.position);
                Real base_G_T = fabs(dot(base_dir_reconnected, vSensor.geometric_normal)) 
                    / distance_squared(prevSensor.position, vSensor.position);
                Real shifted_G_T = fabs(dot(shifted_dir_reconnected, vSensor.geometric_normal)) 
                    / distance_squared(shifted_prevSensor.position, vSensor.position);
                ratio_T *= safe_divide(shifted_G_T, base_G_T);

            }
            if(average(ratio_T)>10.)
            {
                shifted_Wsensor_cache.push_back(make_zero_spectrum());
                sensor_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
                shifted_sampledVSensor_cache.push_back({});
                shifted_pRaster2_cache.push_back({});
                continue;
            }

            if(t>=3){ // ratio of f at the prev vertex before the reconnected vertex
                const Material &mat = scene.materials[prevSensor.material_id];
                Vector3 wi = normalize(sensorPath[t - 3].position - prevSensor.position);
                Vector3 wo = normalize(vSensor.position - prevSensor.position);
                Spectrum base_prev_f = eval(mat, wi, wo, prevSensor, scene.texture_pool) / fabs(dot(wo, prevSensor.geometric_normal));
                const Material &shifted_mat = scene.materials[shifted_prevSensor.material_id];
                wi = normalize(shifted_sensorPath[t-3].position - shifted_prevSensor.position);
                wo = normalize(vSensor.position - shifted_prevSensor.position);
                Spectrum shifted_prev_f = eval(shifted_mat, wi, wo, shifted_prevSensor, scene.texture_pool) / fabs(dot(wo, shifted_prevSensor.geometric_normal));
                ratio_T *= shifted_prev_f / division_protect(base_prev_f);

            }

            {  // retest visibility of the NEE segment
                const PathVertex &v1 = vSensor;
                const PathVertex &v2 = base_sampledVSensor;
                Ray shadow_ray{v1.position, normalize(v2.position - v1.position), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(v2.position, v1.position), shiftedTimeStamp};
                if (occluded_dynamic(scene, shadow_ray))
                {
                    shifted_Wsensor_cache.push_back(make_zero_spectrum());
                    sensor_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
                    shifted_sampledVSensor_cache.push_back({});
                    shifted_pRaster2_cache.push_back({});
                    continue;
                }
            }

            {//difference of the f of the last vertex next to sensor source
                 const Material &matT = scene.materials[vSensor.material_id];
                const Vector3 shifted_wi = normalize(shifted_prevSensor.position - vSensor.position);
                const Vector3 base_wi = normalize(prevSensor.position - vSensor.position);
                const Vector3 wo = normalize(base_sampledVSensor.position - vSensor.position);
                Spectrum shifted_f = eval(matT, shifted_wi, wo, vSensor, scene.texture_pool) / fabs(dot(wo, vSensor.geometric_normal));
                Spectrum base_f = eval(matT, base_wi, wo, vSensor, scene.texture_pool) / fabs(dot(wo, vSensor.geometric_normal));
                shifted_Wsensor = Wsensor * ratio_T * shifted_f /division_protect(base_f); //TODO only record We in Wsensor then later to assembly; then only need to multiply with shifted_f here.
            }
            PathVertex shfted_sampledVSensor(base_sampledVSensor);
            shfted_sampledVSensor.time = shiftedTimeStamp;
            updateGhostIsectInfo(shfted_sampledVSensor, shifted_vSensor, scene, shiftedTimeStamp);

            shifted_Wsensor_cache.push_back(shifted_Wsensor);
            sensor_sm_status_cache.push_back(sensor_sub_status);
            shifted_pRaster2_cache.push_back(base_pRaster2);
            shifted_sampledVSensor_cache.push_back(shfted_sampledVSensor); //TODO right into this info package

            // the ratio of f at the reconnected vertex changed
            if(t<nSensor && t < nSensor_shifted){
                const Material &matT = scene.materials[vSensor.material_id];
                const Vector3 shifted_wi = normalize(shifted_prevSensor.position - vSensor.position);
                const Vector3 base_wi = normalize(prevSensor.position - vSensor.position);
                const Vector3 wo = normalize(sensorPath[t].position - vSensor.position);
                Spectrum shifted_f = eval(matT, shifted_wi, wo, vSensor, scene.texture_pool) / fabs(dot(wo, vSensor.geometric_normal));
                Spectrum base_f = eval(matT, base_wi, wo, vSensor, scene.texture_pool) / fabs(dot(wo, vSensor.geometric_normal));
                ratio_T *= shifted_f / division_protect(base_f);

            }//for later use

        }else if(sensor_sub_status == SHIFT_MAP_STATUS::CONNECTED){
            if(average(ratio_T)>10.)
            {
                shifted_Wsensor_cache.push_back(make_zero_spectrum());
                sensor_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
                shifted_sampledVSensor_cache.push_back({});
                shifted_pRaster2_cache.push_back({});
                continue;
            }
            
            { // retest visibility of the NEE segment
                const PathVertex &v1 = vSensor;
                const PathVertex &v2 = base_sampledVSensor;
                Ray shadow_ray{v1.position, normalize(v2.position - v1.position), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(v2.position, v1.position), shiftedTimeStamp};
                if (occluded_dynamic(scene, shadow_ray))
                {
                    shifted_Wsensor_cache.push_back(make_zero_spectrum());
                    sensor_sm_status_cache.push_back(SHIFT_MAP_STATUS::FAILED);
                    shifted_sampledVSensor_cache.push_back({});
                    shifted_pRaster2_cache.push_back({});
                    continue;
                }
            }
            shifted_Wsensor = Wsensor * ratio_T;
            //retest the ghost intersection of the shifted_sampledSensor
            PathVertex shfted_sampledVSensor(base_sampledVSensor);
            shfted_sampledVSensor.time = shiftedTimeStamp;
            updateGhostIsectInfo(shfted_sampledVSensor, shifted_vSensor, scene, shiftedTimeStamp);
            shifted_Wsensor_cache.push_back(shifted_Wsensor);
            sensor_sm_status_cache.push_back(sensor_sub_status);
            shifted_pRaster2_cache.push_back(base_pRaster2);
            shifted_sampledVSensor_cache.push_back(shfted_sampledVSensor); //TODO right into this info package
        }else{
            assert(false);
        }
        assert(average(ratio_T)==average(ratio_T));
    }

    ///////////////////// for mis weight computation ///////////////////////
    // prepare the accumulated pdfFwd and pdfRevs for both subpaths
    DSubPathState lightPathState(nLight), sensorPathState(nSensor);
    cumprod_pdfs(lightPath, lightPathState);
    cumprod_pdfs(sensorPath, sensorPathState);

    //for the shifted paths
    DSubPathState shifted_lightPathState(nLight_shifted), shifted_sensorPathState(nSensor_shifted);
    cumprod_pdfs(shifted_lightPath, shifted_lightPathState);
    cumprod_pdfs(shifted_sensorPath, shifted_sensorPathState);

    //////////////////doing the assembly here ////////////////////////////////
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
            baseRadiance *= base_weight;
            if (average(baseRadiance) == 0.0)
                continue;

            // render the shifted path
            const SHIFT_MAP_STATUS &light_sub_status = light_sm_status_cache[s_cache];
            const SHIFT_MAP_STATUS &sensor_sub_status = sensor_sm_status_cache[t_cache];
            const Spectrum &shifted_Llight = shifted_Llight_cache[s_cache];
            const Spectrum &shifted_Wsensor = shifted_Wsensor_cache[t_cache];
            const Vector2 &shifted_pRaster2 = shifted_pRaster2_cache[t_cache];

            if (light_sub_status == SHIFT_MAP_STATUS::FAILED || sensor_sub_status == SHIFT_MAP_STATUS::FAILED)
            {
                splatSample(img, pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                continue;
            }

            if (average(shifted_Llight) == 0.0 || average(shifted_Wsensor) == 0.0)
            { // not really neccessary, just for safety
                splatSample(img, pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, true));
                continue;
            }

            const PathVertex &shifted_vLight = shifted_lightPath[s - 1];
            const PathVertex &shifted_vSensor = shifted_sensorPath[t - 1];
            Real shifted_pdf;
            Real shifted_weight = computeMixMISWeight_better_dynamic(shifted_vLight, shifted_vSensor, 
                                                                    shifted_sampledVLight_cache[s_cache], 
                                                                    shifted_sampledVSensor_cache[t_cache],
                                                                    shifted_lightPath, shifted_sensorPath, s, t,
                                                                    shifted_lightPathState, shifted_sensorPathState, 
                                                                    SamplingTechnique::DYNAMIC_TWO_ENDS, scene, shiftedTimeStamp, shifted_pdf);

            Spectrum shiftedRadiance = shifted_Llight * shifted_Wsensor * vDynamicS.transport.contrib / (Real)numSamples; 
            shiftedRadiance /= shiftedVDynamicS.transport.pdfDir;  //if depends on the reconnect position
            shiftedRadiance /= shiftedVDynamicT.transport.pdfDir;
            shiftedRadiance *= ratio_middle * shifted_weight; // using the base MIS weight now.
            splatSample(img, pRaster2, baseRadiance * compute_mis_weight_frame(baseTimeStamp, false, base_pdf, shifted_pdf));
            splatSample(img, shifted_pRaster2, shiftedRadiance * compute_mis_weight_frame(shiftedTimeStamp, false, base_pdf, shifted_pdf));
            
        }
    }
}
/**
 * /////TODO attention to use half the spp here to amortize the sampling in both frames
 * @brief render the finite difference only for the solid dynamic part

 */
