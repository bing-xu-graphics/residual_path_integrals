#pragma once

#include "scene.h"
#include "pcg.h"
#include "../utils.h"

/// for ghost portal sampling
bool sampleGhost(PathVertex &vertex, // output
                 Vector3 &diro,      // output
                 pcg32_state &rng,
                 const Scene &scene,
                 int time)
{
    Real identity_w = next_pcg32_real<Real>(rng);
    Vector2 uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real w = next_pcg32_real<Real>(rng);
    int dynamic_id = sample_dynamic_shape(scene, identity_w);
    const Dynamic &dyn_agent = scene.dynamics[dynamic_id];
    const Shape &ghost = (time == NEW_TIME) ? dyn_agent.old_shape : dyn_agent.base_shape;
    PointAndNormal point_on_ghost =
        sample_point_on_shape(ghost, Vector3{}, uv, w);
    Real pdfAGhost = dynamic_pmf(scene, dynamic_id) *
                     pdf_point_on_shape(ghost, point_on_ghost, Vector3{});
    // without loss of generality, let's always sample hemisphere instead of sphere

    vertex.type = VertexType::Ghost;
    vertex.time = time;
    vertex.dynamic_id = dynamic_id;
    vertex.shape_id = dyn_agent.shape_id;
    vertex.position = point_on_ghost.position;       // p1;
    vertex.geometric_normal = point_on_ghost.normal; // TODO still added
    // no geometry_normal or any other shading information//TODO check nothing here
    //////////////////////// ///for shading or sampling ///TODO maybe remove if not sampling the dir
    vertex.primitive_id = point_on_ghost.primitive_id;
    vertex.st = point_on_ghost.st;
    ShadingInfo shadingInfo = compute_shading_info(ghost, vertex);
    vertex.shading_frame = shadingInfo.shading_frame;
    vertex.uv = shadingInfo.uv;
    Real pdfGhostDir;
    Vector2 dynamic_dir_rdn{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Vector3 dirFromGhost = UniformSampleSphere(dynamic_dir_rdn);
    pdfGhostDir = UniformSpherePdf();

    // output
    diro = dirFromGhost;
    vertex.transport.pdfOriginA = pdfAGhost;
    vertex.transport.pdfDir = pdfGhostDir;
    vertex.transport.pdfFwd = pdfGhostDir;
    vertex.transport.pathLength = 0; // TODO

    vertex.transport.contrib = make_const_spectrum(1.0 / pdfAGhost);
    vertex.hit_ghost_times = 1;
    return true;
}

////below not used now
/*
    The ghost counterpart of the dynamic portal sampling
*/

Real computeGhostPortalMISWeight(const PathVertex *sampledVLight,
                                 const PathVertex *vLight,
                                 const PathVertex *vSensor,
                                 const PathVertex *sampledVSensor)

{
    int hit_ghost_times = (vLight ? vLight->hit_ghost_times : 0) + (vSensor ? vSensor->hit_ghost_times : 0); // each can be both sides
    hit_ghost_times += (sampledVLight ? sampledVLight->hit_ghost_times : 0) + (sampledVSensor ? sampledVSensor->hit_ghost_times : 0);

    if (hit_ghost_times == 0)
    {
        return 0.0;
    }

#ifdef USE_PT_AS_4_GHOST          // using path tracing as the forth sampling strategy
    hit_ghost_times += 1; // TO remove just for testing
#endif
    return 1 / (Real)(hit_ghost_times); // TODO not sure how many times it hit ghost right?? only know it hits. can't decide 1 or 2; only way is to guarantee no penetration
}

void cumprod_pdfs_except_origin(const std::vector<PathVertex> &subPath, DSubPathState &pathState)
{
    assert(pathState.pdf_pairs.size() == subPath.size());
    int n = pathState.pdf_pairs.size();
    Real pdfFwd_product(1.0), pdfRev_product(1.0);
    for (int i = 0; i < n; ++i)
    {
        pdfFwd_product *= subPath[i].transport.pdfFwd;
        pathState.pdf_pairs[i].first = pdfFwd_product;
    }
    for (int i = n - 1; i >= 0; --i)
    {
        pdfRev_product *= subPath[i].transport.pdfRev;
        pathState.pdf_pairs[i].second = pdfRev_product;
    }
}
Real pdfInnerGhostEdge(const Vector3& ghostPos, const Vector3& ghostNormal, const Real& pdfDir, const PathVertex& endV1, const PathVertex& endV2){
    Vector3 vecGhost = endV1.position - endV2.position;
    Real dist = length(vecGhost);
    vecGhost /= dist;
    Real G_ghostEdge = fabs(dot(vecGhost, endV1.geometric_normal) * dot(vecGhost, endV2.geometric_normal)) / dist / dist;
    return pdfDir / fabs(dot(ghostNormal, vecGhost)) * G_ghostEdge;
}

// the ratio of the pdfs so that to change the pdfA of ghost point to pdfA of the actual solid point by connecting the ghost point and the end light/sensor
Real pdfGhostToSolidNext(const Vector3& ghostPos, const Vector3& ghostNormal, const PathVertex& lightOrSensor, const PathVertex& solidVertex){
    Vector3 dirGhostEdge = normalize(lightOrSensor.position - ghostPos);
    Real G_ghost = fabs(dot(ghostNormal, dirGhostEdge)) / distance_squared(ghostPos, lightOrSensor.position);
    Real G_next = fabs(dot(solidVertex.geometric_normal, dirGhostEdge)) / distance_squared(solidVertex.position, lightOrSensor.position);
    return G_next / G_ghost;
}

Real computeGhostPortalMISWeight_better(const PathVertex &vGhost,
                                        const PathVertex &vLight, const PathVertex &vSensor,
                                        const PathVertex &sampledVLight, const PathVertex &sampledVSensor,
                                        const std::vector<PathVertex> &lightPath, const std::vector<PathVertex> &sensorPath,
                                        int S, int T, const DSubPathState &lightSubpathState, const DSubPathState &sensorSubpathState,
                                        const SamplingTechnique &technique, const Scene &scene, int time)
{
    if (technique == SamplingTechnique::GHOST_FROM_LIGHT)
    {
        if (vSensor.hit_ghost_times + sampledVSensor.hit_ghost_times == 0.0)
            return 0.0;
        const DSubPathState &subpathState = sensorSubpathState;
        assert(lightPath.empty());
        // pdfCurrHit: pdf of sampling this path by starting from this ghost object
        Real pdfANextGhost = vGhost.transport.pdfOriginA * pdfGhostToSolidNext(vGhost.position, vGhost.geometric_normal, sampledVLight, sensorPath[1]);
        Real pdfCurrHit = pdfANextGhost;
        pdfCurrHit *= sampledVLight.transport.pdfOriginA;
        pdfCurrHit *= subpathState.pdf_pairs[T - 1].first;
        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        Real pdfSum(0.);             // this will be included in the ghost edge below (the first edge)
        for (int t = 2; t <= T; ++t) // the first one t==1 is light source
        {                            // TODO the last one can just be a envmap void vertex
            // test whether segment [t-1]-[t-2] hits ghost dynamic or not
            int num_ghost_hits = sensorPath[t - 1].hit_ghost_times - sensorPath[t - 2].hit_ghost_times;
            if (num_ghost_hits <= 0)
                continue; // this segment never hits a ghost obj
            // this segment passes the ghost object //or that sensorPath[t-1].ghost_ids.size()>0
            std::vector<int> ghost_shape_ids = sensorPath[t - 1].ghost_ids;
            assert(num_ghost_hits == ghost_shape_ids.size());
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                assert(ghost_shape_id >= 0);
                Vector3 ghostPos = sensorPath[t - 1].ghost_hit_positions[i];
                Vector3 ghostNormal = sensorPath[t - 1].ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                // the light and sensor subpath part
                Real pdfLightPath = safe_divide(subpathState.pdf_pairs[1].second, subpathState.pdf_pairs[t - 2].second);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(subpathState.pdf_pairs[T - 1].first, subpathState.pdf_pairs[t - 1].first);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // pdf to sample the first vertex on lightPath
                if (t == 2)
                {
                    // pdfA of sampling the first solid vertex by directly connecting from the light source to the ghost
                    pdfThisHit *= pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVLight, sensorPath[1]);
                }
                else // t>2
                {    // the pdf to generate the two solid vertices on the ghost edge at the same time
                    Real pdfGhostEdge = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, sensorPath[t - 1], sensorPath[t -2]);
                    pdfThisHit *= pdfGhostEdge;
                }
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum > 0);
            }
            assert(pdfSum > 0);
        }
        // also need to handle the last segment (last vertex->sensor)
        {
            std::vector<int> ghost_shape_ids = sampledVSensor.ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = sampledVSensor.ghost_hit_positions[i];
                Vector3 ghostNormal = sampledVSensor.ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                Real pdfThisHit = pdfAreaThisHit;
                //pdfA to sample the first solid vertex next to the ghost
                pdfThisHit *= pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVSensor, sensorPath[T - 1]); 
                pdfThisHit *= safe_divide(subpathState.pdf_pairs[1].second, subpathState.pdf_pairs[T - 1].second);
                if(T>2){ // the pdfRev at sensorPath[T-2] is possibly changed because we use the sampledVsensor instead the other vertex along the way
                    Vector3 dir_in = normalize(sampledVSensor.position - sensorPath[T - 1].position); 
                    Vector3 dir_out = normalize(sensorPath[T-2].position- sensorPath[T - 1].position);
                    Real pdfWBsdf = pdf_sample_bsdf(scene.materials[sensorPath[T-1].material_id], dir_in, dir_out, sensorPath[T-1], scene.texture_pool, TransportDirection::TO_LIGHT);
                    Real pdfABsdf = pdf_convert_w_to_A(pdfWBsdf, sensorPath[T-1], sensorPath[T-2]);
                    pdfThisHit *= pdfABsdf / sensorPath[T-2].transport.pdfRev;//should be in area measure, but saved it.
                }
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum > 0);
            }
        }
#ifdef USE_PT_AS_4_GHOST
        Real pdfPT = safe_divide(subpathState.pdf_pairs[1].second, subpathState.pdf_pairs[T-1].second);
        if (T >= 3 && pdfPT > 0)
        {
            Real pdfWRevTmp = pdf_bsdf_at_vertex(sensorPath[T - 1], sensorPath[T - 2], sampledVSensor, TransportDirection::TO_LIGHT, scene);
            Real pdfARevTmp = pdf_convert_w_to_A(pdfWRevTmp, sensorPath[T - 1], sensorPath[T - 2]);
            pdfPT *= pdfARevTmp / sensorPath[T - 2].transport.pdfRev;
        }
        pdfPT *= sampledVSensor.transport.pdfOriginA; // the pdf to sample the first vertex from image/sensor sampling; ==cameraPdfA
        pdfPT *= sampledVLight.transport.pdfOriginA;
        pdfSum += MIS(pdfPT);
        assert(pdfPT == pdfPT);
#endif
        return MIS(pdfCurrHit) / pdfSum;
    }
    else if (technique == SamplingTechnique::GHOST_FROM_SENSOR)
    {
        if (vLight.hit_ghost_times + sampledVLight.hit_ghost_times == 0.0)
            return 0.0;
        const DSubPathState &subpathState = lightSubpathState;
        Real pdfCurrHit = vGhost.transport.pdfOriginA * pdfGhostToSolidNext(vGhost.position, vGhost.geometric_normal, sampledVSensor, lightPath[1]);
        pdfCurrHit *= subpathState.pdf_pairs[S - 1].first;
        pdfCurrHit *= sampledVLight.transport.pdfOriginA;
        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        Real pdfSum(0.);
        for (int s = 2; s <= S; ++s)
        {
            int num_ghost_hits = lightPath[s - 1].hit_ghost_times - lightPath[s - 2].hit_ghost_times;
            if (num_ghost_hits <= 0)
                continue;
            std::vector<int> ghost_shape_ids = lightPath[s - 1].ghost_ids;
            assert(num_ghost_hits == ghost_shape_ids.size());
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = lightPath[s - 1].ghost_hit_positions[i];
                Vector3 ghostNormal = lightPath[s - 1].ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                // the light and sensor subpath part
                Real pdfLightPath = safe_divide(subpathState.pdf_pairs[S - 1].first, subpathState.pdf_pairs[s - 1].first);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(subpathState.pdf_pairs[1].second, subpathState.pdf_pairs[s - 2].second);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // pdf to sample the first vertex on lightPath
                if (s == 2)
                {
                   pdfThisHit *= pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVSensor, lightPath[1]);
                }
                else{
                    Real pdfGhostEdge = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, lightPath[s - 1], lightPath[s - 2]);
                    pdfThisHit *= pdfGhostEdge;
                }
                pdfSum += MIS(pdfThisHit);
            }
        }
        {// also need to handle the last segment (last vertex->light)
            std::vector<int> ghost_shape_ids = sampledVLight.ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = sampledVLight.ghost_hit_positions[i];
                Vector3 ghostNormal = sampledVLight.ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                Real pdfThisHit = pdfAreaThisHit * pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVLight, lightPath[S - 1]);
                pdfThisHit *= safe_divide(subpathState.pdf_pairs[1].second, subpathState.pdf_pairs[S - 1].second);
                if(S>2){//the pdfRev at lightPath[S-2] is possibly changed because we use the new sampledVLight.
                    Vector3 dir_in = normalize(sampledVLight.position - lightPath[S - 1].position); 
                    Vector3 dir_out = normalize(lightPath[S-2].position- lightPath[S - 1].position);
                    Real pdfWBsdf = pdf_sample_bsdf(scene.materials[lightPath[S-1].material_id], dir_in, dir_out, lightPath[S-1], scene.texture_pool, TransportDirection::TO_VIEW);
                    Real pdfABsdf = pdf_convert_w_to_A(pdfWBsdf, lightPath[S-1], lightPath[S-2]);
                    pdfThisHit *= pdfABsdf / lightPath[S-2].transport.pdfRev;
                }
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum > 0);
            }
        }
#ifdef USE_PT_AS_4_GHOST
    // add in the uni-directional pt sampling pdf
        Real pdfPT = subpathState.pdf_pairs[S - 1].first;
        pdfPT *= sampledVSensor.transport.pdfOriginA ;
        pdfPT *= pdfGhostToSolidNext(vGhost.position, vGhost.geometric_normal, sampledVSensor, lightPath[1]);
        pdfPT *= sampledVLight.transport.pdfOriginA;
        pdfSum += MIS(pdfPT);
#endif
        return MIS(pdfCurrHit) / pdfSum;
    }
    else if (technique == SamplingTechnique::GHOST_TWO_ENDS)
    {
        if (vLight.hit_ghost_times + vSensor.hit_ghost_times + sampledVLight.hit_ghost_times + sampledVSensor.hit_ghost_times == 0)
            return 0.0;
        Real pdfCurrHit = vGhost.transport.pdfOriginA;
        pdfCurrHit *= sampledVLight.transport.pdfOriginA;
        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        pdfCurrHit *= lightSubpathState.pdf_pairs[S - 1].first;
        pdfCurrHit *= sensorSubpathState.pdf_pairs[T - 1].first;
        // pdf to sample the first segment dir etc..
        Real pdfGhostEdge = pdfInnerGhostEdge(vGhost.position, vGhost.geometric_normal, vGhost.transport.pdfDir, sensorPath[0], sensorPath[1]);
        pdfCurrHit *= pdfGhostEdge;
        Real pdfSum(0.);
        for (int s = 3; s <= S; ++s) // here starts from 3
        {
            int num_ghost_hits = lightPath[s - 1].hit_ghost_times - lightPath[s - 2].hit_ghost_times;
            if (num_ghost_hits <= 0)
                continue;
            std::vector<int> ghost_shape_ids = lightPath[s - 1].ghost_ids;
            assert(num_ghost_hits == ghost_shape_ids.size());
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = lightPath[s - 1].ghost_hit_positions[i];
                Vector3 ghostNormal = lightPath[s - 1].ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                // the light and sensor subpath part
                Real pdfLightPath = safe_divide(lightSubpathState.pdf_pairs[S - 1].first, lightSubpathState.pdf_pairs[s - 1].first);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(lightSubpathState.pdf_pairs[1].second, lightSubpathState.pdf_pairs[s - 2].second);
                pdfSensorPath *= sensorSubpathState.pdf_pairs[T - 1].first;
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                { // make up for the original hit ghost segment; the pdfA to generate sensorPath[1] from sensorPath[0]
                    Vector3 dir_in = normalize(sensorPath[2].position - lightPath[0].position);
                    Vector3 dir_out = normalize(sensorPath[1].position - sensorPath[0].position);
                    Real pdfWNextV2Sensor = pdf_sample_bsdf(scene.materials[sensorPath[0].material_id], dir_in, dir_out, sensorPath[0], scene.texture_pool, TransportDirection::TO_VIEW);
                    Real pdfANextV2Sensor = pdf_convert_w_to_A(pdfWNextV2Sensor, sensorPath[0], sensorPath[1]);
                    pdfThisHit *= pdfANextV2Sensor;
                } 
                // pdf to sample the first vertex on lightPath
                Real pdfGhostEdge = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, lightPath[s - 1], lightPath[s - 2]);
                pdfThisHit *= pdfGhostEdge;
                pdfSum += MIS(pdfThisHit);
            }
        }
        for (int t = 2; t <= T; ++t) // the first one t==1 is light source// here starts from 2
        {
            int num_ghost_hits = sensorPath[t - 1].hit_ghost_times - sensorPath[t - 2].hit_ghost_times;
            if (num_ghost_hits <= 0)
                continue; // this segment never hits a ghost obj
            // this segment passes the ghost object //or that sensorPath[t-1].ghost_ids.size()>0
            std::vector<int> ghost_shape_ids = sensorPath[t - 1].ghost_ids;
            assert(num_ghost_hits == ghost_shape_ids.size());
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = sensorPath[t - 1].ghost_hit_positions[i];
                Vector3 ghostNormal = sensorPath[t - 1].ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                // the light and sensor subpath part
                Real pdfLightPath = safe_divide(sensorSubpathState.pdf_pairs[1].second, sensorSubpathState.pdf_pairs[t - 2].second);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                pdfLightPath *= lightSubpathState.pdf_pairs[S - 1].first;

                Real pdfSensorPath = safe_divide(sensorSubpathState.pdf_pairs[T - 1].first, sensorSubpathState.pdf_pairs[t - 1].first);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                if(t>2){ // make up for the original hit ghost segment; the pdfA to generate sensorPath[0] from sensorPath[1]
                    Vector3 dir_in = normalize(sensorPath[2].position - sensorPath[1].position);
                    Vector3 dir_out = normalize(sensorPath[0].position - sensorPath[1].position);
                    Real pdfWNextV2Light = pdf_sample_bsdf(scene.materials[sensorPath[1].material_id], dir_in, dir_out, sensorPath[1], scene.texture_pool, TransportDirection::TO_LIGHT);
                    Real pdfANextV2Light = pdf_convert_w_to_A(pdfWNextV2Light, sensorPath[1], sensorPath[0]);
                    pdfThisHit *= pdfANextV2Light;
                }
                // pdf to sample the first vertex on lightPath
                Real pdfGhostEdgeThis = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, sensorPath[t - 1], sensorPath[t - 2]);
                pdfThisHit *= pdfGhostEdgeThis;
                pdfSum += MIS(pdfThisHit);
            }

        }
        {
            std::vector<int> ghost_shape_ids = sampledVLight.ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = sampledVLight.ghost_hit_positions[i];
                Vector3 ghostNormal = sampledVLight.ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                Real pdfThisHit = pdfAreaThisHit * pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVLight, lightPath[S - 1]); 
                { // make up for the original hit ghost segment; the pdfA to generate sensorPath[1] from sensorPath[0]
                    Vector3 dir_in = normalize((S>=3? lightPath[2]:sampledVLight).position - sensorPath[0].position);
                    Vector3 dir_out = normalize(sensorPath[1].position - sensorPath[0].position);
                    Real pdfWNextV2Sensor = pdf_sample_bsdf(scene.materials[sensorPath[0].material_id], dir_in, dir_out, sensorPath[0], scene.texture_pool, TransportDirection::TO_VIEW);
                    Real pdfANextV2Sensor = pdf_convert_w_to_A(pdfWNextV2Sensor, sensorPath[0], sensorPath[1]);
                    pdfThisHit *= pdfANextV2Sensor;
                } 
                pdfThisHit *= safe_divide(lightSubpathState.pdf_pairs[1].second, lightSubpathState.pdf_pairs[S - 1].second);
                if(S>2){//the pdfRev at lightPath[S-2] is possibly changed because we use the new sampledVLight.
                    Vector3 dir_in = normalize(sampledVLight.position - lightPath[S - 1].position); 
                    Vector3 dir_out = normalize(lightPath[S-2].position- lightPath[S - 1].position);
                    Real pdfWBsdf = pdf_sample_bsdf(scene.materials[lightPath[S-1].material_id], dir_in, dir_out, lightPath[S-1], scene.texture_pool, TransportDirection::TO_VIEW);
                    Real pdfABsdf = pdf_convert_w_to_A(pdfWBsdf, lightPath[S-1], lightPath[S-2]);
                    pdfThisHit *= pdfABsdf / lightPath[S-2].transport.pdfRev;
                }
                pdfThisHit *= sensorSubpathState.pdf_pairs[T - 1].first;
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
            }
        }
        {
            std::vector<int> ghost_shape_ids = sampledVSensor.ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = sampledVSensor.ghost_hit_positions[i];
                Vector3 ghostNormal = sampledVSensor.ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                Real pdfThisHit = pdfAreaThisHit * pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVSensor, sensorPath[T - 1]);  //pdfA to sample the first solid vertex next to the ghost
                { // make up for the original hit ghost segment; the pdfA to generate sensorPath[0] from sensorPath[1]
                    Vector3 dir_in = normalize((T>=3? sensorPath[2]:sampledVSensor).position - sensorPath[1].position);
                    Vector3 dir_out = normalize(sensorPath[0].position - sensorPath[1].position);
                    Real pdfWNextV2Light = pdf_sample_bsdf(scene.materials[sensorPath[1].material_id], dir_in, dir_out, sensorPath[1], scene.texture_pool, TransportDirection::TO_LIGHT);
                    Real pdfANextV2Light = pdf_convert_w_to_A(pdfWNextV2Light, sensorPath[1], sensorPath[0]);
                    pdfThisHit *= pdfANextV2Light;
                }
                pdfThisHit *= safe_divide(sensorSubpathState.pdf_pairs[1].second, sensorSubpathState.pdf_pairs[T - 1].second);
                if(T>2){ // the pdfRev at sensorPath[T-2] is possibly changed because we use the sampledVsensor instead the other vertex along the way
                    Vector3 dir_in = normalize(sampledVSensor.position - sensorPath[T - 1].position); 
                    Vector3 dir_out = normalize(sensorPath[T-2].position- sensorPath[T - 1].position);
                    Real pdfWBsdf = pdf_sample_bsdf(scene.materials[sensorPath[T-1].material_id], dir_in, dir_out, sensorPath[T-1], scene.texture_pool, TransportDirection::TO_LIGHT);
                    Real pdfABsdf = pdf_convert_w_to_A(pdfWBsdf, sensorPath[T-1], sensorPath[T-2]);
                    pdfThisHit *= pdfABsdf / sensorPath[T-2].transport.pdfRev;//should be in area measure, but saved it.
                }
                pdfThisHit *= lightSubpathState.pdf_pairs[S - 1].first;
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
            }
        }
#ifdef USE_PT_AS_4_GHOST
        Real pdfPT = lightSubpathState.pdf_pairs[S - 1].first * safe_divide(sensorSubpathState.pdf_pairs[1].second, sensorSubpathState.pdf_pairs[T - 1].second);
        if (T >= 3 && pdfPT > 0)
        {
            Real pdfWRevTmp = pdf_bsdf_at_vertex(sensorPath[T - 1], sensorPath[T - 2], sampledVSensor, TransportDirection::TO_LIGHT, scene);
            Real pdfARevTmp = pdf_convert_w_to_A(pdfWRevTmp, sensorPath[T - 1], sensorPath[T - 2]);
            pdfPT *= pdfARevTmp / sensorPath[T - 2].transport.pdfRev;
        }
        pdfPT *= sampledVSensor.transport.pdfOriginA;
        pdfPT *= sampledVLight.transport.pdfOriginA;
        { // make up for the original hit ghost segment; the pdfA to generate sensorPath[0] from sensorPath[1]
            Vector3 dir_in = normalize((T>=3? sensorPath[2]:sampledVSensor).position - sensorPath[1].position);
            Vector3 dir_out = normalize(sensorPath[0].position - sensorPath[1].position);
            Real pdfWNextV2Light = pdf_sample_bsdf(scene.materials[sensorPath[1].material_id], dir_in, dir_out, sensorPath[1], scene.texture_pool, TransportDirection::TO_LIGHT);
            Real pdfANextV2Light = pdf_convert_w_to_A(pdfWNextV2Light, sensorPath[1], sensorPath[0]);
            pdfPT *= pdfANextV2Light;
        } 
        pdfSum += MIS(pdfPT);
#endif
        return MIS(pdfCurrHit) / pdfSum;
    }
    else if (technique == UNI_PT)
    {
        Real pdfPT = sampledVSensor.transport.pdfOriginA * lightSubpathState.pdf_pairs[S - 1].first * sampledVLight.transport.pdfOriginA;
        Real pdfSum = MIS(pdfPT);
        {//the first sensor segment 
            std::vector<int> ghost_shape_ids = lightPath[0].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = lightPath[0].ghost_hit_positions[i];
                Vector3 ghostNormal = lightPath[0].ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                Real pdfThisHit = pdfAreaThisHit * pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVSensor, lightPath[0]); //pdfA to sample the first solid vertex next to the ghost
                pdfThisHit *= lightSubpathState.pdf_pairs[S-1].first;
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
            }
        }
        for(int s=2; s<=S; ++s)
        {
            int num_ghost_hits = lightPath[s - 1].hit_ghost_times - lightPath[s - 2].hit_ghost_times;
            if (num_ghost_hits <= 0)
                continue;
            std::vector<int> ghost_shape_ids = lightPath[s - 1].ghost_ids;
            assert(num_ghost_hits == ghost_shape_ids.size());
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = lightPath[s - 1].ghost_hit_positions[i];
                Vector3 ghostNormal = lightPath[s - 1].ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                // the light and sensor subpath part
                Real pdfLightPath = safe_divide(lightSubpathState.pdf_pairs[S - 1].first, lightSubpathState.pdf_pairs[s - 1].first);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(lightSubpathState.pdf_pairs[0].second, lightSubpathState.pdf_pairs[s - 2].second);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                Real pdfGhostEdge = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, lightPath[s - 1], lightPath[s - 2]);
                pdfThisHit *= pdfGhostEdge;
                pdfSum += MIS(pdfThisHit);
            }
        }
        {//the light segment
            std::vector<int> ghost_shape_ids = sampledVLight.ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = sampledVLight.ghost_hit_positions[i];
                Vector3 ghostNormal = sampledVLight.ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                Real pdfThisHit = pdfAreaThisHit * pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVLight, lightPath[S - 1]); 
                pdfThisHit *= safe_divide(lightSubpathState.pdf_pairs[0].second, lightSubpathState.pdf_pairs[S - 1].second);
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
            }
        }

        return MIS(pdfPT) / pdfSum;
    }
    else
    {
        assert(false);
        return 1.0;
    }
}

// Real compute_mis_ghost_part(const PathVertex *sampledVLight,
//                             const PathVertex *vLight,
//                             const PathVertex *vSensor,
//                             const PathVertex *sampledVSensor)

Real compute_mis_ghost_part(const PathVertex &sampledVLight,
                            const PathVertex &vLight,
                            const PathVertex &vSensor,
                            const PathVertex &sampledVSensor,
                            const SamplingTechnique &technique)

{
    int hit_ghost_times = 0;
    if(technique == SamplingTechnique::GHOST_FROM_LIGHT){
        hit_ghost_times = vSensor.hit_ghost_times + sampledVSensor.hit_ghost_times;
    }else if(technique == SamplingTechnique::GHOST_FROM_SENSOR){
        hit_ghost_times = vLight.hit_ghost_times + sampledVLight.hit_ghost_times;
    }else if(technique == SamplingTechnique::GHOST_TWO_ENDS){
        hit_ghost_times = vLight.hit_ghost_times + vSensor.hit_ghost_times 
            + sampledVLight.hit_ghost_times + sampledVSensor.hit_ghost_times;
    }else{
        assert(false);
    }

    int hit_dynamic_times = vLight.hit_dynamic_times + vSensor.hit_dynamic_times;
    int total_times = hit_ghost_times + hit_dynamic_times;
    return total_times == 0 ? 0.0 : 1 / (Real)(total_times);      // TODO not sure how many times it hit ghost right?? only know it hits. can't decide 1 or 2; only way is to guarantee no penetration
}
