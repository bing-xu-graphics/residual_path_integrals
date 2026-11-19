#pragma once

#include "scene.h"
#include "pcg.h"
#include "../utils.h"
#include "portal_utils_ghost.h"

#define USE_PT_AS_4
#define USE_PT_AS_4_GHOST
#define MIS_MIX 

#define NEE_ONE_SAMPLE_PROB 0.5

Real computePdfAsGhostPath(const PathVertex &vLight, const PathVertex &vSensor,
                                   const PathVertex &sampledVLight, const PathVertex &sampledVSensor,
                                   const std::vector<PathVertex> &lightPath, const std::vector<PathVertex> &sensorPath,
                                   int S, int T, const DSubPathState &lightSubpathState, const DSubPathState &sensorSubpathState,
                                   const SamplingTechnique &technique, const Scene &scene, int time)
{
    
    if (technique == SamplingTechnique::DYNAMIC_FROM_LIGHT)
    {
        const DSubPathState &subpathState = sensorSubpathState;
        Real pdfSum(0.);
        {// the first segment (light->first vertex)
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
                Real pdfThisHit = pdfAreaThisHit * pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVLight, sensorPath[0]);
                pdfThisHit *= subpathState.pdf_pairs[T-1].first;
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum > 0);
            }
        }
        for (int t = 2; t <= T; ++t) // the first one t==1 is the dynamic object
        {                            // TODO the last one can just be a envmap void vertex
            std::vector<int> ghost_shape_ids = sensorPath[t - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue; // this segment never hits a ghost obj
            // this segment passes the ghost object //or that sensorPath[t-1].ghost_ids.size()>0
            
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
                Real pdfLightPath = safe_divide(subpathState.pdf_pairs[0].second, subpathState.pdf_pairs[t - 2].second);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(subpathState.pdf_pairs[T - 1].first, subpathState.pdf_pairs[t - 1].first);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // the pdf to generate the two solid vertices on the ghost edge at the same time
                Real pdfGhostEdge = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, sensorPath[t - 1], sensorPath[t - 2]);
                pdfThisHit *= pdfGhostEdge;
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum > 0);
            }
            assert(pdfSum > 0);
        }
        {// also need to handle the last segment (last vertex->sensor)
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
                pdfThisHit *= safe_divide(subpathState.pdf_pairs[0].second, subpathState.pdf_pairs[T - 1].second);
                if(T>=2){ //TODO check this 1111 the pdfRev at sensorPath[T-2] is possibly changed because we use the sampledVsensor instead the other vertex along the way
                    Vector3 dir_in = normalize(sampledVSensor.position - sensorPath[T - 1].position); 
                    Vector3 dir_out = normalize(sensorPath[T-2].position- sensorPath[T - 1].position);
                    Real pdfWBsdf = pdf_sample_bsdf(scene.materials[sensorPath[T-1].material_id], dir_in, dir_out, sensorPath[T-1], scene.texture_pool, TransportDirection::TO_LIGHT);
                    Real pdfABsdf = pdf_convert_w_to_A(pdfWBsdf, sensorPath[T-1], sensorPath[T-2]);
                    pdfThisHit *= pdfABsdf / sensorPath[T-2].transport.pdfRev;//should be in area measure, but saved it.
                }
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum == pdfSum);
            }
        }
        assert(pdfSum == pdfSum);
        return pdfSum;
    }
    else if(technique == SamplingTechnique::DYNAMIC_FROM_SENSOR)
    {
        const DSubPathState &subpathState = lightSubpathState;
        Real pdfSum(0.);
        {//handle the first segment connected to sensor
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
                pdfThisHit *= pdfGhostToSolidNext(ghostPos, ghostNormal, sampledVSensor, lightPath[0]); 
                pdfThisHit *= subpathState.pdf_pairs[S - 1].first;
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum > 0);
            }
        }
        for (int s = 2; s <= S; ++s)
        {
            std::vector<int> ghost_shape_ids = lightPath[s - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue;
            
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
                Real pdfSensorPath = safe_divide(subpathState.pdf_pairs[0].second, subpathState.pdf_pairs[s - 2].second);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // pdf to sample the first vertex on lightPath
                Real pdfGhostEdge = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, lightPath[s - 1], lightPath[s - 2]);
                pdfThisHit *= pdfGhostEdge;
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
                pdfThisHit *= safe_divide(subpathState.pdf_pairs[0].second, subpathState.pdf_pairs[S - 1].second);
                if(S>=2){//TODO 1111 check this the pdfRev at lightPath[S-2] is possibly changed because we use the new sampledVLight. 
                    Vector3 dir_in = normalize(sampledVLight.position - lightPath[S - 1].position); 
                    Vector3 dir_out = normalize(lightPath[S-2].position- lightPath[S - 1].position);
                    Real pdfWBsdf = pdf_sample_bsdf(scene.materials[lightPath[S-1].material_id], dir_in, dir_out, lightPath[S-1], scene.texture_pool, TransportDirection::TO_VIEW);
                    Real pdfABsdf = pdf_convert_w_to_A(pdfWBsdf, lightPath[S-1], lightPath[S-2]);
                    pdfThisHit *= pdfABsdf / lightPath[S-2].transport.pdfRev;
                }
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
            }
        }
        return pdfSum;
    }
    else if (technique == SamplingTechnique::DYNAMIC_TWO_ENDS)
    {
        Real pdfSum(0.);
        for (int s = 2; s <= S; ++s)  
        {
            std::vector<int> ghost_shape_ids = lightPath[s - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue;

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
                pdfSensorPath *= sensorSubpathState.pdf_pairs[T - 1].first;
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                Real pdfGhostEdge = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, lightPath[s - 1], lightPath[s - 2]);
                pdfThisHit *= pdfGhostEdge;
                pdfSum += MIS(pdfThisHit);
            }
        }
        for (int t = 2; t <= T; ++t) // the first one t==1 is light source// here starts from 2
        {
            std::vector<int> ghost_shape_ids = sensorPath[t - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue; // this segment never hits a ghost obj
            // this segment passes the ghost object //or that sensorPath[t-1].ghost_ids.size()>0
            
            for (int i = 0; i < num_ghost_hits; ++i)
            {
                int ghost_shape_id = ghost_shape_ids[i];
                Vector3 ghostPos = sensorPath[t - 1].ghost_hit_positions[i];
                Vector3 ghostNormal = sensorPath[t - 1].ghost_hit_geom_normals[i];
                PointAndNormal pointAndNormal{ghostPos, ghostNormal};
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_ghost(pdfAreaThisHit, pdfDirThisHit, ghost_shape_id, pointAndNormal, scene, time);
                // the light and sensor subpath part
                Real pdfLightPath = safe_divide(sensorSubpathState.pdf_pairs[0].second, sensorSubpathState.pdf_pairs[t - 2].second);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                pdfLightPath *= lightSubpathState.pdf_pairs[S - 1].first;
                Real pdfSensorPath = safe_divide(sensorSubpathState.pdf_pairs[T - 1].first, sensorSubpathState.pdf_pairs[t - 1].first);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // pdf to sample the first vertex on lightPath
                Real pdfGhostEdgeThis = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, sensorPath[t - 1], sensorPath[t - 2]);
                pdfThisHit *= pdfGhostEdgeThis;
                //modify the edge2light next to the vDynamic and match the corresponding sampling strategy
                {                                                                                                                                       // t==T; only the edge on the light side needs to be modified.
                    Real pdfBSDF_nextDynamicSEdge = pdf_bsdf_at_vertex(lightPath[0], lightPath[1], sensorPath[1], TransportDirection::TO_LIGHT, scene); // the ratio should be in area measure, but the G is the same, so skip the G
                    pdfThisHit *= safe_divide(pdfBSDF_nextDynamicSEdge, lightPath[0].transport.pdfDir); 
                }

                pdfSum += MIS(pdfThisHit);
            }

        }
        {//when the ghost object is on the end segment connected to sampledVLight
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
                pdfThisHit *= safe_divide(lightSubpathState.pdf_pairs[0].second, lightSubpathState.pdf_pairs[S - 1].second) * sensorSubpathState.pdf_pairs[T - 1].first;
                if(S>=2){//the pdfRev at lightPath[S-2] is possibly changed because we use the new sampledVLight. //TODO 1111 check this
                    Vector3 dir_in = normalize(sampledVLight.position - lightPath[S - 1].position); 
                    Vector3 dir_out = normalize(lightPath[S-2].position- lightPath[S - 1].position);
                    Real pdfWBsdf = pdf_sample_bsdf(scene.materials[lightPath[S-1].material_id], dir_in, dir_out, lightPath[S-1], scene.texture_pool, TransportDirection::TO_VIEW);
                    Real pdfABsdf = pdf_convert_w_to_A(pdfWBsdf, lightPath[S-1], lightPath[S-2]);
                    pdfThisHit *= pdfABsdf / lightPath[S-2].transport.pdfRev;
                }
                //no need for modification here
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
            }
        }
        {//when the ghost object is on the end segment connected to sampledVSensor
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
                //modify the edge2light next to the vDynamic and match the corresponding sampling strategy
                {                                                                                                                                       // t==T; only the edge on the light side needs to be modified.
                    Real pdfBSDF_nextDynamicSEdge = pdf_bsdf_at_vertex(lightPath[0], lightPath[1], sensorPath[1], TransportDirection::TO_LIGHT, scene); // the ratio should be in area measure, but the G is the same, so skip the G
                    pdfThisHit *= safe_divide(pdfBSDF_nextDynamicSEdge, lightPath[0].transport.pdfDir);
                }
                pdfThisHit *= safe_divide(sensorSubpathState.pdf_pairs[0].second, sensorSubpathState.pdf_pairs[T - 1].second) * lightSubpathState.pdf_pairs[S - 1].first;
                if(T>=2){ // the pdfRev at sensorPath[T-2] is possibly changed because we use the sampledVsensor instead the other vertex along the way  //TODO 1111 check this
                    Vector3 dir_in = normalize(sampledVSensor.position - sensorPath[T - 1].position); 
                    Vector3 dir_out = normalize(sensorPath[T-2].position- sensorPath[T - 1].position);
                    Real pdfWBsdf = pdf_sample_bsdf(scene.materials[sensorPath[T-1].material_id], dir_in, dir_out, sensorPath[T-1], scene.texture_pool, TransportDirection::TO_LIGHT);
                    Real pdfABsdf = pdf_convert_w_to_A(pdfWBsdf, sensorPath[T-1], sensorPath[T-2]);
                    pdfThisHit *= pdfABsdf / sensorPath[T-2].transport.pdfRev;//should be in area measure, but saved it.
                }
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
            }
        }
        return pdfSum;
    }
    else if (technique == UNI_PT)
    {
        Real pdfSum = 0.;
        {//the ghost intersection is on the first sensor segment 
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
        {//the ghost intersection is on the middle segment
            std::vector<int> ghost_shape_ids = lightPath[s - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue;
            
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
        {//the ghost intersection is on the light segment
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
        return pdfSum;
    }
    else{
        assert(false);
        return 1.0;
    }
}

Real computePdfAsDynamicPath(const PathVertex &vGhost,
                                        const PathVertex &vLight, const PathVertex &vSensor,
                                        const PathVertex &sampledVLight, const PathVertex &sampledVSensor,
                                        const std::vector<PathVertex> &lightPath, const std::vector<PathVertex> &sensorPath,
                                        int S, int T, const DSubPathState &lightSubpathState, const DSubPathState &sensorSubpathState,
                                        const SamplingTechnique &technique, const Scene &scene, int time)
{
    if (technique == SamplingTechnique::GHOST_FROM_LIGHT)
    {
        const DSubPathState &subpathState = sensorSubpathState;
        Real pdfSum = 0.0;
        //account for the pdfs for dynamic portal sampling strategies
        for(int t=2; t<=T; ++t){ 
            if(is_vertex_dynamic(sensorPath[t-1], scene))
            {
                const PathVertex &dynamicHitVertex = sensorPath[t - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sensorPath[t - 2].position - dynamicHitVertex.position), scene, time);
                // the light and sensor subpath part
                Real pdfLightPath = safe_divide(subpathState.pdf_pairs[1].second, subpathState.pdf_pairs[t - 1].second);
                pdfLightPath *= sampledVLight.transport.pdfOriginA; // still sample the light source directly. //TODO if MIS with NEE too
                Real pdfSensorPath = safe_divide(subpathState.pdf_pairs[T - 1].first, subpathState.pdf_pairs[t - 1].first);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                //TODO check here for some potential modification
                if (t>2 && t < T) //if inner vertex, need to modify the pdf //TODO check here
                { 
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, sensorPath[t - 1], sensorPath[t - 2]), sensorPath[t - 2].transport.pdfRev);
                }
                pdfSum += MIS(pdfThisHit);
            }
        }
        return pdfSum;
    }
    else if (technique == SamplingTechnique::GHOST_FROM_SENSOR)
    {
        const DSubPathState &subpathState = lightSubpathState;
        Real pdfSum = 0.0;
        for (int s = 2; s <= S; ++s)
        {
            if (is_vertex_dynamic(lightPath[s - 1], scene))
            {
                const PathVertex &dynamicHitVertex = lightPath[s - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                if (s < S)
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(lightPath[s].position - dynamicHitVertex.position), scene, time);
                else
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sampledVLight.position - dynamicHitVertex.position), scene, time);
                // the light and sensor subpath parts
                Real pdfLightPath = safe_divide(subpathState.pdf_pairs[S - 1].first, subpathState.pdf_pairs[s - 1].first);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(subpathState.pdf_pairs[1].second, subpathState.pdf_pairs[s - 1].second);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // to match the 3rd sampling strategy, the first vertex along the light source is sampled using pdfDirThisHit, instead of bsdf sampling...
                //TODO check here for some potential modification
                if (s < S)
                {
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, lightPath[s - 1], lightPath[s]), lightPath[s].transport.pdfFwd);
                }
                pdfSum += MIS(pdfThisHit);
            }
        }
        return pdfSum;
    }else if (technique == SamplingTechnique::GHOST_TWO_ENDS)
    {
        Real pdfSum = 0.0;
        for (int s = 2; s <= S; ++s)
        {
            if (is_vertex_dynamic(lightPath[s - 1], scene))
            {
                const PathVertex &dynamicHitVertex = lightPath[s - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                if (s < S)
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(lightPath[s].position - dynamicHitVertex.position), scene, time);
                else
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sampledVLight.position - dynamicHitVertex.position), scene, time);

                Real pdfLightPath = safe_divide(lightSubpathState.pdf_pairs[S - 1].first, lightSubpathState.pdf_pairs[s - 1].first);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(lightSubpathState.pdf_pairs[1].second, lightSubpathState.pdf_pairs[s - 1].second);
                pdfSensorPath *= sensorSubpathState.pdf_pairs[T - 1].first * sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                if (s < S)
                { // sensorPath is longer, only need to modify one nextDynamicSEdge; [] normal obj; <> dynamic object; <o> current sampling starting point
                    ///// eye -- [] -- <o> -- [] -- <> -- [] -- light   //currHitPoint
                    ///// eye -- [] -- <> -- [] -- <o> -- [] -- light   //vertualHitPoint
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, lightPath[s - 1], lightPath[s]), lightPath[s].transport.pdfFwd);
                }
                {
                    Vector3 dir_in = normalize((S>=3? lightPath[2]:sampledVLight).position - sensorPath[0].position);
                    Vector3 dir_out = normalize(sensorPath[1].position - sensorPath[0].position);
                    Real pdfWNextV2Sensor = pdf_sample_bsdf(scene.materials[sensorPath[0].material_id], dir_in, dir_out, sensorPath[0], scene.texture_pool, TransportDirection::TO_VIEW);
                    Real pdfANextV2Sensor = pdf_convert_w_to_A(pdfWNextV2Sensor, sensorPath[0], sensorPath[1]);
                    pdfThisHit *= pdfANextV2Sensor;
                }
                pdfSum += MIS(pdfThisHit);
            }
        }
        for (int t = 2; t <= T; ++t)
        {
            if (is_vertex_dynamic(sensorPath[t - 1], scene))
            {
                const PathVertex &dynamicHitVertex = sensorPath[t - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sensorPath[t - 2].position - dynamicHitVertex.position), scene, time);
                // the light and sensor subpath part
                Real pdfSensorPath = safe_divide(sensorSubpathState.pdf_pairs[T - 1].first, sensorSubpathState.pdf_pairs[t - 1].first);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfLightPath = safe_divide(sensorSubpathState.pdf_pairs[1].second, sensorSubpathState.pdf_pairs[t - 1].second);
                pdfLightPath *= sampledVLight.transport.pdfOriginA; // still sample the light source directly. //TODO if MIS with NEE too
                pdfLightPath *= lightSubpathState.pdf_pairs[S - 1].first;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // TODO modify the edge2light next to the vDynamic and match the corresponding sampling strategy
                if(t==2 && T==2){ 
                    // when s==2 and T==2, to make up the pdf of the GhostEdge, it is the pdfA to generate sensorPath[0] from sensorPath[1]
                    Vector3 dir_in = normalize((T>=3? sensorPath[2]:sampledVSensor).position - sensorPath[1].position);
                    Vector3 dir_out = normalize(sensorPath[0].position - sensorPath[1].position);
                    Real pdfWNextV2Light = pdf_sample_bsdf(scene.materials[sensorPath[1].material_id], dir_in, dir_out, sensorPath[1], scene.texture_pool, TransportDirection::TO_LIGHT);
                    Real pdfANextV2Light = pdf_convert_w_to_A(pdfWNextV2Light, sensorPath[1], sensorPath[0]);
                    pdfThisHit *= pdfANextV2Light;
                }else if(t==2 && t<T){
                    //when s==2 and T>2, to make up the pdf of the GhostEdge, it is the pdfDir on the dynamic point to generate sensorPath[0] 
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, sensorPath[t - 1], sensorPath[t - 2]), sensorPath[t - 2].transport.pdfRev);
                }else{
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, sensorPath[t - 1], sensorPath[t - 2]), sensorPath[t - 2].transport.pdfRev);
                    Vector3 dir_in = normalize((T>=3? sensorPath[2]:sampledVSensor).position - sensorPath[1].position);
                    Vector3 dir_out = normalize(sensorPath[0].position - sensorPath[1].position);
                    Real pdfWNextV2Light = pdf_sample_bsdf(scene.materials[sensorPath[1].material_id], dir_in, dir_out, sensorPath[1], scene.texture_pool, TransportDirection::TO_LIGHT);
                    Real pdfANextV2Light = pdf_convert_w_to_A(pdfWNextV2Light, sensorPath[1], sensorPath[0]);
                    pdfThisHit *= pdfANextV2Light;
                }
                pdfSum += MIS(pdfThisHit);
            }
        }
        return pdfSum;
    }
    else if (technique == UNI_PT)
    {
        Real pdfSum = 0.0;
        for (int s = 1; s <= S; ++s) 
        {
            if (is_vertex_dynamic(lightPath[s - 1], scene))
            {
                const PathVertex &dynamicHitVertex = lightPath[s - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                if (s < S)
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(lightPath[s].position - dynamicHitVertex.position), scene, time);
                else
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sampledVLight.position - dynamicHitVertex.position), scene, time);
                // the light and sensor subpath parts
                Real pdfLightPath = safe_divide(lightSubpathState.pdf_pairs[S - 1].first, lightSubpathState.pdf_pairs[s - 1].first);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(lightSubpathState.pdf_pairs[0].second, lightSubpathState.pdf_pairs[s - 1].second);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                if (s > 1 && s < S)//make the right pdf towards light the pdfDir
                {
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, lightPath[s - 1], lightPath[s]), lightPath[s].transport.pdfFwd);
                }
                else if (s == S && s >= 2)
                { // if the last one, need to change the pdfRef of lightPath[s-2]
                    Real pdfRevSecondLast = pdf_sample_bsdf(scene.materials[dynamicHitVertex.material_id],
                                                            normalize(sampledVLight.position - dynamicHitVertex.position),
                                                            normalize(lightPath[s - 2].position - dynamicHitVertex.position),
                                                            dynamicHitVertex, scene.texture_pool);
                    pdfRevSecondLast = pdf_convert_w_to_A(pdfRevSecondLast, dynamicHitVertex, lightPath[s - 2]);
                    pdfThisHit *= safe_divide(pdfRevSecondLast, lightPath[s - 2].transport.pdfRev);
                }
                Real pdfAPrimaryHit = s == 1 ? pdfAreaThisHit : lightPath[0].transport.pdfRev;
                pdfSum += MIS(pdfThisHit) ;
            }
        }
        return pdfSum;
    }else{
        assert(false);
        return 1.0;
    }
}

Real computeMixMISWeight_better_pt(const PathVertex &vLight, const PathVertex &vSensor,
                                   const PathVertex &sampledVLight, const PathVertex &sampledVSensor,
                                   const std::vector<PathVertex> &lightPath, const std::vector<PathVertex> &sensorPath,
                                   int S, int T, const DSubPathState &lightSubpathState, const DSubPathState &sensorSubpathState,
                                   const SamplingTechnique &technique, const Scene &scene, int time,
                                Real& oCurrPdf)
{

        assert(technique == SamplingTechnique::UNI_PT);
        Real pdfPT = sampledVSensor.transport.pdfOriginA * lightSubpathState.pdf_pairs[S - 1].first * sampledVLight.transport.pdfOriginA;
        Real pdfSum = MIS(pdfPT);
        oCurrPdf = pdfPT;
/////////////////////////////////// if the original path is generated by sampling dynamic////////////////////////////
#ifdef MIS_MIX //TODO remember to toggle here if separately testing the dynamic and ghost sampling 
        pdfSum += computePdfAsDynamicPath({},vLight,vSensor,sampledVLight,sampledVSensor,
                                        lightPath,sensorPath, S, T, lightSubpathState, sensorSubpathState,
                                        technique, scene, time); //no vGhost here just use {}
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////// if the original path is generated by sampling ghost////////////////////////////
        pdfSum  += computePdfAsGhostPath(vLight,vSensor,sampledVLight,sampledVSensor,
                                        lightPath,sensorPath, S, T, lightSubpathState, sensorSubpathState,
                                        technique, scene, time);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
        Real weight = MIS(pdfPT) / pdfSum;
        weight = std::clamp(weight, Real(0), Real(1));
        return weight;
}

/**
 * @brief combine the solid dynamic and ghost dynamic together; account the MIS for each other.
 * cover all the starting sampling strategy to be dynamic portal sampling
 */
Real computeMixMISWeight_better_dynamic(const PathVertex &vLight, const PathVertex &vSensor,
                                   const PathVertex &sampledVLight, const PathVertex &sampledVSensor,
                                   const std::vector<PathVertex> &lightPath, const std::vector<PathVertex> &sensorPath,
                                   int S, int T, const DSubPathState &lightSubpathState, const DSubPathState &sensorSubpathState,
                                   const SamplingTechnique &technique, const Scene &scene, int time,
                                   Real& oCurrPdf)
{
    if (technique == SamplingTechnique::DYNAMIC_FROM_LIGHT)
    {
        const DSubPathState &subpathState = sensorSubpathState;
        // TODO check the 3-vertex path case
        //  pdfCurrHit: pdf of sampling this path by starting from this dynamic object and go towards two directions.
        Real pdfCurrHit = lightPath[0].transport.pdfOriginA; // 1) pdf of sampling this dynamic object
        pdfCurrHit *= sampledVLight.transport.pdfOriginA;    // 2) pdf of the light subpath: directly sample the light source
        pdfCurrHit *= subpathState.pdf_pairs[T - 1].first;   // 3) pdf of the sensor subpath: accumulated in the pdfFwd.
        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        //  need to loop over all the other possibilities
        oCurrPdf = pdfCurrHit;
        
        if(pdfCurrHit==0.) return 0.;
        
        assert(pdfCurrHit == pdfCurrHit);
        Real pdfSum = MIS(pdfCurrHit);
        for (int t = 2; t <= T; ++t) // these are handled by SamplingTechnique::DYNAMIC_TWO_ENDS; so MIS with it
        {
            if (is_vertex_dynamic(sensorPath[t - 1], scene))
            {
                const PathVertex &dynamicHitVertex = sensorPath[t - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sensorPath[t - 2].position - dynamicHitVertex.position), scene, time);
                // the light and sensor subpath part
                Real pdfLightPath = safe_divide(subpathState.pdf_pairs[0].second, subpathState.pdf_pairs[t - 1].second);
                pdfLightPath *= sampledVLight.transport.pdfOriginA; // still sample the light source directly. //TODO if MIS with NEE too
                Real pdfSensorPath = safe_divide(subpathState.pdf_pairs[T - 1].first, subpathState.pdf_pairs[t - 1].first);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // to match the 3rd sampling strategy, the first vertex along the light source is sampled using pdfDirThisHit, instead of bsdf sampling...
                if (t < T)
                { // if t== T then this case is sampled by SamplingTechnique::DYNAMIC_FROM_SENSOR; no need to change it.
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, sensorPath[t - 1], sensorPath[t - 2]), sensorPath[t - 2].transport.pdfRev);
                }else{
                    Real pdfWRev_duetoEnd = pdf_bsdf_at_vertex(sensorPath[t-1], sensorPath[t-2], sampledVSensor,TransportDirection::TO_LIGHT, scene);
                    Real pdfARev_duetoEnd = pdf_convert_w_to_A(pdfWRev_duetoEnd, sensorPath[t-1], sensorPath[t-2]);
                    pdfThisHit *= safe_divide(pdfARev_duetoEnd, sensorPath[t - 2].transport.pdfRev);
                }
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum == pdfSum);
            }
        }
#ifdef USE_PT_AS_4
        // add in the uni-directional pt sampling pdf
        Real pdfPT = safe_divide(subpathState.pdf_pairs[0].second, subpathState.pdf_pairs[T - 1].second); // TODO times the pdf of sampling the camera
        // pdfRev at pdf_pairs[T-2] changed due to sampledVsensor =>[T-1] =>[T-2]
        if (T >= 2 && pdfPT > 0)
        {
            Real pdfWRevTmp = pdf_bsdf_at_vertex(sensorPath[T - 1], sensorPath[T - 2], sampledVSensor, TransportDirection::TO_LIGHT, scene);
            Real pdfARevTmp = pdf_convert_w_to_A(pdfWRevTmp, sensorPath[T - 1], sensorPath[T - 2]);
            pdfPT *= pdfARevTmp / sensorPath[T - 2].transport.pdfRev;
        }
        pdfPT *= sampledVSensor.transport.pdfFwd; // the pdf to sample the first vertex from image/sensor sampling; ==cameraPdfA
        pdfPT *= sampledVLight.transport.pdfOriginA;
        pdfSum += MIS(pdfPT);
        assert(pdfPT == pdfPT);
#endif
#ifdef MIS_MIX
        pdfSum += computePdfAsGhostPath(vLight, vSensor,sampledVLight, sampledVSensor,
                    lightPath, sensorPath, S, T, lightSubpathState, sensorSubpathState,
                    technique, scene, time);
#endif
        assert(pdfSum == pdfSum);
        if(pdfSum==0.) return 0.;
        return MIS(pdfCurrHit) / pdfSum;
    }
    else if (technique == SamplingTechnique::DYNAMIC_FROM_SENSOR)
    {
        const DSubPathState &subpathState = lightSubpathState;
        // pdfCurrHit: pdf of sampling this path by starting from this dynamic object
        Real pdfCurrHit = sensorPath[0].transport.pdfOriginA; // 1) pdf of sampling this dynamic object
        pdfCurrHit *= subpathState.pdf_pairs[S - 1].first;
        pdfCurrHit *= sampledVLight.transport.pdfOriginA; // to sample the lightsource. //we can also omit this because this part is the same

        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        oCurrPdf = pdfCurrHit;
        if(pdfCurrHit==0.)return 0.;
        assert(pdfCurrHit == pdfCurrHit);

        Real pdfSum = MIS(pdfCurrHit);
        for (int s = 2; s <= S; ++s)
        {
            if (is_vertex_dynamic(lightPath[s - 1], scene))
            {
                const PathVertex &dynamicHitVertex = lightPath[s - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                if (s < S)
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(lightPath[s].position - dynamicHitVertex.position), scene, time);
                else
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sampledVLight.position - dynamicHitVertex.position), scene, time);

                // the light and sensor subpath parts
                Real pdfLightPath = safe_divide(subpathState.pdf_pairs[S - 1].first, subpathState.pdf_pairs[s - 1].first);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(subpathState.pdf_pairs[0].second, subpathState.pdf_pairs[s - 1].second);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // to match the 3rd sampling strategy, the first vertex along the light source is sampled using pdfDirThisHit, instead of bsdf sampling...
                if (s < S)
                {
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, lightPath[s - 1], lightPath[s]), lightPath[s].transport.pdfFwd);
                }else{
                    Real pdfWRev_duetoEnd = pdf_bsdf_at_vertex(lightPath[s-1], lightPath[s-2], sampledVLight, TransportDirection::TO_VIEW, scene);
                    Real pdfARev_duetoEnd = pdf_convert_w_to_A(pdfWRev_duetoEnd, lightPath[s-1], lightPath[s-2]);
                    pdfThisHit *= safe_divide(pdfARev_duetoEnd, lightPath[s - 2].transport.pdfRev);
                }
                pdfSum += MIS(pdfThisHit);
            }
        }
#ifdef USE_PT_AS_4
        // add in the uni-directional pt sampling pdf
        Real pdfPT = subpathState.pdf_pairs[S - 1].first;
        pdfPT *= sampledVSensor.transport.pdfFwd;
        pdfPT *= sampledVLight.transport.pdfOriginA;
        pdfSum += MIS(pdfPT);
#endif
#ifdef MIS_MIX
        pdfSum += computePdfAsGhostPath(vLight, vSensor,sampledVLight, sampledVSensor,
                    lightPath, sensorPath, S, T, lightSubpathState, sensorSubpathState,
                    technique, scene, time);
#endif
        return MIS(pdfCurrHit) / pdfSum;
    }
    else if (technique == SamplingTechnique::DYNAMIC_TWO_ENDS)
    {
        Real pdfCurrHit = sensorPath[0].transport.pdfOriginA;
        pdfCurrHit *= sampledVLight.transport.pdfOriginA;
        pdfCurrHit *= lightSubpathState.pdf_pairs[S - 1].first;
        pdfCurrHit *= sensorSubpathState.pdf_pairs[T - 1].first;
        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        oCurrPdf = pdfCurrHit;
        assert(pdfCurrHit == pdfCurrHit);
        if(pdfCurrHit==0.) return 0.;
        Real pdfSum = MIS(pdfCurrHit);
        for (int s = 2; s <= S; ++s)
        {
            if (is_vertex_dynamic(lightPath[s - 1], scene))
            {
                const PathVertex &dynamicHitVertex = lightPath[s - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                if (s < S)
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(lightPath[s].position - dynamicHitVertex.position), scene, time);
                else
                    pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sampledVLight.position - dynamicHitVertex.position), scene, time);

                // the light subpath
                Real pdfLightPath = safe_divide(lightSubpathState.pdf_pairs[S - 1].first, lightSubpathState.pdf_pairs[s - 1].first);
                pdfLightPath *= sampledVLight.transport.pdfOriginA;
                Real pdfSensorPath = safe_divide(lightSubpathState.pdf_pairs[0].second, lightSubpathState.pdf_pairs[s - 1].second);
                pdfSensorPath *= sensorSubpathState.pdf_pairs[T - 1].first * sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                if (s < S)
                { // sensorPath is longer, only need to modify one nextDynamicSEdge; [] normal obj; <> dynamic object; <o> current sampling starting point
                    ///// eye -- [] -- <o> -- [] -- <> -- [] -- light   //currHitPoint
                    ///// eye -- [] -- <> -- [] -- <o> -- [] -- light   //vertualHitPoint
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, lightPath[s - 1], lightPath[s]), lightPath[s].transport.pdfFwd);
                }
                pdfSum += MIS(pdfThisHit);
            }
        }
        for (int t = 2; t <= T; ++t)
        {
            if (is_vertex_dynamic(sensorPath[t - 1], scene))
            {
                const PathVertex &dynamicHitVertex = sensorPath[t - 1];
                Real pdfAreaThisHit(1.0), pdfDirThisHit(1.0);
                pdf_directly_sample_dynamic(pdfAreaThisHit, pdfDirThisHit, dynamicHitVertex, normalize(sensorPath[t - 2].position - dynamicHitVertex.position), scene, time);
                // the light and sensor subpath part
                Real pdfSensorPath = safe_divide(sensorSubpathState.pdf_pairs[T - 1].first, sensorSubpathState.pdf_pairs[t - 1].first);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfLightPath = safe_divide(sensorSubpathState.pdf_pairs[0].second, sensorSubpathState.pdf_pairs[t - 1].second);
                pdfLightPath *= sampledVLight.transport.pdfOriginA; // still sample the light source directly. //TODO if MIS with NEE too
                pdfLightPath *= lightSubpathState.pdf_pairs[S - 1].first;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                // TODO modify the edge2light next to the vDynamic and match the corresponding sampling strategy
                if (t < T)
                { // lightPath is loner, need to modify two nextDynamicSEdge
                    ///// eye -- [] -- <> -- [] -- <o> -- [] -- light   //currHitPoint
                    ///// eye -- [] -- <o> -- [] -- <> -- [] -- light   //virtualHitPoint
                    pdfThisHit *= safe_divide(pdf_convert_w_to_A(pdfDirThisHit, sensorPath[t - 1], sensorPath[t - 2]), sensorPath[t - 2].transport.pdfRev);
                    Real pdfBSDF_nextDynamicSEdge = pdf_bsdf_at_vertex(lightPath[0], lightPath[1], sensorPath[1], TransportDirection::TO_LIGHT, scene); // the ratio should be in area measure, but the G is the same, so skip the G
                    pdfThisHit *= safe_divide(pdfBSDF_nextDynamicSEdge, lightPath[0].transport.pdfDir);
                }
                else
                {                                                                                                                                       // t==T; only the edge on the light side needs to be modified.
                    Real pdfBSDF_nextDynamicSEdge = pdf_bsdf_at_vertex(lightPath[0], lightPath[1], sensorPath[1], TransportDirection::TO_LIGHT, scene); // the ratio should be in area measure, but the G is the same, so skip the G
                    pdfThisHit *= safe_divide(pdfBSDF_nextDynamicSEdge, lightPath[0].transport.pdfDir);
                }
                pdfSum += MIS(pdfThisHit);
            }
        }

#ifdef USE_PT_AS_4
        Real pdfPT = lightSubpathState.pdf_pairs[S - 1].first * safe_divide(sensorSubpathState.pdf_pairs[0].second, sensorSubpathState.pdf_pairs[T - 1].second);
        if (T >= 2 && pdfPT > 0)
        {
            Real pdfWRevTmp = pdf_bsdf_at_vertex(sensorPath[T - 1], sensorPath[T - 2], sampledVSensor, TransportDirection::TO_LIGHT, scene);
            Real pdfARevTmp = pdf_convert_w_to_A(pdfWRevTmp, sensorPath[T - 1], sensorPath[T - 2]);
            pdfPT *= pdfARevTmp / sensorPath[T - 2].transport.pdfRev;
        }
        pdfPT *= sampledVSensor.transport.pdfFwd;
        pdfPT *= sampledVLight.transport.pdfOriginA;
        // modify the edge2light next to the vDynamic
        {
            Real pdfBSDF_nextDynamicSEdge = pdf_bsdf_at_vertex(lightPath[0], lightPath[1], sensorPath[1], TransportDirection::TO_LIGHT, scene); // the ratio should be in area measure, but the G is the same, so skip the G
            pdfPT *= safe_divide(pdfBSDF_nextDynamicSEdge, lightPath[0].transport.pdfDir);
        }
        pdfSum += MIS(pdfPT);
#endif
#ifdef MIS_MIX
        Real pdfGhost = computePdfAsGhostPath(vLight, vSensor,sampledVLight, sampledVSensor,
                    lightPath, sensorPath, S, T, lightSubpathState, sensorSubpathState,
                    technique, scene, time);
        if(pdfGhost > 0)
            pdfSum += MIS(pdfGhost);
#endif
        return MIS(pdfCurrHit) / pdfSum;
    }
    else if (technique == UNI_PT)
    {//TODO attention , no pdfLight for all here.....
        assert(false);
        return computeMixMISWeight_better_pt(vLight, vSensor, sampledVLight, sampledVSensor,
                                         lightPath, sensorPath, S, T, lightSubpathState, sensorSubpathState,
                                         technique, scene, time, oCurrPdf);
    }
    else
    {
        assert(false);
        int hit_dynamic_times = vLight.hit_dynamic_times + vSensor.hit_dynamic_times - 1;
        Real weight = hit_dynamic_times > 0 ? 1 / (Real)hit_dynamic_times : 0.0;
        return weight;
    }
}

/**
 * @brief combine the solid dynamic and ghost dynamic together; account the MIS for each other.
 * cover all the starting sampling strategy to be ghost portal sampling
 */
Real computeMixMISWeight_better_ghost(const PathVertex &vGhost,
                                        const PathVertex &vLight, const PathVertex &vSensor,
                                        const PathVertex &sampledVLight, const PathVertex &sampledVSensor,
                                        const std::vector<PathVertex> &lightPath, const std::vector<PathVertex> &sensorPath,
                                        int S, int T, const DSubPathState &lightSubpathState, const DSubPathState &sensorSubpathState,
                                        const SamplingTechnique &technique, const Scene &scene, int time,
                                        Real& oCurrPdf)
{
    if (technique == SamplingTechnique::GHOST_FROM_LIGHT)
    {
        const DSubPathState &subpathState = sensorSubpathState;
        assert(lightPath.empty());
        // pdfCurrHit: pdf of sampling this path by starting from this ghost object
        Real pdfANextGhost = vGhost.transport.pdfOriginA * pdfGhostToSolidNext(vGhost.position, vGhost.geometric_normal, sampledVLight, sensorPath[1]);
        Real pdfCurrHit = pdfANextGhost;
        pdfCurrHit *= sampledVLight.transport.pdfOriginA;
        pdfCurrHit *= subpathState.pdf_pairs[T - 1].first;
        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        oCurrPdf = pdfCurrHit;
        Real pdfSum(0.);             // this will be included in the ghost edge below (the first edge)
        for (int t = 2; t <= T; ++t) // the first one t==1 is light source
        {                            // TODO the last one can just be a envmap void vertex
            std::vector<int> ghost_shape_ids = sensorPath[t - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue; // this segment never hits a ghost obj
            // this segment passes the ghost object //or that sensorPath[t-1].ghost_ids.size()>0
            
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
                { // pdfA of sampling the first solid vertex by directly connecting from the light source to the ghost
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
            }
        }
#ifdef USE_PT_AS_4_GHOST
        Real pdfPT = safe_divide(subpathState.pdf_pairs[1].second, subpathState.pdf_pairs[T-1].second);
        if (T >= 3 && pdfPT > 0)
        {// the pdfRev at sensorPath[T-2] is possibly changed because we use the sampledVsensor instead the other vertex along the way
            Real pdfWRevTmp = pdf_bsdf_at_vertex(sensorPath[T - 1], sensorPath[T - 2], sampledVSensor, TransportDirection::TO_LIGHT, scene);
            Real pdfARevTmp = pdf_convert_w_to_A(pdfWRevTmp, sensorPath[T - 1], sensorPath[T - 2]);
            pdfPT *= pdfARevTmp / sensorPath[T - 2].transport.pdfRev;
        }
        pdfPT *= sampledVSensor.transport.pdfOriginA; // the pdf to sample the first vertex from image/sensor sampling; ==cameraPdfA
        pdfPT *= sampledVLight.transport.pdfOriginA;
        pdfSum += MIS(pdfPT);
        assert(pdfPT == pdfPT);
#endif
#ifdef MIS_MIX
        Real pdfDynamicPaths = computePdfAsDynamicPath(vGhost,vLight,vSensor,sampledVLight,sampledVSensor,
                                        lightPath,sensorPath, S, T, lightSubpathState, sensorSubpathState,
                                        technique, scene, time);
        pdfSum += pdfDynamicPaths;
#endif
        if(pdfSum ==0.0) return 0.0;
        Real weight = MIS(pdfCurrHit) / pdfSum;
        weight = std::clamp(weight, Real(0), Real(1));
        assert( weight <= 1.0);
        return weight;
    }
    else if (technique == SamplingTechnique::GHOST_FROM_SENSOR)
    {
        const DSubPathState &subpathState = lightSubpathState;
        Real pdfCurrHit = vGhost.transport.pdfOriginA * pdfGhostToSolidNext(vGhost.position, vGhost.geometric_normal, sampledVSensor, lightPath[1]);
        pdfCurrHit *= subpathState.pdf_pairs[S - 1].first;
        pdfCurrHit *= sampledVLight.transport.pdfOriginA;
        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        oCurrPdf = pdfCurrHit;
        Real pdfSum(0.);
        for (int s = 2; s <= S; ++s)
        {
            std::vector<int> ghost_shape_ids = lightPath[s - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue;
            
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
                assert(pdfSum > 0);
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
#ifdef MIS_MIX
        pdfSum += computePdfAsDynamicPath(vGhost,vLight,vSensor,sampledVLight,sampledVSensor,
                                        lightPath,sensorPath, S, T, lightSubpathState, sensorSubpathState,
                                        technique, scene, time);
#endif
        if(pdfSum ==0.0) return 0.0;
        Real weight = MIS(pdfCurrHit) / pdfSum;
        weight = std::clamp(weight, Real(0), Real(1)); //these are just to handle numerical issues
        return weight;
    }
    else if (technique == SamplingTechnique::GHOST_TWO_ENDS)
    {
        Real pdfCurrHit = vGhost.transport.pdfOriginA;
        pdfCurrHit *= sampledVLight.transport.pdfOriginA;
        pdfCurrHit *= sampledVSensor.transport.pdfDir;
        pdfCurrHit *= lightSubpathState.pdf_pairs[S - 1].first;
        pdfCurrHit *= sensorSubpathState.pdf_pairs[T - 1].first;
        // pdf to sample the first segment dir etc..
        Real pdfGhostEdge = pdfInnerGhostEdge(vGhost.position, vGhost.geometric_normal, vGhost.transport.pdfDir, sensorPath[0], sensorPath[1]);
        pdfCurrHit *= pdfGhostEdge;
        oCurrPdf = pdfCurrHit;
        Real pdfSum(0.);
        for (int s = 3; s <= S; ++s) // here starts from 3
        {
            std::vector<int> ghost_shape_ids = lightPath[s - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue;
            
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
                    Vector3 dir_in = normalize(lightPath[2].position - sensorPath[0].position);
                    Vector3 dir_out = normalize(sensorPath[1].position - sensorPath[0].position);
                    Real pdfWNextV2Sensor = pdf_sample_bsdf(scene.materials[sensorPath[0].material_id], dir_in, dir_out, sensorPath[0], scene.texture_pool, TransportDirection::TO_VIEW);
                    Real pdfANextV2Sensor = pdf_convert_w_to_A(pdfWNextV2Sensor, sensorPath[0], sensorPath[1]);
                    pdfThisHit *= pdfANextV2Sensor;
                } 
                // pdf to sample the first vertex on lightPath
                Real pdfGhostEdge = pdfInnerGhostEdge(ghostPos, ghostNormal, pdfDirThisHit, lightPath[s - 1], lightPath[s - 2]);
                pdfThisHit *= pdfGhostEdge;
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum == pdfSum);
            }
        }
        for (int t = 2; t <= T; ++t) // the first one t==1 is light source// here starts from 2
        {
            std::vector<int> ghost_shape_ids = sensorPath[t - 1].ghost_ids;
            int num_ghost_hits = ghost_shape_ids.size();
            if (num_ghost_hits <= 0)
                continue; // this segment never hits a ghost obj
            // this segment passes the ghost object //or that sensorPath[t-1].ghost_ids.size()>0
            
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
                assert(pdfSum == pdfSum);
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
                assert(pdfSum == pdfSum);
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
                pdfThisHit*= lightSubpathState.pdf_pairs[S - 1].first;
                pdfThisHit *= sampledVLight.transport.pdfOriginA * sampledVSensor.transport.pdfDir;
                pdfSum += MIS(pdfThisHit);
                assert(pdfSum == pdfSum);
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
        assert(pdfSum == pdfSum);
#endif
#ifdef MIS_MIX
        pdfSum += computePdfAsDynamicPath(vGhost,vLight,vSensor,sampledVLight,sampledVSensor,
                                        lightPath,sensorPath, S, T, lightSubpathState, sensorSubpathState,
                                        technique, scene, time);
#endif
        if(pdfSum == 0.0) return 0.0;
        assert(pdfSum == pdfSum);
        assert(pdfCurrHit == pdfCurrHit);
        Real weight = MIS(pdfCurrHit) / pdfSum;
        weight = std::clamp(weight, Real(0), Real(1));
        return weight;
    }
    else if (technique == UNI_PT)
    {
        return computeMixMISWeight_better_pt(vLight, vSensor, sampledVLight, sampledVSensor,
                                    lightPath, sensorPath, S, T, lightSubpathState, sensorSubpathState,
                                    technique, scene, time, oCurrPdf);
    }
    else
    {
        assert(false);
        return 1.0;
    }
}
