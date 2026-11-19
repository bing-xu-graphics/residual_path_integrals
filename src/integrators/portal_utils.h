#pragma once

#include "scene.h"
#include "pcg.h"
#include "../utils.h"
/**
 * @brief sample a point and then a direction on the (solid) dynamic object
 *
 * @param vertex output the sampled vertex; filled it with all the neccessary information needed as a shading point
 * @param diro output the sampled outgoing dir
 * @param pdfDir output the pdf of outgoing dir
 * @param rng
 * @param scene
 */

void sampleDynamic(PathVertex &vertex,
                   Vector3 &diro,
                   //    Real &pdfDir,
                   pcg32_state &rng,
                   const Scene &scene,
                   int time)
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

    //////////////////
    // sample the outgoing direction
    Real pdfDynamicDir;
    Vector2 dynamic_dir_rdn{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Vector3 dirFromDynamic = sample_dir_from_shape(dynamic, vertex, dynamic_dir_rdn, pdfDynamicDir, time); // sample the hemisphere

    // output
    diro = dirFromDynamic;
    vertex.hit_dynamic_times = 1;
    vertex.transport.pdfOriginA = pdfDynamic;
    vertex.transport.pdfDir = pdfDynamicDir;
    vertex.transport.pdfFwd = 1.0; // TODO change this; this should be the pdf to generate this dynamic vertex along random walk direction
    vertex.transport.pathLength = 1;
    vertex.transport.contrib = make_const_spectrum(1.0 / pdfDynamic); // FIXME attention
    return;
}

/**
 * @brief scattering sampling for the input vertex; currently only supports surface
 *
 * @param vertex
 * @param dir_in
 * @param dir_out
 * @param pdfW
 * @param rng
 * @param scene
 * @return true
 * @return false
 */
bool sampleScattering(PathVertex &vertex,
                      const Vector3 &dir_in,
                      Vector3 &dir_out,
                      Real &pdfW,
                      pcg32_state &rng,
                      const Scene &scene)
{
    int material_id = -1;
    assert(vertex.type == VertexType::Surface);

    assert(vertex.material_id >= 0);
    const Material &mat = scene.materials[vertex.material_id];
    Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
    std::optional<BSDFSampleRecord> bsdf_sample_ =
        sample_bsdf(mat,
                    dir_in,
                    vertex,
                    scene.texture_pool,
                    bsdf_rnd_param_uv,
                    bsdf_rnd_param_w);
    if (!bsdf_sample_)
    {
        // BSDF sampling failed. Abort the loop.
        return false;
    }
    const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
    Vector3 dir_bsdf = bsdf_sample.dir_out;
    vertex.roughness = bsdf_sample.roughness;
    Real pdf_bsdf = pdf_sample_bsdf(
        mat, dir_in, dir_bsdf, vertex, scene.texture_pool);
    if (pdf_bsdf <= 0)
        return false;
    pdfW = pdf_bsdf;
    dir_out = dir_bsdf;
    return true;
}

/**
 * @brief generate two subpaths starting from the middle (solid) dynamic object;
 * 1) sample the surface of the (solid) dynamic object; 2d random number
 * 2) sample a initial direction for one subpath; without loss of generality, use this for dir2Light
 * 3) sample another direction dir2Sensor, conditional on dir2Light using bsdf sampling
 * //TODO this should be optimized by sampling a 4d/6d cached importance-bound-radiance guiding system.
 *
 * @param lightPath
 * @param sensorPath
 * @param rng
 * @param scene
 * @return true
 * @return false
 */
bool genDynamicSubpaths(std::vector<PathVertex> &lightPath,
                        std::vector<PathVertex> &sensorPath,
                        pcg32_state &rng,
                        const Scene &scene,
                        int time)
{
    int maxDepth = scene.options.max_depth;
    // sample a point on the dynamic surface first and then sample a 2d direction
    PathVertex vDynamic;
    Vector3 dir2Light, dir2Sensor;                        // the initial two directions for light and sensor subpaths
    sampleDynamic(vDynamic, dir2Light, rng, scene, time); // without loss of generality
    // sample the other dir using bsdf
    Real pdfDir2Light = vDynamic.transport.pdfDir, pdfDir2Sensor;
    bool isOtherDir = sampleScattering(vDynamic, dir2Light, dir2Sensor, pdfDir2Sensor, rng, scene);
    if (!isOtherDir)
    {
        printf("[WARN]  sampleScattering() failed once!\n");
        return false;
    }
    Spectrum throughputS = make_const_spectrum(1.0) / vDynamic.transport.pdfOriginA; // TODO attention not including  / vDynamic.transport.pdfDir for the throughput
    Spectrum throughputT = make_const_spectrum(1.0) / vDynamic.transport.pdfOriginA;
    vDynamic.transport.contrib = throughputS; // FIXME attention
    vDynamic.hit_dynamic_times = 1;           // coz this one is sampled from dynamic object.
    lightPath.push_back(vDynamic);
    vDynamic.transport.pdfDir = pdfDir2Sensor;
    vDynamic.transport.pdfFwd = pdfDir2Sensor;
    sensorPath.push_back(vDynamic);                                                                                                           // update the pdfDir and copy into the path vector                                                                                                     // TODO attention store the vDynamic in boteh subpaths
    RayDifferential ray_diff;                                                                                                                 // dummy_variable
    int k = maxDepth - 1;                                                                                                                     // numVertices - 2 save one spot for camera, and save one for direct sampling of the light
    randomWalk(lightPath, make_const_spectrum(1.0), pdfDir2Light, dir2Light, ray_diff, k, 1, TransportDirection::TO_LIGHT, rng, scene, time); // TODO check the throughput doens't include the dynamic contribution
    randomWalk(sensorPath, make_const_spectrum(1.0), pdfDir2Sensor, dir2Sensor, ray_diff, k, 1, TransportDirection::TO_VIEW, rng, scene, time);
    // TODO send some info out?? dir2Light, pdfDir2Light stuff
    return true;
}

enum SamplingTechnique
{
    DYNAMIC_FROM_LIGHT = 0,
    DYNAMIC_FROM_SENSOR = 1,
    DYNAMIC_TWO_ENDS = 2,
    UNI_PT = 3,
    GHOST_FROM_LIGHT = 4,
    GHOST_FROM_SENSOR = 5,
    GHOST_TWO_ENDS = 6
};

/**
 * @brief record the information package for the dynamic subpaths or ghost subpath
 *
 */
struct DSubPathState
{
    Vector3 mOrigin;                                 // path origin
    Vector3 mDirection;                              // where to go next
    Spectrum mThroughput = make_const_spectrum(1.0); // the throughput of the path so far
    int mPathLength : 30;                            // number of path segments, including this
    int mIsFiniteLight : 1;                          // just generate by finite light
    int mSpecularPath : 1;                           // All scattering events so far were specular

    // TODO currently redundant with PathVertex
    int hit_dynamic_times = 0;
    int hit_ghost_times = 0;
    int time = 0;                                 // match the ray from which the intersection is generated
    std::vector<std::pair<Real, Real>> pdf_pairs; // the pdf2Light and pdf2Sensor at dynamic intersection points; size == # of intersections with dynamic objects
    DSubPathState(int n_pairs) : pdf_pairs(n_pairs, std::make_pair(1.0, 1.0)) {}
};

// continued product of the pdfFwd and pdfRev along a subpath, pdfFwd is multiplied along the random walk direction
//  while the pdfRev is multiplied along the reverse random walk direction
void cumprod_pdfs(const std::vector<PathVertex> &subPath, DSubPathState &pathState)
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

bool is_vertex_dynamic(const PathVertex &vertex, const Scene &scene)
{
    assert(vertex.shape_id >= 0);
    return is_dynamic(scene.shapes[vertex.shape_id]);
}

// the pdf to directly sample this dynmaic object
void pdf_directly_sample_dynamic(Real &pdfArea, Real &pdfDir,
                                 const PathVertex &vertex, const Vector3 &sampled_dir, const Scene &scene, int time)
{
    assert(is_vertex_dynamic(vertex, scene));
    assert(vertex.dynamic_id >= 0);
    const Dynamic &dyn_agent = scene.dynamics[vertex.dynamic_id];
    const Shape &dynamic = (time == NEW_TIME) ? dyn_agent.base_shape : dyn_agent.old_shape;
    PointAndNormal point_on_dynamic{vertex.position, vertex.geometric_normal};
    pdfArea = dynamic_pmf(scene, vertex.dynamic_id) *
              pdf_point_on_shape(dynamic, point_on_dynamic, Vector3{}, time);
    assert(pdfArea != INFINITY);
    pdfDir = pdf_dir_from_shape(dynamic, vertex, sampled_dir, time);
}

void pdf_directly_sample_ghost(Real &pdfArea, Real &pdfDir, // output
                               const int &shape_id, const PointAndNormal &pointAndNormal,
                               const Scene &scene, int time)
{
    int dynamic_id = get_dynamic_id(scene.ghost_shapes[shape_id - scene.shapes.size()]);
    assert(dynamic_id >= 0);
    const Dynamic &dyn_agent = scene.dynamics[dynamic_id];
    const Shape &ghost = (time == NEW_TIME) ? dyn_agent.old_shape : dyn_agent.base_shape;
    pdfArea = dynamic_pmf(scene, dynamic_id) *
              pdf_point_on_shape(ghost, pointAndNormal, Vector3{}); // TODO: ignored the point_on_shape
    pdfDir = UniformSpherePdf();
}
// convert pdfW from current to sample the toSample vertex into the pdfA of sampling tosample vertex
Real pdf_convert_w_to_A(Real pdfW, const PathVertex &current, const PathVertex &tosample)
{
    Real cosTheta = fabs(dot(tosample.geometric_normal, normalize(tosample.position - current.position))); // the theta at the vertex to be sampled
    Real dist2 = distance_squared(current.position, tosample.position);
    return pdfW * cosTheta / dist2;
}

Real safe_divide(const Real& value1, const Real& value2)
{
    if(value1 == value2) return 1.0;
    if (value2 == 0.0)
        return value1 / 1e-30;
    else
        return value1 / value2;
}

// (S-1) is the index of vLight, (T-1) is the index of vSensor
Real computePortalMISWeight_better(const PathVertex &vLight, const PathVertex &vSensor,
                                   const PathVertex &sampledVLight, const PathVertex &sampledVSensor,
                                   const std::vector<PathVertex> &lightPath, const std::vector<PathVertex> &sensorPath,
                                   int S, int T, const DSubPathState &lightSubpathState, const DSubPathState &sensorSubpathState,
                                   const SamplingTechnique &technique, const Scene &scene, int time)
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
        pdfPT *= sampledVSensor.transport.pdfOriginA; // the pdf to sample the first vertex from image/sensor sampling; ==cameraPdfA
        pdfPT *= sampledVLight.transport.pdfOriginA;
        pdfSum += MIS(pdfPT);
        assert(pdfPT == pdfPT);
#endif
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
                }
                pdfSum += MIS(pdfThisHit);
            }
        }
#ifdef USE_PT_AS_4
        // add in the uni-directional pt sampling pdf
        Real pdfPT = subpathState.pdf_pairs[S - 1].first;
        pdfPT *= sampledVSensor.transport.pdfOriginA;
        pdfPT *= sampledVLight.transport.pdfOriginA;
        pdfSum += MIS(pdfPT);
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
                }//TODO check here!!!!!
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
                }//TODO check here!!!!!
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
        pdfPT *= sampledVSensor.transport.pdfOriginA;
        pdfPT *= sampledVLight.transport.pdfOriginA;
        // modify the edge2light next to the vDynamic
        {
            Real pdfBSDF_nextDynamicSEdge = pdf_bsdf_at_vertex(lightPath[0], lightPath[1], sensorPath[1], TransportDirection::TO_LIGHT, scene); // the ratio should be in area measure, but the G is the same, so skip the G
            pdfPT *= safe_divide(pdfBSDF_nextDynamicSEdge, lightPath[0].transport.pdfDir);
        }
        pdfSum += MIS(pdfPT);
#endif
        return MIS(pdfCurrHit) / pdfSum;
    }
    else if (technique == UNI_PT)
    {//TODO attention , no pdfLight for all here.....

        Real pdfPT = sampledVSensor.transport.pdfOriginA * lightSubpathState.pdf_pairs[S - 1].first;
        Real pdfSum = MIS(pdfPT);

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
                Real pdfSensorPath = safe_divide(lightSubpathState.pdf_pairs[0].second, lightSubpathState.pdf_pairs[s - 1].second);
                pdfSensorPath *= sampledVSensor.transport.pdfDir;
                Real pdfThisHit = pdfAreaThisHit * pdfLightPath * pdfSensorPath;
                if (s > 1 && s < S)
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
        return MIS(pdfPT) / pdfSum;
    }
    else
    {
        assert(false);
        int hit_dynamic_times = vLight.hit_dynamic_times + vSensor.hit_dynamic_times - 1;
        Real weight = hit_dynamic_times > 0 ? 1 / (Real)hit_dynamic_times : 0.0;
        return weight;
    }
}

Real computePortalMISWeight(const PathVertex &vLight, const PathVertex &vSensor)
{
    int hit_dynamic_times = vLight.hit_dynamic_times + vSensor.hit_dynamic_times - 1;
#ifdef USE_PT_AS_4          // using path tracing as the forth sampling strategy
    hit_dynamic_times += 1; // TO remove just for testing
#endif
    Real weight = hit_dynamic_times > 0 ? 1 / (Real)hit_dynamic_times : 0.0;
    return weight;
}

// this will combine the pdfs of solid dynamic objects and ghost objects
Real computeDynamicPathMISWeight()
{
    return 1.0;
}

/// combine the naive mis of solid dynamic objects and ghost objects
Real compute_mis_solid_part(const PathVertex &sampledVLight,
                            const PathVertex &vLight,
                            const PathVertex &vSensor,
                            const PathVertex &sampledVSensor,
                            const SamplingTechnique &technique)
{
    int hit_dynamic_times = vLight.hit_dynamic_times + vSensor.hit_dynamic_times - 1;
    
    int hit_ghost_times;
    if (technique == SamplingTechnique::DYNAMIC_FROM_LIGHT)
    {
        hit_ghost_times = vSensor.hit_ghost_times 
            + sampledVLight.ghost_ids.size() 
            + sampledVSensor.ghost_ids.size();
    }
    else if (technique == SamplingTechnique::DYNAMIC_FROM_SENSOR)
    {
        hit_ghost_times = vLight.hit_ghost_times 
            + sampledVLight.ghost_ids.size() 
            + sampledVSensor.ghost_ids.size(); 
    }
    else if (technique == SamplingTechnique::DYNAMIC_TWO_ENDS){
        hit_ghost_times = vLight.hit_ghost_times + vSensor.hit_ghost_times
                        + sampledVLight.ghost_ids.size() + sampledVSensor.ghost_ids.size();
    }
     
    int total_hit_times = hit_dynamic_times + hit_ghost_times;
    Real weight = total_hit_times > 0 ? 1 / (Real)total_hit_times : 0.0;
    return weight;
}
