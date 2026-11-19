#pragma once

#include "scene.h"
#include "pcg.h"
#include "utils.h"
/**
 * @brief this should be similar in both path tracing and light tracing;
 * the goal is to find all the interactions along the way util bounces == max_depth
 * and record the vertex as well as the pdfFwd, pdfRev
 * Some note: we only store the last light vertex if the direction is towards the light; and in the Light Vertex, 
 * throughput is multipled with the emission and the needed G.
 */
inline void randomWalk(std::vector<PathVertex> &path,
                       Spectrum throughput, // the accumulated throughput * L
                       Real pdf,            // the accumulated forward pdf from the prev vertex
                       const Vector3 &initDir,
                       RayDifferential &ray_diff,
                       int numVertices,
                       int startIndex,
                       const TransportDirection &mode,
                       pcg32_state &rng,
                       const Scene &scene,
                       int time) // TODO remove the default time
{
    Vector3 nextDir(initDir);
    Real eta_scale = Real(1);
    for (int idx = startIndex; idx < numVertices; ++idx)
    {
        if( idx == numVertices - 1 && mode == TransportDirection::TO_VIEW) break; //check this 
        PathVertex &prev = path[idx - 1];
        Ray ray{prev.position, nextDir, get_intersection_epsilon(scene), infinity<Real>(), time};

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        if (!vertex_)
        {
            if (mode == TransportDirection::TO_LIGHT) // if TO_CAMERA then just break
            {
                if (has_envmap(scene)) //TODO
                { // Radiance check if hit envmap

                    const Light &envmap = get_envmap(scene);
                    Spectrum Emission = emission(envmap,
                                        -ray.dir, // pointing outwards from light
                                        ray_diff.spread,
                                        PointAndNormal{}, // dummy parameter for envmap
                                        scene);

                    PathVertex vertex;
                    vertex.type = VertexType::Light;
                    vertex.hit_dynamic_times = prev.hit_dynamic_times;
                    vertex.hit_ghost_times = (prev.type == VertexType::Camera || prev.type == VertexType::Light) ? 0 : prev.hit_ghost_times;
                    vertex.transport.contrib = throughput * Emission;  //TODO whether to include the Emission
                    vertex.ghost_ids = ray.ghost_shape_ids;
                    vertex.ghost_hit_positions = ray.ghost_hit_positions;
                    vertex.ghost_hit_geom_normals = ray.ghost_hit_geom_normals;
                    vertex.roughness = 1.0;//always set the light source to be rough
                    //set the pdfs                                                
                    Real G1 = 1, G2 = 1;//for envmap
                    PointAndNormal light_point{Vector3{0, 0, 0}, -ray.dir}; // pointing outwards from light
                    Real pdfALight = light_pmf(scene, scene.envmap_light_id) *
                            pdf_point_on_light(envmap, light_point, prev.position, scene);
                    vertex.transport.pdfFwd = pdf * G1; //set the pdfFwd
                    vertex.transport.pdfOriginA = pdfALight;//set the pdfA; its supposed to be pdf * G1, here just carry this for mis weight computation
                    Real pdfBsdfW_rev = pdf_dir_from_light(envmap, light_point, -ray.dir, scene);
                    prev.transport.pdfRev = pdfBsdfW_rev * G2; //TODO do i need to set this?
                    path.push_back(vertex);
                }
            }
            break;
        }

        // hit something
        PathVertex vertex = *vertex_;
        vertex.hit_dynamic_times = prev.hit_dynamic_times;
        vertex.hit_ghost_times = path.size() < 2 ? 0 : prev.hit_ghost_times; // never get from the first vertex; due to that we never give the prev a ghost
        vertex.time = time;
        vertex.roughness = 1.0;//TODO remove this
        Vector3 wi = -ray.dir;
        const Material &mat = scene.materials[vertex.material_id];
        //////////////////////////////////////////////////////////
        // TODO if light just break the path
        if (is_light(scene.shapes[vertex.shape_id])) /////////TODO check this.
        {
            if(mode == TransportDirection::TO_VIEW) break;
            vertex.type = VertexType::Light;
            vertex.transport.pathLength = prev.transport.pathLength + 1;
            Spectrum Emission = emission(vertex, -ray.dir, scene);
            
            //set the pdfs
            int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
            assert(light_id >= 0);
            const Light &hitLight = scene.lights[light_id];
            PointAndNormal light_point{vertex.position, vertex.geometric_normal};
            Real pdfALight = light_pmf(scene, light_id) *
                      pdf_point_on_light(hitLight, light_point, prev.position, scene);

            // convert pdfW to pdfA
            Real G1 = fabs(dot(wi, vertex.geometric_normal)) / length_squared(prev.position - vertex.position);
            Real G2 = fabs(dot(wi, prev.geometric_normal)) / length_squared(prev.position - vertex.position);
            vertex.transport.pdfFwd = pdf * G1; //set the pdfFwd
            //pdfRev would be directly sample from the light source. //but do we actually need this pdfRev?
            vertex.transport.pdfOriginA = pdfALight; //set the pdfA; its supposed to be pdf * G1, here just carry this for mis weight computation
            vertex.transport.contrib = throughput * G1 * Emission; // not including the Emission here 

            Real pdfBsdfW_rev = pdf_dir_from_light(hitLight, light_point, -wi, scene);
            prev.transport.pdfRev = pdfBsdfW_rev * G2;  //TODO do i need to set this?

            vertex.ghost_ids = ray.ghost_shape_ids;
            vertex.ghost_hit_positions = ray.ghost_hit_positions;
            vertex.ghost_hit_geom_normals = ray.ghost_hit_geom_normals;

            path.push_back(vertex); 
            break;
        }

        if( idx == numVertices - 1) break;

        if (is_dynamic(scene.shapes[vertex.shape_id]))
        {
            vertex.hit_dynamic_times += 1;
            vertex.dynamic_id = get_dynamic_id(scene.shapes[vertex.shape_id]);
        }

        { // accumulate the num_ghost_hits
            if (prev.type != VertexType::Ghost){ //TODO check prev actually could never be Ghost right?
                vertex.hit_ghost_times += ray.num_ghost_hits;
                vertex.ghost_ids = ray.ghost_shape_ids;
                vertex.ghost_hit_positions = ray.ghost_hit_positions;
                vertex.ghost_hit_geom_normals = ray.ghost_hit_geom_normals;
            }
            ray.num_ghost_hits = 0;
        }
        Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);

        std::optional<BSDFSampleRecord> bsdf_sample_ =
            sample_bsdf(mat,
                        wi,
                        vertex,
                        scene.texture_pool,
                        bsdf_rnd_param_uv,
                        bsdf_rnd_param_w,
                        mode);
        if (!bsdf_sample_) // if no valid sample; would rather not push in the vertex???
            break;
        const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
        Vector3 dirBSDF = bsdf_sample.dir_out;
        // Update ray differentials & eta_scale
        if (bsdf_sample.eta == 0)
        {
            ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
        }
        else
        {
            ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
            eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
        }

        Spectrum f = eval(mat, wi, dirBSDF, vertex, scene.texture_pool); // this includes the cosTheta
        if(f.x<0||f.y<0||f.z<0) break;//don't about the hell of this numerical texture error

        Real pdfBsdfW = pdf_sample_bsdf(mat, wi, dirBSDF, vertex, scene.texture_pool);
        if (pdfBsdfW <= 0)
            break;
        {
            vertex.type = VertexType::Surface;
            vertex.roughness = bsdf_sample.roughness;
            vertex.transport.contrib = throughput;
            vertex.transport.pathLength = prev.transport.pathLength + 1;
            // convert pdfW to pdfA
            Real G1 = fabs(dot(wi, vertex.geometric_normal)) / length_squared(prev.position - vertex.position);
            vertex.transport.pdfFwd = pdf * G1;
            Real pdfBsdfW_rev = pdf_sample_bsdf(mat, dirBSDF, wi, vertex, scene.texture_pool);
            Real G2 = fabs(dot(wi, prev.geometric_normal)) / length_squared(prev.position - vertex.position);
            prev.transport.pdfRev = pdfBsdfW_rev * G2;
            path.push_back(vertex);
        }

        // update information for loop
        pdf = pdfBsdfW; // record to store in next
        throughput *= f / pdfBsdfW;
        nextDir = dirBSDF;
    }
}

inline void randomWalk_easy(std::vector<PathVertex> &path,
                       Spectrum throughput, // the accumulated throughput * L
                       Real pdf,            // the accumulated forward pdf from the prev vertex
                       const Vector3 &initDir,
                       RayDifferential &ray_diff,
                       int numVertices,
                       int startIndex,
                       const TransportDirection &mode,
                       pcg32_state &rng,
                       const Scene &scene,
                       int time) // TODO remove the default time
{
    Vector3 nextDir(initDir);
    Real eta_scale = Real(1);
    for (int idx = startIndex; idx < numVertices; ++idx)
    {
        PathVertex &prev = path[idx - 1];
        Ray ray{prev.position, nextDir, get_intersection_epsilon(scene), infinity<Real>(), time};

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        if (!vertex_)
        {
            break;
        }

        // hit something
        PathVertex vertex = *vertex_;
        vertex.hit_dynamic_times = prev.hit_dynamic_times;
        vertex.hit_ghost_times = path.size() < 2 ? 0 : prev.hit_ghost_times; // never get from the first vertex; due to that we never give the prev a ghost
        vertex.time = time;
        vertex.roughness = 1.0;//TODO remove this
        Vector3 wi = -ray.dir;
        const Material &mat = scene.materials[vertex.material_id];
        //////////////////////////////////////////////////////////
        // TODO if light just break the path
        if (is_light(scene.shapes[vertex.shape_id])) /////////TODO check this.
        {
            break;
        }

        if (is_dynamic(scene.shapes[vertex.shape_id]))
        {
            vertex.hit_dynamic_times += 1;
            vertex.dynamic_id = get_dynamic_id(scene.shapes[vertex.shape_id]);
        }

        { // accumulate the num_ghost_hits
            if (prev.type != VertexType::Ghost){ //TODO check prev actually could never be Ghost right?
                vertex.hit_ghost_times += ray.num_ghost_hits;
                vertex.ghost_ids = ray.ghost_shape_ids;
                vertex.ghost_hit_positions = ray.ghost_hit_positions;
                vertex.ghost_hit_geom_normals = ray.ghost_hit_geom_normals;
            }
            ray.num_ghost_hits = 0;
        }
        Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);

        std::optional<BSDFSampleRecord> bsdf_sample_ =
            sample_bsdf(mat,
                        wi,
                        vertex,
                        scene.texture_pool,
                        bsdf_rnd_param_uv,
                        bsdf_rnd_param_w,
                        mode);
        if (!bsdf_sample_) // if no valid sample; would rather not push in the vertex???
            break;
        const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
        Vector3 dirBSDF = bsdf_sample.dir_out;
        // Update ray differentials & eta_scale
        if (bsdf_sample.eta == 0)
        {
            ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
        }
        else
        {
            ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
            eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
        }

        Spectrum f = eval(mat, wi, dirBSDF, vertex, scene.texture_pool); // this includes the cosTheta
        if(f.x<0||f.y<0||f.z<0) break;//don't about the hell of this numerical texture error
        Real pdfBsdfW = pdf_sample_bsdf(mat, wi, dirBSDF, vertex, scene.texture_pool);
        if (pdfBsdfW <= 0)
            break;
        {
            vertex.type = VertexType::Surface;
            vertex.roughness = bsdf_sample.roughness;
            vertex.transport.contrib = throughput;
            vertex.transport.pathLength = prev.transport.pathLength + 1;
            // convert pdfW to pdfA
            Real G1 = fabs(dot(wi, vertex.geometric_normal)) / length_squared(prev.position - vertex.position);
            vertex.transport.pdfFwd = pdf * G1;
            Real pdfBsdfW_rev = pdf_sample_bsdf(mat, dirBSDF, wi, vertex, scene.texture_pool);
            Real G2 = fabs(dot(wi, prev.geometric_normal)) / length_squared(prev.position - vertex.position);
            prev.transport.pdfRev = pdfBsdfW_rev * G2;
            path.push_back(vertex);
        }

        // update information for loop
        pdf = pdfBsdfW; // record to store in next
        throughput *= f / pdfBsdfW;
        nextDir = dirBSDF;
    }
}

/**
 * @brief initialize the light subpaths used in bidirectional method: bdpt or vcm
 * @param scene
 * @param k at most k vertices
 * ///TODO stratify the starting point on the light source!!!
 */

#pragma once

#include "scene.h"
#include "pcg.h"
#include "parallel.h"
#include "progress_reporter.h"
#include "transform.h"

//////////////////////////////////////////////////////////////////////
// a clever trick pulled from mitsuba
template <typename Type>
class ScopedAssignment
{
public:
    // ScopedAssignment Public Methods
    ScopedAssignment(Type *target = nullptr, Type value = Type())
        : target(target)
    {
        if (target)
        {
            backup = *target;
            *target = value;
        }
    }
    ~ScopedAssignment()
    {
        if (target)
            *target = backup;
    }
    ScopedAssignment(const ScopedAssignment &) = delete;
    ScopedAssignment &operator=(const ScopedAssignment &) = delete;
    ScopedAssignment &operator=(ScopedAssignment &&other)
    {
        if (target)
            *target = backup;
        target = other.target;
        backup = other.backup;
        other.target = nullptr;
        return *this;
    }

private:
    Type *target, backup;
};
//////////////////////////////////////////////////////////////////////
// given the sampled light point, output Le * G similar to directIllumination function
Spectrum directIlluminationGivenLight(const PathVertex &vertex, 
                                      PathVertex &sampledVLight, //to be modified
                                      const Scene &scene, 
                                      int time, 
                                      bool testGhost = false, 
                                      bool testVisibility = false)
{
    const Light &light = scene.lights[sampledVLight.shape_id]; // the light_id is saved inside
    PointAndNormal point_on_light{sampledVLight.position, sampledVLight.geometric_normal};
    int shadow_hit_ghost_times = 0;
    Real G = 0;
    Vector3 dir2Light;
    if (!is_envmap(light))
    {
        dir2Light = normalize(sampledVLight.position - vertex.position);
        Ray shadow_ray{vertex.position, dir2Light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(sampledVLight.position, vertex.position),
                       time};
        if (!testGhost)
        { // occlusion test
            if (testVisibility)
            {
                if (!occluded(scene, shadow_ray))
                {
                    G = fabs(dot(-dir2Light, point_on_light.normal)) /
                        distance_squared(point_on_light.position, vertex.position);
                    G *= fabs(dot(dir2Light, vertex.geometric_normal));
                }
            }
            else
            {
                G = fabs(dot(-dir2Light, point_on_light.normal)) /
                    distance_squared(point_on_light.position, vertex.position);
                G *= fabs(dot(dir2Light, vertex.geometric_normal));
            }
        }
        else
        { // use intersection test
            std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, RayDifferential{});
            if (!vertex_) // not occluded
            {
                G = fabs(dot(-dir2Light, point_on_light.normal)) /
                    distance_squared(point_on_light.position, vertex.position);
                G *= fabs(dot(dir2Light, vertex.geometric_normal));

                shadow_hit_ghost_times = shadow_ray.num_ghost_hits;
                sampledVLight.ghost_ids = shadow_ray.ghost_shape_ids;
                sampledVLight.ghost_hit_positions = shadow_ray.ghost_hit_positions;
                sampledVLight.ghost_hit_geom_normals = shadow_ray.ghost_hit_geom_normals;
                shadow_ray.num_ghost_hits = 0;
            }
        }
    }
    else
    {
        dir2Light = -point_on_light.normal;
        Ray shadow_ray{vertex.position, dir2Light,
                       get_shadow_epsilon(scene),
                       infinity<Real>() /* envmaps are infinitely far away */, time};
        if (!testGhost)
        {
            if (testVisibility)
            {
                if (!occluded(scene, shadow_ray))
                {
                    // We integrate envmaps using the solid angle measure, so the geometry term is 1.
                    G = 1 * (vertex.type == VertexType::Ghost) ? fabs(dot(dir2Light, vertex.geometric_normal)) : max(dot(dir2Light, vertex.geometric_normal), Real(0));
                    // TODO for envmap
                }
            }
            else
            {
                G = 1 * (vertex.type == VertexType::Ghost) ? fabs(dot(dir2Light, vertex.geometric_normal)) : max(dot(dir2Light, vertex.geometric_normal), Real(0));
            }
        }
        else
        {
            std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, RayDifferential{});
            if (!vertex_)
            {
                G = 1 * (vertex.type == VertexType::Ghost) ? fabs(dot(dir2Light, vertex.geometric_normal)) : max(dot(dir2Light, vertex.geometric_normal), Real(0));
                shadow_hit_ghost_times = shadow_ray.num_ghost_hits;
                sampledVLight.ghost_ids = shadow_ray.ghost_shape_ids;
                sampledVLight.ghost_hit_positions = shadow_ray.ghost_hit_positions;
                sampledVLight.ghost_hit_geom_normals = shadow_ray.ghost_hit_geom_normals;
                shadow_ray.num_ghost_hits = 0;
            }
        }
    }
    if (G > 0 && sampledVLight.transport.pdfOriginA > 0)
    {
        Spectrum L = emission(light, -dir2Light, Real(0), point_on_light, scene);
        Spectrum C1 = G * L / sampledVLight.transport.pdfOriginA;
        return C1;
    }
    return make_zero_spectrum();
}

Spectrum directIllumination(const PathVertex &vertex,
                            Vector3 &dir_light,
                            PathVertex &sampledVLight, // to retrieve the information from the sampled vertex
                            pcg32_state &rng,
                            const Scene &scene,
                            int time, // TODO check all the time
                            bool testGhost = false)
{
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light =
        sample_point_on_light(light, vertex.position, light_uv, shape_w, scene);

    Vector3 dir2Light;
    Real G = 0;
    int shadow_hit_ghost_times = 0;
    if (is_envmap(light))
    {
        dir2Light = -point_on_light.normal;
        Ray shadow_ray{vertex.position, dir2Light,
                       get_shadow_epsilon(scene),
                       infinity<Real>() /* envmaps are infinitely far away */, time};
        if (!testGhost)
        {
            if (!occluded(scene, shadow_ray))
            {
                // We integrate envmaps using the solid angle measure, so the geometry term is 1.
                G = 1 * (vertex.type == VertexType::Ghost) ? fabs(dot(dir2Light, vertex.geometric_normal)) : max(dot(dir2Light, vertex.geometric_normal), Real(0));
            }
        }
        else
        {
            std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, RayDifferential{});
            if (!vertex_)
            {
                G = 1 * (vertex.type == VertexType::Ghost) ? fabs(dot(dir2Light, vertex.geometric_normal)) : max(dot(dir2Light, vertex.geometric_normal), Real(0));
                shadow_hit_ghost_times = shadow_ray.num_ghost_hits;
                sampledVLight.ghost_ids = shadow_ray.ghost_shape_ids;
                sampledVLight.ghost_hit_positions = shadow_ray.ghost_hit_positions;
                sampledVLight.ghost_hit_geom_normals = shadow_ray.ghost_hit_geom_normals;
                shadow_ray.num_ghost_hits = 0;
            }
        }
    
    }else{
        dir2Light = normalize(point_on_light.position - vertex.position); // + 1e-1 * dir2Light
        Ray shadow_ray{vertex.position, dir2Light,
                       get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) *
                           distance(point_on_light.position, vertex.position),
                       time};
        if (!testGhost)
        { // occlusion tests
            if (!occluded(scene, shadow_ray))
            {
                G = fabs(dot(-dir2Light, point_on_light.normal)) /
                    distance_squared(point_on_light.position, vertex.position);
                G *= fabs(dot(dir2Light, vertex.geometric_normal));
            }
        }
        else
        { // use intersection test
            std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, RayDifferential{});
            
            if (!vertex_) // not occluded
            {
                G = fabs(dot(-dir2Light, point_on_light.normal)) /
                    distance_squared(point_on_light.position, vertex.position);
                G *= fabs(dot(dir2Light, vertex.geometric_normal));

                shadow_hit_ghost_times = shadow_ray.num_ghost_hits;
                sampledVLight.ghost_ids = shadow_ray.ghost_shape_ids;
                sampledVLight.ghost_hit_positions = shadow_ray.ghost_hit_positions;
                sampledVLight.ghost_hit_geom_normals = shadow_ray.ghost_hit_geom_normals;
                shadow_ray.num_ghost_hits = 0;
            }
        }
    }
    
    Real pdfLight = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, point_on_light, vertex.position, scene);

    sampledVLight.type = VertexType::Light;
    sampledVLight.time = time;
    sampledVLight.hit_ghost_times = shadow_hit_ghost_times;
    sampledVLight.position = point_on_light.position;
    sampledVLight.geometric_normal = point_on_light.normal;
    sampledVLight.shape_id = light_id; // FIXME
    dir_light = dir2Light;

    if (G > 0 && pdfLight > 0)
    {
        Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
        Spectrum C1 = G * L / pdfLight; // TODO attention doesn't include f
        sampledVLight.transport.contrib = L / pdfLight;
        sampledVLight.transport.pathLength = 1; // TODO record how many vertices are there in the path
        sampledVLight.transport.pdfOriginA = pdfLight;
        return C1;
    }
    return make_zero_spectrum(); // forgot this, debug for long!!!
}

Real pdf_bsdf_at_vertex(const PathVertex &vertex,
                        const PathVertex &next,
                        const PathVertex &prev,
                        const TransportDirection &mode,
                        const Scene &scene)
{
    const Material &mat = scene.materials[vertex.material_id];
    const Vector3 &wi = normalize(prev.position - vertex.position);
    const Vector3 &wo = normalize(next.position - vertex.position);
    Real pdfBsdfW = pdf_sample_bsdf(mat, wi, wo, vertex, scene.texture_pool);
    return pdfBsdfW;
}

Real pdfDirLightA(const PathVertex &vLight, const PathVertex &vRef, const Light &light, const Scene &scene)
{
    Vector3 dirFromLight = vRef.position - vLight.position;
    Real dist = length(dirFromLight);
    dirFromLight /= dist;
    Real pdfLightDirW = pdf_dir_from_light(
        light, PointAndNormal{vLight.position, vLight.geometric_normal},
        dirFromLight, scene);
    Real G = fabs(dot(dirFromLight, vRef.geometric_normal)) / dist / dist;
    return pdfLightDirW * G;
}

inline Real MIS(Real pdf)
{
    return pdf;
}
