
#pragma once

#include "scene.h"
#include "pcg.h"
#include "parallel.h"
#include "progress_reporter.h"
#include "transform.h"
#include "utils.h"

/**Specifies the roughness threshold for classifying materials as 'diffuse', in contrast to 'specular',
 *			for the purposes of constructing paths pairs .This value should usually be somewhere between 0.0005 and 0.01.
 *			If the result image has noise similar to standard path tracing, increasing or decreasing this value may sometimes help. This implementation assumes that this value is small.\default{\code{0.001}}
 */
#define SM_THRESHOLD 0.001

#define CONNECTIBLE_DISTANCE_THRESHOLD 0.01

/**
 * @brief for mis the fd calculated from new2old and old2new
 *
 * @param timestamp
 * @return Real
 */
Real compute_mis_weight_frame(int timestamp, bool failed_sm, Real basePdf, Real shiftedPdf)
{   
    Real mis_weight = safe_divide(basePdf,basePdf + shiftedPdf) ;
    return timestamp == NEW_TIME ? (failed_sm ? 1.0 : mis_weight) : (failed_sm ? -1.0 : -mis_weight);
}

Real compute_mis_weight_frame(int timestamp, bool failed_sm)
{   
    return timestamp == NEW_TIME ? (failed_sm ? 1.0 : 0.5) : (failed_sm ? -1.0 : -0.5);
}

enum struct SHIFT_MAP_STATUS
{
    DIVERGED,
    CONNECTED, ///< Connected, allows using BSDF values from the base path.
    RECENTLY_CONNECTED, ///< Connected, but different incoming direction so needs a BSDF evaluation.
    FAILED, // if failed, the weight for forward and reverse shift mapping (for two frames) is 0.5; otherwise it's 1.0 which is exactly finite difference for two random paths
    NEVER_DIVERGED
};

inline Vector3 division_protect(Vector3 v)
{
    return Vector3{max(v.x, EPS_COSINE), max(v.y, EPS_COSINE), max(v.z, EPS_COSINE)};
}

bool connectVertices(const PathVertex &vertex1, PathVertex &vertex2, const Scene &scene, int time, bool fillHits = false)
{
    Vector3 dir12 = vertex2.position - vertex1.position;
    Real dist = length(dir12);
    assert(dist > 1e-3);
    dir12 /= dist;
    if (dot(dir12, vertex1.geometric_normal) <= 0)
        return false;
    if (dot(-dir12, vertex2.geometric_normal) <= 0)
        return false; // TODO change to eps_cos??
    // visibility test
    if (fillHits)
    {
        Ray shadow_ray{vertex1.position, dir12, get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) * dist, time};
        std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, RayDifferential{});
        if (vertex_)
            return false; // hit something
        vertex2.hit_ghost_times = vertex1.hit_ghost_times + shadow_ray.num_ghost_hits;
        // TODO what about the dynamic
        assert(vertex2.shape_id >= 0);
        vertex2.hit_dynamic_times = vertex1.hit_dynamic_times + (is_dynamic(scene.shapes[vertex2.shape_id]));
    }
    else
    {
        Ray shadow_ray{vertex1.position, dir12, get_shadow_epsilon(scene),
                       (1 - get_shadow_epsilon(scene)) * dist, time};
        if (occluded(scene, shadow_ray))
            return false; // not connectible
    }

    return true;
}

bool connectVertices(const PathVertex &vertex1, const PathVertex &vertex2, const Scene &scene, int time)
{
    Vector3 dir12 = vertex2.position - vertex1.position;
    Real dist = length(dir12);
    assert(dist > 1e-3);
    if (dist < 1e-3)
        return false;
    dir12 /= dist;
    if (dot(dir12, vertex1.geometric_normal) <= 0)
        return false;
    if (dot(-dir12, vertex2.geometric_normal) <= 0)
        return false; // TODO change to eps_cos??
                      // visibility test
    Ray shadow_ray{vertex1.position, dir12, get_shadow_epsilon(scene),
                   (1 - get_shadow_epsilon(scene)) * dist, time};
    if (occluded(scene, shadow_ray))
        return false; // not connectible
    return true;
}

/*
 fill in the informaiton of ghost intersections of the targetV, the ray from sourceV to targetV is already checked to be visible
 only need to test the intersection with ghost objects
*/
void updateGhostIsectInfo(PathVertex& targetV, const PathVertex& sourceV, const Scene& scene, int time)
{
    Vector3 dir12 = targetV.position - sourceV.position;
    Real dist = length(dir12);
    dir12 /= dist;
    Ray shadow_ray{sourceV.position, dir12, get_shadow_epsilon(scene),
                   (1 - get_shadow_epsilon(scene)) * dist, time};
    bool hitGhost = intersect_ghost(scene, shadow_ray);
    
    {
        targetV.hit_ghost_times = sourceV.hit_ghost_times + shadow_ray.num_ghost_hits; //this one doens't even matter any more
        // TODO what about the dynamic
        targetV.ghost_ids = shadow_ray.ghost_shape_ids;
        targetV.ghost_hit_positions = shadow_ray.ghost_hit_positions;
        targetV.ghost_hit_geom_normals = shadow_ray.ghost_hit_geom_normals;
    }
}
// if this segment is created not by tracing, but by reconnection, let's recompute the pdfFwd of current and pdfRev of prev
bool updatePdfsOfSegment(const PathVertex &prevprev, PathVertex &prev, 
                        PathVertex &current, const PathVertex &next, 
                        const Scene &scene, int time)
{
    //1)update pdfFwd of current
    Real pdfWCurrent = pdf_bsdf_at_vertex(prev, current, prevprev, TransportDirection::TO_VIEW, scene);
    Real pdfACurrent = pdf_convert_w_to_A(pdfWCurrent, prev, current);
    if(pdfACurrent<=0) return false;
    current.transport.pdfFwd = pdfACurrent;
    //2)update pdfRev of prev
    if(next.shape_id>=0){ //otherwise there's the last one of non-end vertex, just leave it to be 1.0
        Real pdfWPrev = pdf_bsdf_at_vertex(current, prev, next, TransportDirection::TO_LIGHT, scene);
        Real pdfAPrev = pdf_convert_w_to_A(pdfWPrev, current, prev);
        if(pdfAPrev<=0) return false;
        prev.transport.pdfRev = pdfAPrev;
    }
    
    return true;
}
 
//when the segment is the segment right next to the dynamic point. 
bool updatePdfsOfSegment(PathVertex &shifted_vDynamic, PathVertex &shifted_vDynamic_next, 
                        const PathVertex &nextnext,
                        const Scene &scene, int time)
{
    //1)update pdfFwd of shifted_vDynamic_next
    Real pdfWCurrent = shifted_vDynamic.transport.pdfDir;
    Real pdfACurrent = pdf_convert_w_to_A(pdfWCurrent, shifted_vDynamic, shifted_vDynamic_next);
    assert(pdfACurrent == pdfACurrent);
    if(pdfACurrent<=0) return false;
    shifted_vDynamic_next.transport.pdfFwd = pdfACurrent;
    //2)update pdfRev of shifted_vDynamic
    if(nextnext.shape_id>=0){ //otherwise there's the last one of non-end vertex, just leave it to be 1.0
        Real pdfWPrev = pdf_bsdf_at_vertex(shifted_vDynamic_next, shifted_vDynamic, nextnext, TransportDirection::TO_LIGHT, scene);
        Real pdfAPrev = pdf_convert_w_to_A(pdfWPrev, shifted_vDynamic_next, shifted_vDynamic);
        if(pdfAPrev<=0) return false;
        shifted_vDynamic.transport.pdfRev = pdfAPrev;
    }
    return true;
}

PathVertex get_vertex_at_time(const PathVertex &vertex, int timestamp, const Scene &scene)
{
    PathVertex shifted(vertex);
    ////////TODO change the position and normal of the shifted at the new timestamp, also the shading info???
    const Dynamic &dynamic = scene.dynamics[vertex.dynamic_id]; // TODO make sure that the dyamic_id is set in all the vertex
    const Matrix4x4 xform = timestamp == OLD_TIME ? dynamic.trans2old : dynamic.inv_trans2old;
    shifted.time = timestamp;
    shifted.position = xform_point(xform, vertex.position);
    shifted.geometric_normal = xform_normal(inverse(xform), vertex.geometric_normal);
    const Shape &shifted_shape = timestamp == OLD_TIME ? dynamic.old_shape : dynamic.base_shape;

    ShadingInfo shadingInfo = compute_shading_info(shifted_shape, shifted);
    shifted.shading_frame = shadingInfo.shading_frame;
    shifted.uv = shadingInfo.uv;
    return shifted;
}

PathVertex get_ghost_vertex_at_time(const PathVertex &vertex, int timestamp, const Scene &scene)
{
    PathVertex shifted(vertex);
    ////////TODO change the position and normal of the shifted at the new timestamp, also the shading info???
    const Dynamic &dynamic = scene.dynamics[vertex.dynamic_id]; // TODO make sure that the dyamic_id is set in all the vertex
    const Matrix4x4 xform = timestamp == NEW_TIME ? dynamic.trans2old : dynamic.inv_trans2old;

    shifted.time = timestamp;
    shifted.position = xform_point(xform, vertex.position);
    shifted.geometric_normal = xform_normal(inverse(xform), vertex.geometric_normal);
    const Shape &shifted_shape = timestamp == NEW_TIME ? dynamic.old_shape : dynamic.base_shape;

    ShadingInfo shadingInfo = compute_shading_info(shifted_shape, shifted);
    shifted.shading_frame = shadingInfo.shading_frame;
    shifted.uv = shadingInfo.uv;
    return shifted;
}

/// Classification of vertices into diffuse and glossy.

enum MATType {
    VERTEX_TYPE_GLOSSY,     ///< "Specular" vertex that requires the random seed replay or half-vector duplication shift.
    VERTEX_TYPE_DIFFUSE     ///< "Non-specular" vertex that is rough enough for the reconnection shift.
};

/// @modified from gradient-domain path tracing
/// Returns the vertex type (diffuse / glossy) of a vertex, for the purposes of determining
/// the shifting strategy. copied from gradient domain path tracing 
///
/// A bare classification by roughness alone is not good for multi-component BSDFs since they
/// may contain a diffuse component and a perfect specular component. If the base path
/// is currently working with a sample from a BSDF's smooth component, we don't want to care
/// about the specular component of the BSDF right now - we want to deal with the smooth component.
///
/// For this reason, we vary the classification a little bit based on the situation.
/// This is perfectly valid, and should be done.
MATType getVertexRoughnessType(const PathVertex &vertex, const Scene &scene) {
    // Return the lowest roughness value of the components of the vertex's BSDF.
    // If 'bsdfType' does not have a delta component, do not take perfect speculars (zero roughness) into account in this.
    Real lowest_roughness = std::numeric_limits<Real>::infinity();
    const Material &mat = scene.materials[vertex.material_id];
    bool found_smooth = false;
    bool found_dirac = false;
    Real component_roughness = vertex.roughness;
    if(component_roughness <= SM_THRESHOLD)
    {
        return VERTEX_TYPE_GLOSSY;
    }
    else
    {
        return VERTEX_TYPE_DIFFUSE;
    }

}

/**the reconnection creteria :test whether the shifted_prev can be reconnected to base_vertex
 * 1) shifted_prev and base_vertex should be at least 1e-3 distant apart.
 * 2)  shifted_prev and base_vertex should all be rough
 * 3) base_vertex must be static
*/
bool connectability_condition(const PathVertex& base_vertex, const PathVertex& base_prev, const PathVertex& shifted_prev, const Scene& scene){
    if (getVertexRoughnessType(base_vertex, scene) == MATType::VERTEX_TYPE_GLOSSY){
        return false;
    } 
    if (getVertexRoughnessType(shifted_prev, scene) == MATType::VERTEX_TYPE_GLOSSY){
        return false;
    } 
    if (getVertexRoughnessType(base_prev, scene) == MATType::VERTEX_TYPE_GLOSSY){
        return false;
    } 
    
    if(distance(base_prev.position, base_vertex.position) < CONNECTIBLE_DISTANCE_THRESHOLD){
        return false;
    }   
    if(distance(shifted_prev.position, base_vertex.position) < CONNECTIBLE_DISTANCE_THRESHOLD){
        return false;
    } 

    return true;
}
