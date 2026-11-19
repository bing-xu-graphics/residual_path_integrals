#pragma once

#include "lajolla.h"
#include "frame.h"
#include "ray.h"
#include "spectrum.h"
#include "vector.h"

#include <optional>

struct Scene;

///////////////////////embree behavior modifier//////////////

#include <embree3/rtcore.h>

struct Payload{
    unsigned int num_ghost_hits = 0; //seems for occlusion ,can directly use this
    std::vector<Real> dist_t;//record all the distance t for all the ghost intersection to illiminate those farther than the cloestHit 
    std::vector<int> ghost_shape_ids;
    std::vector<Vector3> ghost_hit_positions;
    std::vector<Vector3> ghost_hit_geom_normals;
};

//self defined ray payload for the intersectioncontex
struct IntersectContext{
    RTCIntersectContext rtc_context;
    //some other payloads below:
    void* payload;  
    IntersectContext(){
        payload = new Payload();
    }
     ~IntersectContext() {
        delete static_cast<Payload *>(payload);
    }
};

inline void InitIntersectionContext(struct IntersectContext* context)
{
  rtcInitIntersectContext(&context->rtc_context);
}

////////////////////////////////////////////////////

//information for light transport
enum struct VertexType { Camera, Light, Surface, Medium, Ghost, DUMMY };

struct TransportInfo{
    
    Spectrum contrib = make_zero_spectrum(); //the accumulated throughput/contribution L * (brdf * G)from 0 to n-1 * We 
    int pathLength; //number of segments between source and current vertex
    /* pdfFwd is for the pdfW of how this current vertex is sampled from the prev. (in the direction of random walk)
       pdfRev is for the pdfW of how this current vertex is sampled from the next. (reverse random walk direction)*/
    Real pdfFwd = 1.0, pdfRev = 1.0; //accumulated pdf from 0 to n-1 vertices
    ///////////// for incremental rendering 
    Real pdfOriginA = 0, pdfDir = 0; //for the inital sampled vertex in the middle

    //@copied from vcm
    Real dVCM = 0; // MIS quantity used for vertex connection and merging
    Real dVC = 0; // MIS quantity used for vertex connection

    //constructor
    TransportInfo() = default;
    TransportInfo(const Spectrum &contrib, int pathLength, Real pdfFwd, Real pdfRev, Real pdfOriginA, Real pdfDir)
        : contrib(contrib), pathLength(pathLength), pdfFwd(pdfFwd), pdfRev(pdfRev), pdfOriginA(pdfOriginA), pdfDir(pdfDir)
    {
    }
    //copy constructor
    TransportInfo(const TransportInfo &t)
        : contrib(t.contrib), pathLength(t.pathLength), pdfFwd(t.pdfFwd), 
        pdfRev(t.pdfRev), pdfOriginA(t.pdfOriginA), pdfDir(t.pdfDir), dVCM(t.dVCM), dVC(t.dVC)
    {
    }
    
};

/// An "PathVertex" represents a vertex of a light path.
/// We store the information we need for computing any sort of path contribution & sampling density.
struct PathVertex {
    //all local information
    Vector3 position;
    Vector3 geometric_normal; // always face at the same direction at shading_frame.n
    Frame shading_frame;
    Vector2 st; // A 2D parametrization of the surface. Irrelavant to UV mapping.
                // for triangle this is the barycentric coordinates, which we use
                // for interpolating the uv map.
    Vector2 uv; // The actual UV we use for texture fetching.
    // For texture filtering, stores approximatedly min(abs(du/dx), abs(dv/dx), abs(du/dy), abs(dv/dy))
    Real uv_screen_size;
    Real mean_curvature; // For ray differential propagation.
    Real ray_radius; // For ray differential propagation.
    int shape_id = -1;
    int primitive_id = -1; // For triangle meshes. This indicates which triangle it hits.
    int material_id = -1; //can use this to retrieve the BSDF

    // If the path vertex is inside a medium, these two IDs
    // are the same.
    int interior_medium_id = -1;
    int exterior_medium_id = -1;

    /// for some light transport method like BDPT
    struct TransportInfo transport;
    //information related to incremental rendering
    VertexType type = VertexType::DUMMY;//light or sensor or surface or medium interaction
    int dynamic_id = -1;
    std::vector<int> ghost_ids; //record all the shape_id of the ghost objects by the segment ending here; (along the random walk direction)
    std::vector<Vector3> ghost_hit_positions;
    std::vector<Vector3> ghost_hit_geom_normals;
    int hit_dynamic_times = 0;
    int hit_ghost_times = 0;
    int time = 0;//match the ray from which the intersection is generated
    Real roughness = 0; //for the roughness of the material sampled at this vertex

    //constructor
    PathVertex() = default;
    PathVertex(const Vector3 &position,
               const Vector3 &geometric_normal,
               const Frame &shading_frame,
               const Vector2 &st,
               const Vector2 &uv,
               Real uv_screen_size,
               Real mean_curvature,
               Real ray_radius,
               int shape_id,
               int primitive_id,
               int material_id,
               int interior_medium_id,
               int exterior_medium_id,
               TransportInfo transport,
               VertexType type,
               int dynamic_id,
               std::vector<int> ghost_ids,
               std::vector<Vector3> ghost_hit_positions,
               std::vector<Vector3> ghost_hit_geom_normals,
               int hit_dynamic_times,
               int hit_ghost_times,
               int time,
               Real roughness)
        : position(position),
          geometric_normal(geometric_normal),
          shading_frame(shading_frame),
          st(st),
          uv(uv),
          uv_screen_size(uv_screen_size),
          mean_curvature(mean_curvature),
          ray_radius(ray_radius),
          shape_id(shape_id),
          primitive_id(primitive_id),
          material_id(material_id),
          interior_medium_id(interior_medium_id),
          exterior_medium_id(exterior_medium_id),
          transport(transport),
          type(type),
          dynamic_id(dynamic_id),
          ghost_hit_positions(ghost_hit_positions),
          ghost_hit_geom_normals(ghost_hit_geom_normals),
          hit_dynamic_times(hit_dynamic_times),
          hit_ghost_times(hit_ghost_times),
          time(time),
          roughness(roughness)
    {
    }
    //copy constructor
    PathVertex(const PathVertex &v)
        : position(v.position),
          geometric_normal(v.geometric_normal),
          shading_frame(v.shading_frame),
          st(v.st),
          uv(v.uv),
          uv_screen_size(v.uv_screen_size),
          mean_curvature(v.mean_curvature),
          ray_radius(v.ray_radius),
          shape_id(v.shape_id),
          primitive_id(v.primitive_id),
          material_id(v.material_id),
          interior_medium_id(v.interior_medium_id),
          exterior_medium_id(v.exterior_medium_id),
          transport(v.transport),
          type(v.type),
          dynamic_id(v.dynamic_id),
          ghost_ids(v.ghost_ids),
          ghost_hit_positions(v.ghost_hit_positions),
          ghost_hit_geom_normals(v.ghost_hit_geom_normals),
          hit_dynamic_times(v.hit_dynamic_times),
          hit_ghost_times(v.hit_ghost_times),
          time(v.time),
          roughness(v.roughness)
    {
    }
};

std::optional<PathVertex> intersect(const Scene &scene,
                                    Ray &ray,
                                    const RayDifferential &ray_diff = RayDifferential{});

/// Test is a ray segment intersect with anything in a scene.
bool occluded(const Scene &scene, Ray &ray);

/// Computes the emission at a path vertex v, with the viewing direction
/// pointing outwards of the intersection.
Spectrum emission(const PathVertex &v,
                  const Vector3 &view_dir,
                  const Scene &scene);

std::optional<PathVertex> intersect_dynamic(const Scene &scene,
                                    Ray &ray,
                                    const RayDifferential &ray_diff = RayDifferential{});

bool occluded_dynamic(const Scene &scene, Ray &ray);

bool intersect_ghost(const Scene &scene, Ray &ray);
