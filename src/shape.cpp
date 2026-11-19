#include "shape.h"
#include "intersection.h"
#include "point_and_normal.h"
#include "ray.h"
#include <embree3/rtcore.h>
#include "utils.h"

struct register_embree_op
{
    uint32_t operator()(const Sphere &sphere) const;
    uint32_t operator()(const TriangleMesh &mesh) const;

    const RTCDevice &device;
    const RTCScene &scene;
};

struct register_embree_byID_op
{
    uint32_t operator()(const Sphere &sphere) const;
    uint32_t operator()(const TriangleMesh &mesh) const;

    const RTCDevice &device;
    const RTCScene &scene;
    const uint32_t &geomID;
};

struct sample_point_on_shape_op
{
    PointAndNormal operator()(const Sphere &sphere) const;
    PointAndNormal operator()(const TriangleMesh &mesh) const;

    const Vector3 &ref_point;
    const Vector2 &uv; // for selecting a point on a 2D surface
    const Real &w;     // for selecting triangles
    int time = 1;
};

struct surface_area_op
{
    Real operator()(const Sphere &sphere) const;
    Real operator()(const TriangleMesh &mesh) const;
};

struct pdf_point_on_shape_op
{
    Real operator()(const Sphere &sphere) const;
    Real operator()(const TriangleMesh &mesh) const;

    const PointAndNormal &point_on_shape;
    const Vector3 &ref_point;
    int time = 1;
};

struct init_sampling_dist_op
{
    void operator()(Sphere &sphere) const;
    void operator()(TriangleMesh &mesh) const;
};

struct compute_shading_info_op
{
    ShadingInfo operator()(const Sphere &sphere) const;
    ShadingInfo operator()(const TriangleMesh &mesh) const;

    const PathVertex &vertex;
};

#include "shapes/sphere.inl"
#include "shapes/triangle_mesh.inl"

uint32_t register_embree(const Shape &shape, const RTCDevice &device, const RTCScene &scene)
{
    return std::visit(register_embree_op{device, scene}, shape);
}

PointAndNormal sample_point_on_shape(const Shape &shape,
                                     const Vector3 &ref_point,
                                     const Vector2 &uv,
                                     Real w)
{
    return std::visit(sample_point_on_shape_op{ref_point, uv, w}, shape);
}

Real pdf_point_on_shape(const Shape &shape,
                        const PointAndNormal &point_on_shape,
                        const Vector3 &ref_point)
{
    return std::visit(pdf_point_on_shape_op{point_on_shape, ref_point}, shape);
}

Real surface_area(const Shape &shape)
{
    return std::visit(surface_area_op{}, shape);
}

void init_sampling_dist(Shape &shape)
{
    return std::visit(init_sampling_dist_op{}, shape);
}

ShadingInfo compute_shading_info(const Shape &shape, const PathVertex &vertex)
{
    return std::visit(compute_shading_info_op{vertex}, shape);
}

///////////for incremental rendering

Vector3  sample_dir_from_shape(const Shape &shape,
                              const PathVertex &vertex,
                              const Vector2 &uv,
                              Real &pdf,
                              int time)
{
    Vector3 dirLocal = sample_cos_hemisphere(uv, pdf); // in local coordinates
    dirLocal.z = max(dirLocal.z, (Real)EPS_COSINE);
    Frame frame = vertex.shading_frame; // FIXME what to do?
    Vector3 dirWorld = to_world(frame, dirLocal);
    return dirWorld;
}

Real pdf_dir_from_shape(const Shape &shape,
                        const PathVertex &vertex,
                        const Vector3 &sampled_dir,
                        int time)
{
    Vector3 dirLocal = to_local(vertex.shading_frame, sampled_dir);
    return c_INVPI * max(dirLocal.z, (Real)EPS_COSINE);
}

///////////some counterpart with bounded timestamp t
PointAndNormal sample_point_on_shape(const Shape &shape,
                                     const Vector3 &ref_point,
                                     const Vector2 &uv,
                                     Real w, int timestamp)
{
    return std::visit(sample_point_on_shape_op{ref_point, uv, w}, shape);
}

Real pdf_point_on_shape(const Shape &shape,
                        const PointAndNormal &point_on_shape,
                        const Vector3 &ref_point, int timestamp)
{
    return std::visit(pdf_point_on_shape_op{point_on_shape, ref_point, timestamp}, shape);
}

#include "dynamic.h"
#include "shapes/dynamic.inl"

uint32_t register_embree_dynamic(const Shape &shape,
                                 const RTCDevice &device,
                                 const RTCScene &scene,
                                 const int time_steps,
                                 Dynamic &dynamic)
{
    return std::visit(register_embree_dynamic_op{device, scene, time_steps, dynamic}, shape);
}

#include "shapes/ghost.inl"

uint32_t register_embree_ghost_dynamic(const Shape &shape,
                                       const RTCDevice &device,
                                       const RTCScene &scene,
                                       const RTCScene &ghost_scene,
                                       const int time_steps,
                                       Dynamic &dynamic)
{
    uint32_t ghost_geomID =  std::visit(register_embree_ghost_dynamic_op{device, scene, ghost_scene, time_steps, dynamic}, shape);
    //prepare the ghost id
    return ghost_geomID;
}

///////////// sampling helpers/////

PointAndNormal sample_point_on_dynamic(const Dynamic &dynamic,
                                       const Vector3 &ref_point,
                                       const Vector2 &uv,
                                       Real w,
                                       int time)
{
    if (time == NEW_TIME)
    {
        const Shape &shape = dynamic.base_shape;
        return std::visit(sample_point_on_shape_op{ref_point, uv, w}, shape);
    }
    else if (time == OLD_TIME)
    {
        const Shape &shape = dynamic.old_shape;
        return std::visit(sample_point_on_shape_op{ref_point, uv, w}, shape);
    }
    else
    {
        printf("[ERROR] unsupported time!!!!\n");
        assert(false);
        return PointAndNormal{};
    }
}

Real pdf_point_on_dynamic(const Dynamic &dynamic,
                          const PointAndNormal &point_on_shape,
                          const Vector3 &ref_point,
                          int time)
{
    if (time == NEW_TIME)
    {
        const Shape &shape = dynamic.base_shape;
        return std::visit(pdf_point_on_shape_op{point_on_shape, ref_point}, shape);
    }
    else if (time == OLD_TIME)
    {
        const Shape &shape = dynamic.old_shape;
        return std::visit(pdf_point_on_shape_op{point_on_shape, ref_point}, shape);
    }
    else
    {
        printf("[ERROR] unsupported time!!!!\n");
        assert(false);
        return 0.;
    }
}

PointAndNormal sample_point_on_ghost(const Dynamic &dynamic,
                                     const Vector3 &ref_point,
                                     const Vector2 &uv,
                                     Real w,
                                     int time)
{
    if (time == NEW_TIME)
    {
        const Shape &shape = dynamic.old_shape;
        return std::visit(sample_point_on_shape_op{ref_point, uv, w}, shape);
    }
    else if (time == OLD_TIME)
    {
        const Shape &shape = dynamic.base_shape;
        return std::visit(sample_point_on_shape_op{ref_point, uv, w}, shape);
    }
    else
    {
        printf("[ERROR] sample_point_on_ghost() unsupported time!!!!\n");
        assert(false);
        return PointAndNormal{};
    }
}

Real pdf_point_on_ghost(const Dynamic &dynamic,
                        const PointAndNormal &point_on_shape,
                        const Vector3 &ref_point,
                        int time)
{
    if (time == NEW_TIME)
    {
        const Shape &shape = dynamic.old_shape;
        return std::visit(pdf_point_on_shape_op{point_on_shape, ref_point}, shape);
    }
    else if (time == OLD_TIME)
    {
        const Shape &shape = dynamic.base_shape;
        return std::visit(pdf_point_on_shape_op{point_on_shape, ref_point}, shape);
    }
    else
    {
        printf("[ERROR] pdf_point_on_ghost() unsupported time!!!!\n");
        assert(false);
        return 0.;
    }
}
