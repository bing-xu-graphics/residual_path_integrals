#include "utils.h"
Real light_power_op::operator()(const DiffuseAreaLight &light) const
{
    return luminance(light.intensity) * surface_area(scene.shapes[light.shape_id]) * c_PI;
}

PointAndNormal sample_point_on_light_op::operator()(const DiffuseAreaLight &light) const
{
    const Shape &shape = scene.shapes[light.shape_id];
    return sample_point_on_shape(shape, ref_point, rnd_param_uv, rnd_param_w);
}

Vector3 sample_dir_from_light_op::operator()(const DiffuseAreaLight &light) const
{

    Real pdf;
    Vector3 dirLocal = sample_cos_hemisphere(rnd_param_uv, pdf);//in local coordinates
    dirLocal.z = max(dirLocal.z, (Real)EPS_COSINE);
    Frame frame(point_on_light.normal);
    return to_world(frame, dirLocal);
}

Real pdf_dir_from_light_op::operator()(const DiffuseAreaLight &light) const{
    Frame frame(point_on_light.normal);
    Vector3 dirLocal = to_local(frame, light_dir);
    return c_INVPI * max(dirLocal.z, (Real)EPS_COSINE);
}

Real pdf_point_on_light_op::operator()(const DiffuseAreaLight &light) const
{
    return pdf_point_on_shape(
        scene.shapes[light.shape_id], point_on_light, ref_point);
}

Spectrum emission_op::operator()(const DiffuseAreaLight &light) const
{
    if (dot(point_on_light.normal, view_dir) <= 0)
    {
        return make_zero_spectrum();
    }
    return light.intensity;
}

void init_sampling_dist_op::operator()(DiffuseAreaLight &light) const
{
}
