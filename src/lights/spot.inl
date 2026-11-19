#include "utils.h"

Real light_power_op::operator()(const Spot &light) const
{
    return luminance(light.I) * 2.0 * M_PI * (1.0 - 0.5 * (light.cosFalloffStart + light.cosTotalWidth));
}

PointAndNormal sample_point_on_light_op::operator()(const Spot &light) const
{

    return PointAndNormal{light.pLight, light.normal};
}

Real pdf_point_on_light_op::operator()(const Spot &light) const
{
    return 1.0;
}

Spectrum emission_op::operator()(const Spot &light) const
{
    Real cosTheta = dot(view_dir, light.normal);
    if(cosTheta < light.cosTotalWidth) return make_zero_spectrum();
    if(cosTheta >= light.cosFalloffStart) return light.I;
    //compute the falloff inside spotlight cone
    Real delta = (cosTheta - light.cosTotalWidth) / (light.cosFalloffStart - light.cosTotalWidth);
    Real falloff = (delta * delta) * (delta * delta);

    return light.I * falloff;
}

//only needed in light tracing and bdpt
Vector3 sample_dir_from_light_op::operator()(const Spot &light) const
{
    return {0.,0.,0.}; //TODO; just placeholder
}

Real pdf_dir_from_light_op::operator()(const Spot &light) const{
    return 1.0;//TODO; just placeholder
}

void init_sampling_dist_op::operator()(Spot &light) const
{
}
