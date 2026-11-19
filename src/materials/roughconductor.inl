#include "../microfacet.h"

Spectrum eval_op::operator()(const RoughConductor &bsdf) const {

    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // The half-vector is a crucial component of the microfacet models.
    // Since microfacet assumes that the surface is made of many small mirrors/glasses,
    // The "average" between input and output direction determines the orientation
    // of the mirror our ray hits (since applying reflection of dir_in over half_vector
    // gives us dir_out). Microfacet models build all sorts of quantities based on the
    // half vector. It's also called the "micro normal".
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return make_zero_spectrum();
    }

    Spectrum Ks = eval(
        bsdf.specular_reflectance, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x = eval(
        bsdf.m_alphaU, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_y = eval(
        bsdf.m_alphaV, vertex.uv, vertex.uv_screen_size, texture_pool);
        
    // Clamp roughness to avoid numerical issues.
    alpha_x = std::clamp(alpha_x, Real(0.0001), Real(1));
    alpha_y = std::clamp(alpha_y, Real(0.0001), Real(1));

    // For D, they use an anisotropic Trowbridge-Reitz distribution (aka GGX)
    Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);
    if(D == 0) return make_zero_spectrum();
    /* fresnel factor*/
    Real h_dot_in = dot(half_vector, dir_in);
    Spectrum F = fresnel_conductor(h_dot_in, bsdf.eta, bsdf.m_k) * Ks;

    // For G, they use the Smith masking term corresponds to GTR2
    Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
    Real G_out = smith_masking_gtr2(to_local(frame, dir_out), alpha_x, alpha_y);
    Real G = G_in * G_out;

    return F * D * G / (4 * dot(dir_in, vertex.shading_frame.n));
    
}

Real pdf_sample_bsdf_op::operator()(const RoughConductor &bsdf) const {
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }

    Spectrum S = eval(
        bsdf.specular_reflectance, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x = eval(
        bsdf.m_alphaU, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_y = eval(
        bsdf.m_alphaV, vertex.uv, vertex.uv_screen_size, texture_pool);
        
    // Clamp roughness to avoid numerical issues.
    alpha_x = std::clamp(alpha_x, Real(0.0001), Real(1));
    alpha_y = std::clamp(alpha_y, Real(0.0001), Real(1));
    
    Real D = GTR2(to_local(frame, half_vector), alpha_x, alpha_y);

    if (luminance(S) <= 0) {
        return 0;
    }

    Real G_in = smith_masking_gtr2(to_local(frame, dir_in), alpha_x, alpha_y);
    return D * G_in / (4 * dot(dir_in, vertex.shading_frame.n));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const RoughConductor &bsdf) const {
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // We use the reflectance to choose between sampling the dielectric or diffuse layer.
    Spectrum Ks = eval(
        bsdf.specular_reflectance, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x = eval(
        bsdf.m_alphaU, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_y = eval(
        bsdf.m_alphaV, vertex.uv, vertex.uv_screen_size, texture_pool);
        
    // Clamp roughness to avoid numerical issues.
    alpha_x = std::clamp(alpha_x, Real(0.0001), Real(1));
    alpha_y = std::clamp(alpha_y, Real(0.0001), Real(1));

    Vector3 local_dir_in = to_local(vertex.shading_frame, dir_in);
 
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

// Transform the micro normal to world space
    Vector3 half_vector = to_world(vertex.shading_frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    Real roughness = 0.5f * alpha_x + alpha_y;
    return BSDFSampleRecord{
        reflected, Real(0) /* eta */, roughness /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const RoughConductor &bsdf) const {
    return bsdf.specular_reflectance;
}
