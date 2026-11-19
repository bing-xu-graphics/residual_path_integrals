#pragma once

#include "vector.h"
#include "frame.h"
/// Sample direction in the upper hemisphere with cosine-proportional pdf
/** The returned PDF is with respect to solid angle measure */
#define EPS_COSINE (Real)1e-10f
#define EPS_RAY    (Real)1e-3f

////////////////////
// for incremental rendering
#define NEW_TIME 1
#define OLD_TIME 0

enum RAY_VISIBLE_MASK{
    RAY_VISIBLE_DYNAMIC = 1,
    RAY_VISIBLE_GHOST = 1<<1,
    RAY_VISIBLE_MIX = RAY_VISIBLE_DYNAMIC | RAY_VISIBLE_GHOST
};

/// @brief ///copied from pbrt3
/// @param u 
/// @return /
inline Vector2 ConcentricSampleDisk(const Vector2 &u) {
    // Map uniform random numbers to $[-1,1]^2$
    Vector2 uOffset = 2.0 * u - Vector2{1.0,1.0};

    // Handle degeneracy at the origin
    if (uOffset.x == 0 && uOffset.y == 0) return Vector2{0.,0.};

    // Apply concentric mapping to point
    Real theta, r;
    if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
        r = uOffset.x;
        theta = M_PI_4 * (uOffset.y / uOffset.x);
    } else {
        r = uOffset.y;
        theta = M_PI_2 - M_PI_4 * (uOffset.x / uOffset.y);
    }
    return r * Vector2{std::cos(theta), std::sin(theta)};
}

inline Vector3 CosineSampleHemisphere(const Vector2 &u) {
    Vector2 d = ConcentricSampleDisk(u);
    Real z = std::sqrt(std::max((Real)0, 1 - d.x * d.x - d.y * d.y));
    return Vector3(d.x, d.y, z);
}
inline Real CosineHemispherePdf(Real cosTheta) { return cosTheta * M_1_PI; }

inline Vector2 squareToUniformDiskConcentric(const Vector2 &rnd_param){
    Real r1 = 2.0 * rnd_param[0] - 1.0;
    Real r2 = 2.0 * rnd_param[1] - 1.0;

    Real phi, r;
    if(r1 == 0 && r2 == 0){
        r = phi = 0;
    }else if(r1 * r1 > r2 * r2){
        r = r1;
        phi = (M_PI/4.0) * (r2/r1);
    }else{
        r = r2;
        phi = (M_PI/2.0) - (r1/r2) * (M_PI/4.0);
    }
    Real cosPhi, sinPhi;

    cosPhi = std::cos(phi);
    sinPhi = std::sin(phi);
    return Vector2{r*cosPhi, r *sinPhi};
}

inline Vector3 sample_cos_hemisphere(const Vector2 &rnd_param, Real &pdf) {
    Real phi = c_TWOPI * rnd_param[0];
    Real tmp = sqrt(std::clamp(1 - rnd_param[1], Real(0), Real(1)));
    Vector3 ret{
        cos(phi) * tmp, sin(phi) * tmp,
        sqrt(std::clamp(rnd_param[1], Real(0), Real(1)))
    };
    pdf = ret.z * c_INVPI;
    return ret;

}

inline Vector3 UniformSampleSphere(const Vector2 &u) {
    Real z = 1 - 2 * u[0];
    Real r = sqrt(std::max((Real)0, (Real)1 - z * z));
    Real phi = 2 * c_PI * u[1];
    Vector3 ret{r * std::cos(phi), r * std::sin(phi), z};
    return normalize(ret);
}

inline Real UniformSpherePdf() { return M_1_PI * 0.25; }
