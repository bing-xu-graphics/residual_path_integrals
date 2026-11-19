#pragma once

#include "vector.h"
#include "spectrum.h"

#include <string>
#include <cstring>
#include <vector>

/// A N-channel image stored in a contiguous vector
/// The storage format is HWC -- outer dimension is height
/// then width, then channels.
template <typename T>
struct Image
{
    Image() {}
    Image(int w, int h) : width(w), height(h)
    {
        data.resize(w * h);
        memset(data.data(), 0, sizeof(T) * data.size());
    }

    T &operator()(int x)
    {
        return data[x];
    }

    const T &operator()(int x) const
    {
        return data[x];
    }

    T &operator()(int x, int y)
    {
        return data[y * width + x];
    }

    const T &operator()(int x, int y) const
    {
        return data[y * width + x];
    }

    int width;
    int height;
    std::vector<T> data;
};

template <typename T>
inline Image<T> operator+(const Image<T> &img0, const Image<T> &img1)
{
    assert(img0.width == img1.width && img0.height == img1.height);
    Image<T> res(img0.width, img0.height);
    for (int y = 0; y < img0.height; y++)
    {
        for (int x = 0; x < img0.width; x++)
        {
            res(x, y) = img0(x, y) + img1(x, y);
        }
    }
    return res;
}

template <typename T>
inline Image<T> operator-(const Image<T> &img0, const Image<T> &img1)
{
    assert(img0.width == img1.width && img0.height == img1.height);
    Image<T> res(img0.width, img0.height);
    for (int y = 0; y < img0.height; y++)
    {
        for (int x = 0; x < img0.width; x++)
        {
            res(x, y) = img0(x, y) - img1(x, y);
        }
    }
    return res;
}

template <typename T>
inline Image<T> operator*(const Image<T> &img0, const Real &s)
{
    Image<T> res(img0.width, img0.height);
    for (int y = 0; y < img0.height; y++)
    {
        for (int x = 0; x < img0.width; x++)
        {
            res(x, y) = img0(x, y)* s;
        }
    }
    return res;
}

template <typename T>
inline Image<T> operator*(const Real &s, const Image<T> &img0)
{
    return img0 *s;
}

using Image1 = Image<Real>;
using Image3 = Image<Vector3>;

/// Read from an 1 channel image. If the image is not actually
/// single channel, the first channel is used.
/// Supported formats: JPG, PNG, TGA, BMP, PSD, GIF, HDR, PIC
Image1 imread1(const fs::path &filename);
/// Read from a 3 channels image.
/// If the image only has 1 channel, we set all 3 channels to the same color.
/// If the image has more than 3 channels, we truncate it to 3.
/// Undefined behavior if the image has 2 channels (does that even happen?)
/// Supported formats: JPG, PNG, TGA, BMP, PSD, GIF, HDR, PIC
Image3 imread3(const fs::path &filename);

/// Save an image to a file.
/// Supported formats: PFM & exr
void imwrite(const fs::path &filename, const Image3 &image);

inline Image3 to_image3(const Image1 &img)
{
    Image3 out(img.width, img.height);
    std::transform(img.data.cbegin(), img.data.cend(), out.data.begin(),
                   [](Real v)
                   { return Vector3(v, v, v); });
    return out;
}

inline Image1 to_image1(const Image3 &img)
{
    Image1 out(img.width, img.height);
    std::transform(img.data.cbegin(), img.data.cend(), out.data.begin(),
                   [](const Vector3 &v)
                   { return average(v); });
    return out;
}

inline void splatSample(Image3 &img, Image1 &spps, const Vector2 &pRaster2, Spectrum v)
{
    int xx = floor(pRaster2.x);
    int yy = floor(pRaster2.y);
    img(xx, yy) = img(xx, yy) + v;
}

inline void splatSample(Image3 &img, const Vector2 &pRaster2, Spectrum v)
{
    int xx = floor(pRaster2.x);
    int yy = floor(pRaster2.y);
    assert(average(v) == average(v));
    assert(!isinf(average(v)));
    img(xx, yy) = img(xx, yy) + v;
}
