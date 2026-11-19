#include "../image.h"
#include <cstdio>

int main(int argc, char *argv[]) {
    Image3 img(32, 24);
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            img(x, y) = Vector3{
                3 * (img.width * y + x) / Real(1024),
                5 * (img.width * y + x) / Real(1024),
                7 * (img.width * y + x) / Real(1024)};
        }
    }

    // round trip test
    imwrite("test.exr", img);
    Image3 rimg = imread3("test.exr");
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            Vector3 target = Vector3{
                3 * (img.width * y + x) / Real(1024),
                5 * (img.width * y + x) / Real(1024),
                7 * (img.width * y + x) / Real(1024)};
            Vector3 diff = rimg(x, y) - target;
            if (length(diff) > Real(1e-2)) {
                printf("FAIL\n");
                return 1;                
            }
        }
    }

    std::remove("test.exr");

    Image3 residual(2, 1);
    residual(0) = Vector3{-123.456789, 0.000123456789, 123.456789};
    residual(1) = Vector3{-0.25, 0.0, 0.25};
    imwrite("signed-residual.exr", residual);
    Image3 restored = imread3("signed-residual.exr");
    std::remove("signed-residual.exr");
    for (size_t pixel = 0; pixel < residual.data.size(); ++pixel) {
        if (length(residual.data[pixel] - restored.data[pixel]) > 1e-5) {
            printf("FAIL: signed full-float EXR round trip\n");
            return 1;
        }
    }

    printf("SUCCESS\n");
    return 0;
}
