#include "image.h"
#include <cmath>
#include <iomanip>
#include <stdexcept>

int main(int argc, char *argv[]) {
    try {
        if (argc < 3 || argc > 4)
            throw std::runtime_error("Usage: image_compare new.exr old.exr [difference.exr]");
        Image3 first = imread3(argv[1]);
        Image3 second = imread3(argv[2]);
        if (first.width != second.width || first.height != second.height)
            throw std::runtime_error("Image dimensions differ.");
        Image3 difference = first - second;
        double squared_error = 0;
        double maximum = 0;
        double minimum_value = infinity<double>();
        double maximum_value = -infinity<double>();
        for (size_t pixel = 0; pixel < difference.data.size(); ++pixel) {
            for (int channel = 0; channel < 3; ++channel) {
                double value = difference.data[pixel][channel];
                if (!std::isfinite(value)) throw std::runtime_error("Non-finite image value.");
                squared_error += value * value;
                maximum = std::max(maximum, std::abs(value));
                minimum_value = std::min(minimum_value, value);
                maximum_value = std::max(maximum_value, value);
            }
        }
        std::cout << std::setprecision(12)
                  << "rmse=" << std::sqrt(squared_error / (3 * difference.data.size()))
                  << " max_abs=" << maximum << " min=" << minimum_value
                  << " max=" << maximum_value << "\n";
        if (argc == 4) imwrite(argv[3], difference);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
