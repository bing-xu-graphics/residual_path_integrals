#include "parsers/parse_scene.h"
#include "parallel.h"
#include "image.h"
#include "render.h"
#include <embree3/rtcore.h>
#include <thread>
#include <stdexcept>

int main(int argc, char *argv[]) {
    RTCDevice device = nullptr;
    bool initialized = false;
    try {
        int threads = max(1, static_cast<int>(std::thread::hardware_concurrency()));
        int spp = 0;
        int resolution = 0;
        int frame = 0;
        std::string mode;
        fs::path input;
        fs::path output = "residual.exr";
        for (int argument = 1; argument < argc; ++argument) {
            std::string option = argv[argument];
            if (option == "--help" || option == "-h") {
                std::cout << "Usage: glacier scene.xml [-o image.exr] [-t threads]\n"
                             "  [--integrator residual|path] [--frame old|new]\n"
                             "  [--spp samples] [--resolution width_pixels]\n";
                return 0;
            }
            if (option == "-o" || option == "-t" || option == "--spp" ||
                option == "--resolution" || option == "--integrator" || option == "--frame") {
                if (++argument >= argc) throw std::runtime_error("Missing value for " + option);
                std::string value = argv[argument];
                if (option == "-o") output = value;
                else if (option == "--integrator") {
                    if (value != "residual" && value != "path")
                        throw std::runtime_error("Integrator must be residual or path.");
                    mode = value;
                } else if (option == "--frame") {
                    if (value != "old" && value != "new")
                        throw std::runtime_error("Frame must be old or new.");
                    frame = value == "new" ? 1 : 2;
                } else {
                    size_t consumed;
                    int number = std::stoi(value, &consumed);
                    if (consumed != value.size() || number <= 0)
                        throw std::runtime_error("Expected a positive integer for " + option);
                    if (option == "-t") threads = number;
                    if (option == "--spp") spp = number;
                    if (option == "--resolution") resolution = number;
                }
            } else if (!option.empty() && option[0] == '-') {
                throw std::runtime_error("Unknown option: " + option);
            } else {
                if (!input.empty()) throw std::runtime_error("Specify exactly one scene.");
                input = option;
            }
        }
        if (input.empty()) throw std::runtime_error("Specify a scene XML; use --help for usage.");
        if (output.extension() != ".exr") throw std::runtime_error("Output must have an .exr extension.");
        input = fs::absolute(input);
        output = fs::absolute(output);
        device = rtcNewDevice(nullptr);
        if (!device) throw std::runtime_error("Could not initialize Embree.");
        {
            Scene scene = parse_scene(input, device);
            if (!mode.empty()) scene.options.integrator = mode == "path" ? Integrator::Path : Integrator::FD_INCREMENTAL;
            if (spp) scene.options.samples_per_pixel = spp;
            if (frame) scene.options.frame_idx = frame;
            if (resolution) {
                Real scale = Real(resolution) / scene.camera.width;
                Real scaled_height = std::round(scene.camera.height * scale);
                if (scaled_height > std::numeric_limits<int>::max())
                    throw std::runtime_error("Requested resolution is too large.");
                Real fov = 360.0 / c_PI * std::atan(scene.camera.width /
                                                  (2.0 * scene.camera.image_plane_distance));
                scene.camera = Camera(scene.camera.cam_to_world, fov, resolution,
                                      max(1, int(scaled_height)), scene.camera.cam_normal,
                                      scene.camera.filter, scene.camera.medium_id,
                                      scene.camera.focalDistance);
            }
            parallel_init(threads);
            initialized = true;
            Image3 result = render(scene);
            for (const auto &pixel : result.data)
                if (!std::isfinite(pixel.x) || !std::isfinite(pixel.y) || !std::isfinite(pixel.z))
                    throw std::runtime_error("Render contains non-finite pixels.");
            fs::create_directories(output.parent_path());
            imwrite(output, result);
            std::cout << "Image written to " << output << "\n";
            parallel_cleanup();
            initialized = false;
        }
        rtcReleaseDevice(device);
        return 0;
    } catch (const std::exception &error) {
        if (initialized) parallel_cleanup();
        if (device) rtcReleaseDevice(device);
        std::cerr << error.what() << "\n";
        return 1;
    }
}

