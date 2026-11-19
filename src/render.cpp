#include "render.h"
#include "intersection.h"
#include "material.h"
#include "parallel.h"
#include "pcg.h"
#include "progress_reporter.h"
#include "scene.h"
#include "integrators/path_sampling.h"
#include "integrators/portal_utils.h"
#include "integrators/portal_pt.h"
#include "integrators/portal_pt_ghost.h"
#include "integrators/incremental_sm_solid.h"
#include "integrators/path_tracing_frame.h"
#include <stdexcept>
#include <mutex>
#include <limits>

Image3 path_render(const Scene &scene)
{
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);
    constexpr int tile_size = 16;

    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    const int seed2 = static_cast<int>(0x853c49e6748fea9bULL) + 42;

    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile)
                 {

        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Spectrum radiance = make_zero_spectrum();
                int spp = scene.options.samples_per_pixel;
                for (int s = 0; s < spp; s++) {
                    uint64_t steam_id = (uint64_t(x) * h + y) * spp + s + 4096;
                    pcg32_state rng_per_sample = init_pcg32(steam_id, seed2);
                    radiance += path_tracing_frame(scene, x, y, rng_per_sample,
                                                   scene.options.frame_idx == 1 ? NEW_TIME : OLD_TIME);
                }
                img(x, y) = radiance / Real(spp);
            }
        }
        reporter.update(1); },
                 Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    return img;
}

Image3 render_fd_by_incremental(const Scene &scene)
{
    int w = scene.camera.width, h = scene.camera.height;
    Image3 result(w, h);
    std::mutex merge_mutex;
    if (scene.options.max_depth < 2)
    {
        printf("[ERROR] portal: max depth must be at least 2!\n");
        return result;
    }
    printf("[INFO] Rendering the finite difference using incremental renderer\n");
    int64_t granularity = scene.options.granularity;
    printf("[INFO] granularity = %lld\n", static_cast<long long>(granularity));
    int spp = scene.options.samples_per_pixel;
    int total_num_samples = spp * w * h;
    int64_t num_blocks = (total_num_samples + granularity - 1) / granularity;
    printf("[INFO] total blocks = %lld\n", static_cast<long long>(num_blocks));
    int maxDepth = scene.options.max_depth;
    ProgressReporter reporter(num_blocks);
    parallel_for([&](const int &blockid) {
        Image3 img(w, h);
        int num_samples = min(granularity, total_num_samples - (blockid * granularity));
        for (int s = 0; s < num_samples; ++s)
        {
            int num_pixel = min(blockid * granularity + s, (int64_t)total_num_samples) / spp;
            int x = num_pixel % w;
            int y = num_pixel / w;

            int numSamples = total_num_samples;
            std::vector<PathVertex> lightPath, sensorPath;
            lightPath.reserve(maxDepth);
            sensorPath.reserve(maxDepth);

            int seed = 0x853c49e6748fea9bULL;
            pcg32_state rng1;

            { // Path mapping version
                rng1 = init_pcg32(min(blockid * granularity + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                PathVertex vDynamic_22;
                sampleDynamicPos(vDynamic_22, rng1, scene, NEW_TIME); // without loss of generality
                pcg32_state shifted_rng1(rng1);
                Spectrum Llight_22 = prepare_dynamic_path_1_from_light(vDynamic_22, lightPath, sensorPath, rng1, scene, numSamples, NEW_TIME);
                const PathVertex &shiftedVDynamic_22 = get_vertex_at_time(vDynamic_22, OLD_TIME, scene);
                if (average(Llight_22) > 0)
                    compute_shifted_dynamic_path_1_from_light(vDynamic_22, shiftedVDynamic_22, lightPath, sensorPath, rng1, shifted_rng1, scene, numSamples, NEW_TIME, OLD_TIME, Llight_22, img);
            }
            {
                rng1 = init_pcg32(min(blockid * granularity + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                PathVertex vDynamic_11;
                sampleDynamicPos(vDynamic_11, rng1, scene, OLD_TIME);
                pcg32_state shifted_rng1(rng1);
                Spectrum Llight_11 = prepare_dynamic_path_1_from_light(vDynamic_11, lightPath, sensorPath, rng1, scene, numSamples, OLD_TIME);
                const PathVertex &shiftedVDynamic_11 = get_vertex_at_time(vDynamic_11, NEW_TIME, scene);
                if (average(Llight_11) > 0)
                    compute_shifted_dynamic_path_1_from_light(vDynamic_11, shiftedVDynamic_11, lightPath, sensorPath, rng1, shifted_rng1, scene, numSamples, OLD_TIME, NEW_TIME, Llight_11, img);
            }
            seed += 1000;
            {
                rng1 = init_pcg32(min(blockid * granularity + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();

                PathVertex vDynamic_22;
                sampleDynamicPos(vDynamic_22, rng1, scene, NEW_TIME);
                pcg32_state shifted_rng1(rng1); //copy before generating the dir
                Vector3 dynamic_dir2Light_22;
                sampleDynamicDir(vDynamic_22, dynamic_dir2Light_22, rng1, scene, NEW_TIME);

                prepare_dynamic_path_3_two_ends(vDynamic_22, dynamic_dir2Light_22, lightPath, sensorPath, rng1, scene, numSamples, NEW_TIME);
                const PathVertex &shiftedVDynamic_22 = get_vertex_at_time(vDynamic_22, OLD_TIME, scene);
                compute_shifted_dynamic_path_3_two_ends_dp(vDynamic_22, shiftedVDynamic_22, lightPath, sensorPath, rng1, shifted_rng1, scene, numSamples, NEW_TIME, OLD_TIME, img);
            }

            {
                rng1 = init_pcg32(min(blockid * granularity + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();

                PathVertex vDynamic_11;
                sampleDynamicPos(vDynamic_11, rng1, scene, OLD_TIME); // without loss of generality

                pcg32_state shifted_rng1(rng1); //copy before generating the dir
                Vector3 dynamic_dir2Light_11;
                sampleDynamicDir(vDynamic_11, dynamic_dir2Light_11, rng1, scene, OLD_TIME);

                prepare_dynamic_path_3_two_ends(vDynamic_11, dynamic_dir2Light_11, lightPath, sensorPath, rng1, scene, numSamples, OLD_TIME);
                const PathVertex &shiftedVDynamic_11 = get_vertex_at_time(vDynamic_11, NEW_TIME, scene);
                compute_shifted_dynamic_path_3_two_ends_dp(vDynamic_11, shiftedVDynamic_11, lightPath, sensorPath, rng1, shifted_rng1, scene, numSamples, OLD_TIME, NEW_TIME, img);
            }
            pcg32_state rng2;
            {
                rng2 = init_pcg32(min(blockid * (granularity - 1) + s, (int64_t)total_num_samples), seed);
                RAY_VISIBLE_MASK ray_mask = RAY_VISIBLE_MIX;
                Spectrum pt_radiance = render_pt_frame_selector_ghost(scene, x, y, ray_mask, rng2, NEW_TIME) / Real(scene.options.samples_per_pixel);
                img(x, y) = img(x, y) + pt_radiance;

                rng2 = init_pcg32(min(blockid * (granularity - 1) + s, (int64_t)total_num_samples), seed);
                pt_radiance = render_pt_frame_selector_ghost(scene, x, y, ray_mask, rng2, OLD_TIME) / Real(scene.options.samples_per_pixel);
                img(x, y) = img(x, y) - pt_radiance;
            }
            seed += 1000;
            {// 2) for s==1 and t>=2 resample the light vertex and re-generate the sensorSubapth
                rng2 = init_pcg32(min(blockid * (granularity - 1) + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                PathVertex vGhost;
                Vector3 ghost_dir2Light;
                sampleGhost(vGhost, ghost_dir2Light, rng2, scene, NEW_TIME);
                Spectrum Llight = prepare_ghost_path_1_from_light(vGhost, lightPath, sensorPath, rng2, scene, numSamples, NEW_TIME);
                if (average(Llight) > 0){
                    render_ghost_path_1_from_light(vGhost, Llight, lightPath, sensorPath, rng2, scene, numSamples, NEW_TIME, img);
                }   
            }
            {
                rng2 = init_pcg32(min(blockid * (granularity - 1) + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                
                PathVertex vGhost;
                Vector3 ghost_dir2Light;
                sampleGhost(vGhost, ghost_dir2Light, rng2, scene, OLD_TIME);
                Spectrum Llight = prepare_ghost_path_1_from_light(vGhost, lightPath, sensorPath, rng2, scene, numSamples, OLD_TIME);
                if (average(Llight) > 0){
                    Llight *= -1.0;
                    render_ghost_path_1_from_light(vGhost, Llight, lightPath, sensorPath, rng2, scene, numSamples, OLD_TIME, img);
                }   
            }
            seed += 1000;
            {
                rng2 = init_pcg32(min(blockid * (granularity - 1) + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                
                PathVertex vGhost;
                Vector3 ghost_dir2Light;
                sampleGhost(vGhost, ghost_dir2Light, rng2, scene, NEW_TIME);
                prepare_ghost_path_3_two_ends(vGhost, ghost_dir2Light, lightPath, sensorPath, rng2, scene, numSamples, NEW_TIME);
                render_ghost_path_3_two_ends(vGhost, ghost_dir2Light, lightPath, sensorPath, rng2, scene, numSamples, NEW_TIME, img, false);
            } 
            {
                rng2 = init_pcg32(min(blockid * (granularity - 1) + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                
                PathVertex vGhost;
                Vector3 ghost_dir2Light;
                sampleGhost(vGhost, ghost_dir2Light, rng2, scene, OLD_TIME);
                prepare_ghost_path_3_two_ends(vGhost, ghost_dir2Light, lightPath, sensorPath, rng2, scene, numSamples, OLD_TIME);
                render_ghost_path_3_two_ends(vGhost, ghost_dir2Light, lightPath, sensorPath, rng2, scene, numSamples, OLD_TIME, img, true);
            } 
            seed += 1000;
            pcg32_state rng3;
            {
                rng3 = init_pcg32(min(blockid * granularity + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                PathVertex vDynamic; Vector3 dynamic_dir2Light;
                sampleDynamic(vDynamic, dynamic_dir2Light, rng3, scene, NEW_TIME); // without loss of generality
                Vector2 pRaster2;
                Spectrum We = prepare_dynamic_path_2_from_sensor(vDynamic, lightPath, sensorPath, rng3, scene, numSamples, NEW_TIME, pRaster2);
                if (average(We) > 0){
                    render_dynamic_path_2_from_sensor(lightPath, sensorPath, rng3, scene, numSamples, NEW_TIME, We, pRaster2, img);
                }
            }
            {//// GHOST_FROM_SENSOR; OLD_TIME
                rng3 = init_pcg32(min(blockid * granularity + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                PathVertex vGhost;
                Vector3 ghost_dir2Light;
                sampleGhost(vGhost, ghost_dir2Light, rng3, scene, OLD_TIME);
                Vector2 pRaster2{-1, -1};
                Spectrum We = prepare_ghost_path_2_from_sensor(vGhost, lightPath, sensorPath, rng3, scene, numSamples, OLD_TIME, pRaster2);
                if (average(We) > 0){
                    We *= -1.0;
                    render_ghost_path2_from_sensor(vGhost, We, pRaster2, lightPath, sensorPath, rng3, scene, numSamples, OLD_TIME, img);
                }
            }
            seed += 1000;
            {//DYNAMIC_FROM_SENSOR OLD_TIME
                rng3 = init_pcg32(min(blockid * granularity + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                PathVertex vDynamic; Vector3 dynamic_dir2Light;
                sampleDynamic(vDynamic, dynamic_dir2Light, rng3, scene, OLD_TIME); // without loss of generality
                Vector2 pRaster2;
                Spectrum We = prepare_dynamic_path_2_from_sensor(vDynamic, lightPath, sensorPath, rng3, scene, numSamples, OLD_TIME, pRaster2);
                if (average(We) > 0){
                    We *= -1.0;
                    render_dynamic_path_2_from_sensor(lightPath, sensorPath, rng3, scene, numSamples, OLD_TIME, We, pRaster2, img);
                }
            }
            {
                rng3 = init_pcg32(min(blockid * granularity + s, (int64_t)total_num_samples), seed);
                lightPath.clear();
                sensorPath.clear();
                PathVertex vGhost;
                Vector3 ghost_dir2Light;
                sampleGhost(vGhost, ghost_dir2Light, rng3, scene, NEW_TIME);
                Vector2 pRaster2{-1, -1};
                Spectrum We = prepare_ghost_path_2_from_sensor(vGhost, lightPath, sensorPath, rng3, scene, numSamples, NEW_TIME, pRaster2);
                if (average(We) > 0){
                    render_ghost_path2_from_sensor(vGhost, We, pRaster2, lightPath, sensorPath, rng3, scene, numSamples, NEW_TIME, img);
                }
            }

        }
        {
            std::lock_guard<std::mutex> lock(merge_mutex);
            for (size_t pixel = 0; pixel < result.data.size(); ++pixel)
                result.data[pixel] += img.data[pixel];
        }
        reporter.update(1);
    },
                 num_blocks);
    reporter.done();
    return result;
}

Image3 render(const Scene &scene)
{
    if (!scene.media.empty())
        throw std::runtime_error("This release supports surface transport only.");
    if (scene.camera.width <= 0 || scene.camera.height <= 0 ||
        static_cast<int64_t>(scene.camera.width) * scene.camera.height >
            std::numeric_limits<int>::max() / max(1, scene.options.samples_per_pixel))
        throw std::runtime_error("Invalid resolution or sample budget exceeds 32-bit indexing.");
    if (scene.options.samples_per_pixel <= 0 || scene.options.granularity <= 0)
        throw std::runtime_error("Sample count and granularity must be positive.");
    if (scene.options.max_depth < 2)
        throw std::runtime_error("maxDepth must be at least 2.");
    if (scene.options.frame_idx != 1 && scene.options.frame_idx != 2)
        throw std::runtime_error("frame_idx must be 1 (new) or 2 (old).");
    if (scene.options.integrator == Integrator::Path)
        return path_render(scene);
    if (scene.dynamics.empty())
        return Image3(scene.camera.width, scene.camera.height);
    return render_fd_by_incremental(scene);
}
