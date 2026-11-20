<h1 align="center"> Residual Path Integrals for Re-rendering</h1>

<p align="center">
  <strong>🏆 EGSR 2024 Best Paper Award</strong><br>
  <a href="https://bingxu.tech/papers/residual_path_integral_bing_final.pdf">Paper</a> ·
  <a href="https://arxiv.org/abs/2406.16302">arXiv</a> ·
  <a href="https://bingxu.tech/projects/2024_residual/page.html">Project page</a>
</p>

---

This is a release of the CPU research code used to generate the results in
[Residual Path Integrals for Re-rendering](https://bingxu.tech/projects/2024_residual/page.html),
by Bing Xu, Tzu-Mao Li, Iliyan Georgiev, Trevor Hedstrom, and Ravi Ramamoorthi (2024).

The renderer estimates a **signed, linear RGB residual**
`D = I_new - I_old` between two scene states. Given an accurately rendered old image,
the updated image is `I_old + D`. Negative residual values are valid; do not clamp,
take the absolute value, or tone-map before adding the images.

Included: the residual renderer, a path tracer, the original Cornell box and
staircase scenes, and a small example added for this release. Generated images
are not included.


Notice that 
- This code was from 2024 for the paper results. The paper's timings apply to the original implementation. I guess now at 2026, we encourage using coding agents to optimize this code.

- I spent a lot of time debugging to ensure unbiasedness. Some files and helper routines remain from those intermediate debugging steps and may not be needed for the final residual rendering.

- It currently contains two example scenes, I'll look for and upload the rest of the scenes we created back then soon.

## 🏎️ Download

Scene meshes (`.obj`) and textures (`.jpg`) are stored with Git LFS. Install
[Git LFS](https://git-lfs.com/) before cloning:

```sh
git lfs install
git clone https://github.com/bing-xu-graphics/residual_path_integrals.git
cd residual_path_integrals
git lfs pull
```

For an existing checkout, run `git lfs install` and `git lfs pull` from the
repository root. Use a Git clone rather than a source ZIP, which may contain
LFS pointers instead of the assets. Scene XML, material definitions, and license
files remain in regular Git. Rendered outputs are not included.

## 🏎️ Build

Requirements: a C++17 compiler, CMake 3.20 or newer, and **Embree 3** with its runtime
dependencies. Embree 4 may not API-compatible with this source.
Download an [Embree 3 release](https://github.com/RenderKit/embree/releases/tag/v3.13.5)
or build it from source.
Third-party image and XML utilities are already included in `src/3rdparty/`.
No GPU, Python package, profiler, or graphics-window library is required to render.

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEMBREE_ROOT=/path/to/embree3
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

On Windows, run CMake with an installed Visual Studio C++ toolchain and Windows SDK:

```powershell
cmake -S . -B build -A x64 -DEMBREE_ROOT=C:/libraries/embree3
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

Executables are `build/glacier` and `build/image_compare` with a single-configuration
generator, or `build/Release/glacier.exe` and `build/Release/image_compare.exe` with
Visual Studio. Commands below use the single-configuration paths; substitute the
Windows paths as needed.

If CMake cannot find Embree, set `EMBREE_INCLUDE_DIR` (containing `embree3/`) and
`EMBREE_LIBRARY` (the full library filename). Windows DLLs are copied from the
library directory's sibling `bin/`; override `EMBREE_RUNTIME_DIR` if necessary.
On Linux/macOS, Embree and TBB must be on the runtime library search path.
Use a native Embree build on Apple Silicon. Only Windows x64 has been tested for
this release.

## 🏎️ Render the difference image

```sh
./build/glacier scenes/cornell-box/scene.xml --spp 16 -o output/residual.exr
```

The output is a signed **32-bit floating-point EXR**. Parent output directories
are created automatically. For a fast preview:

```sh
./build/glacier scenes/cornell-box/scene.xml --resolution 64 --spp 4 -t 1 -o output/preview.exr
```

I recommend using an EXR viewer that supports negative values, such as
[tev](https://github.com/Tom94/tev). A viewer that clamps values to zero will hide
the negative part of the residual. 


## 🏎️ Run the paper scenes

The scene files, meshes, and textures are included:

| Scene | Residual input | Original resolution | Residual sample budget |
| --- | --- | --- | --- |
| Cornell box (historical `001_cbox`) | `scenes/paper-cornell-box/residual.xml` | 512 × 512 | 4 |
| Wooden staircase (historical `002_staircase`) | `scenes/paper-staircase/residual.xml` | 720 × 1280 | 4 |
| New minimal example (not a paper scene) | `scenes/cornell-box/scene.xml` | 256 × 256 | 16 |

Run the original residual configurations without changing their resolution or
sample budget:

```sh
./build/glacier scenes/paper-cornell-box/residual.xml -t 4 -o output/paper-cornell-box/residual.exr
./build/glacier scenes/paper-staircase/residual.xml -t 4 -o output/paper-staircase/residual.exr
```

For a quick test, append `--resolution 32 --spp 1`. Omit these options to use the
original settings. The active scene settings are unchanged. Their sampler declarations may
print a warning: the renderer uses independent samples, as the original code did.
Output is always signed EXR, regardless of the XML film format.

To render the residual and record the command and elapsed time, use the Python
script (standard library only):

```sh
python tools/render_scene.py --bin-dir build --scene paper-cornell-box
python tools/render_scene.py --bin-dir build --scene paper-staircase
```

On Windows use `--bin-dir build/Release`. The script uses the XML residual sample
count and four threads. Change these with `--spp`, `--resolution`, and `--threads`.
The residual image, command, and elapsed time are saved in `output/<scene>/`.

Each paper scene also preserves the historical `004_frame_1_gt.xml` as
`reference-old.xml` and `006_frame_2_gt.xml` as `reference-new.xml`;
`residual.xml` comes from `002_incremental.xml`. Their original PT budgets are
10,240 samples for Cornell box and 2,048 for staircase:

```sh
./build/glacier scenes/paper-cornell-box/reference-old.xml -t 4 -o output/paper-cornell-box/historical-old.exr
./build/glacier scenes/paper-cornell-box/reference-new.xml -t 4 -o output/paper-cornell-box/historical-new.exr
./build/glacier scenes/paper-staircase/reference-old.xml -t 4 -o output/paper-staircase/historical-old.exr
./build/glacier scenes/paper-staircase/reference-new.xml -t 4 -o output/paper-staircase/historical-new.exr
```

These reference XML files may describe a different edit from `residual.xml`.
We have not reproduced the paper figures with the cleaned executable.

## 🏎️ Path-tracing baseline

Render a conventional path-traced image:

```sh
./build/glacier scenes/cornell-box/scene.xml --integrator path --frame new --spp 256 -o output/new.exr
```

Use `--frame old` to render the previous scene state instead.

## 🏎️ Scene changes and controls

A scene describes the **new** geometry. A shape's optional `dynamic` transform
maps its new world-space geometry to its old state:

```xml
<shape type="obj">
    <string name="filename" value="meshes/moving-box.obj"/>
    <transform name="toWorld">
        <scale value="0.5"/>
        <translate x="0.6" y="0.65" z="2"/>
    </transform>
    <ref id="white"/>
    <dynamic>
        <transform name="toWorld">
            <translate x="-1" y="0" z="0"/>
        </transform>
    </dynamic>
</shape>
```

For meshes, the old transform is `toWorldOld * toWorld`. Unmarked shapes are
static. Keep the camera and lights fixed, and use rigid object transformations.
Deformation, topology changes, camera/light animation, and material edits have
not been tested in this release. The renderer rejects participating-media scenes.

| Control | Meaning |
| --- | --- |
| `--integrator residual\|path` | Override the XML integrator |
| `--frame old\|new` | Select the PT state; residual always uses both |
| `--spp N` | Positive sample-budget override |
| `--resolution N` | Output width, with height scaled to preserve aspect ratio |
| `-t N` | Worker count; use 1 for deterministic diagnostic comparisons |
| `-o path.exr` | Output filename; default is `residual.exr` |
| `--help` | Print usage |

XML integrator types are `residual` and `path`; `fd_incremental` is retained as
an alias. Integrator options are `maxDepth` (at least 2; default 5), `rrDepth`,
`granularity` (positive residual samples per work block), and `frame_idx`
(1 = new, 2 = old, for PT). Pixel sampling is independent. Sample seeds are fixed
in `src/render.cpp`; record any seed changes when doing new experiments.

## 🏎️ Code layout

Back in 2024, I spent a lot of time debugging this code to identify and remove
sources of bias. Some files and helper routines remain from those intermediate
debugging steps and are not needed for the final residual rendering.

| File | Responsibility |
| --- | --- |
| `src/render.cpp` | Rendering loops, sample counts, random streams, and image accumulation |
| `src/integrators/path_tracing_frame.h` | Conventional PT at a selected scene state |
| `src/integrators/path_sampling.h` | Shared random walks, direct illumination, and PDFs |
| `src/integrators/incremental_basepath.h` | Dynamic-vertex and ghost-edge base-path generation |
| `src/integrators/incremental_sm_solid.h` | Active dynamic-path mappings and paired contributions |
| `src/integrators/incremental_shiftmapping_utils.h` | Mapping, reconnection, and PDF helpers |
| `src/integrators/portal_pt*.h` | Dynamic/ghost contribution evaluation used by the residual estimator |
| `src/integrators/portal_utils*.h` | Sampling distributions and multiple-importance weights |
| `src/dynamic.h`, `src/shapes/dynamic.inl`, `src/shapes/ghost.inl` | Two-state geometry and ghost registration |
| `src/parsers/parse_scene.cpp`, `src/intersection.cpp` | XML scene loading and state-aware ray queries |

This paper was from 2024 for the paper results. I guess now at 2026, we encourage using coding agents to optimize this code. The paper's timings apply to the original implementation.

The residual entry point is `render_fd_by_incremental()`. Sampling, mapping
Jacobians, rejection rules, and MIS weights must agree when adding a technique.
Some helper headers also contain alternative formulations that are not exposed
as rendering modes.

Residual splats accumulate into a private image per active work block, then merge
under a lock. This removes concurrent writes to the same pixel, at a memory cost
of roughly `24 * width * height` bytes per active worker plus the output image.
Thread scheduling can change the floating-point summation order. Assertions
remain enabled in release builds.

## 🏎️ Validation
CTest runs seven C++ tests. To run the rendering checks, use Python 3
with its standard library:

```sh
python tools/validate.py --bin-dir build --output output/validation
```

On Windows use `--bin-dir build/Release`. The script saves EXRs, generated test
scenes, every command and its output (`commands.json`), and `metrics.json`.
It checks single-thread repeatability, one-versus-four-thread agreement, an exact
zero when no dynamic geometry exists, finite output, and both residual signs.
It also measures the residual for an identity transform and compares a moving
scene against a higher-sample PT difference.


## 🏎️ License for the assets etc

The paper scenes keep their original active settings and asset filenames.
The top-level scene folders and the three XML filenames were renamed, and the
commented-out moving-chair blocks were removed from the staircase XML files.
The original repository credits its downloaded
scenes to the Mitsuba scene collection, Wenzel Jakob, the Cornell Program of
Computer Graphics, Marko Dabrovic, Eric Veach, Jonas Pilo, and Bernhard Vogl.
The Wooden Staircase is by **Wig42**, distributed via
[Benedikt Bitterli's rendering resources](https://benedikt-bitterli.me/resources/)
under [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/); its original
[scene license](scenes/paper-staircase/LICENSE.txt) is retained. Third-party scene
assets are not relicensed by the code's MIT license.

Original contributions and the inherited
[lajolla renderer](https://github.com/BachiLi/lajolla_public) are covered by
[LICENSE](LICENSE), retaining professor Tzu-Mao Li's upstream MIT attribution.
Third-party components retain their own notices in their source:
pugixml, stb_image, tinyexr (including embedded OpenEXR notices), miniz, and tinyply.
The renderer uses PCG random-number generation and draws on pbrt, Mitsuba, and
SmallVCM. The pbrt-derived threading code retains its
[BSD notice](licenses/pbrt.txt), retrieved from the
[upstream license](https://github.com/mmp/pbrt-v3/blob/master/LICENSE.txt).
Embree and its runtime dependencies are installed separately.


## 🏎️ Citation
When using this work, please cite: 

```bibtex
@article{xu2024residual,
  title   = {Residual Path Integrals for Re-rendering},
  author  = {Xu, Bing and Li, Tzu-Mao and Georgiev, Iliyan and
             Hedstrom, Trevor and Ramamoorthi, Ravi},
  journal = {Computer Graphics Forum},
  volume  = {43},
  number  = {4},
  year    = {2024},
  doi     = {10.1111/cgf.15152}
}
```
