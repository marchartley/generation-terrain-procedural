The OpenGL lifecycle smoke test uses the real viewer and requires a desktop OpenGL 4.5 compatibility context. Run it from the repository root so terrain definitions, shaders, and assets resolve correctly.

```sh
mkdir -p /tmp/terrain-smoke-build
/opt/Qt/6.12.0/gcc_64/bin/qmake tests/opengl-lifecycle.pro -o /tmp/terrain-smoke-build/Makefile
make -C /tmp/terrain-smoke-build -j4
/tmp/terrain-smoke-build/opengl-lifecycle-smoke
```

Adjust the qmake path for your Qt installation. A successful run prints `PASS` and exits with status 0. The test checks construction without a context, shader compilation and uniform updates from a GUI callback, restoration of the current context, independent buffers for copied meshes, shader reload, water level changes, all four terrain representations, textured framebuffer capture, and GPU resource rebuilding without replacing CPU terrain, and cleanup when the viewer is destroyed. It generates its heightmap and screenshot in a temporary directory.

Pass a filename as the first argument to retain the screenshot, for example `/tmp/terrain-smoke-build/opengl-lifecycle-smoke /tmp/terrain-frame.png`. The test enables `TERRAIN_GL_DEBUG=1` for driver diagnostics; the same variable can be used when launching the interactive application.

The gizmo regression creates and destroys renderers concurrently on eight worker threads before OpenGL is available, then renders the prepared objects, their centers, translation axes and rotation rings. It also exercises the first visible main mouse gizmo frame and mouse press/release events. Before the fix, this first-render test reproduces a null `Shader::this` in `Shader::ensureResources()` without requiring a race.

The smoke test also exercises all four “convert to all” actions on a sloping,
fractional-height terrain, checks the resulting elevations, and compares the
voxel, layered and implicit framebuffer images. Marching Cubes uses the same
unfiltered sampling path for those three views; raw layer boxes remain available
when Marching Cubes is disabled.

The terrain conversion regression runs without a display or OpenGL context:

```sh
/opt/Qt/6.12.0/gcc_64/bin/qmake tests/terrain-conversions.pro -o build/Desktop_Qt_6_12_0_Debug/Makefile.conversions CONFIG+=debug CONFIG-=release
make -C build/Desktop_Qt_6_12_0_Debug -f Makefile.conversions -j4
OMP_NUM_THREADS=4 ./build/Desktop_Qt_6_12_0_Debug/terrain-conversions
```

Use an existing compatible build directory to reuse the application's object
files, or create a separate build directory and adjust these paths. The checks
cover zero/unit/fractional heights, slopes, single-sample columns, caves, material
independence, box-kernel volume conservation, coordinate transforms, requested
sampling density, implicit cache resizing, terrain above height 40, replacement
of old layer contents, and serialization of converted sampled fields.

Conversion conventions:

- Voxel fields have samples at integer coordinates, positive values inside and
  zero on the surface. Heights and layer thicknesses use the same terrain units.
- `getVoxelized(dimensions, scale)` takes an extent in terrain units and a scale
  in samples per unit. Elevations are never normalized to the destination depth.
- Default layer reconstruction preserves vertical interfaces without filtering.
  `voxelize(height, kernelSize)` with a positive kernel size explicitly requests
  the volume-fraction box reconstruction described in *Arches*.
- Conversion to a heightmap extracts the highest surface and necessarily loses
  caves, overhangs, and materials. A voxel field has no material labels; converting
  it to layers supplies bedrock for solid intervals. Features smaller than the
  sampling interval cannot always survive rasterization.


Rendering performance can be measured with the same deterministic 64×64×32
terrain in each view:

```sh
/opt/Qt/6.12.0/gcc_64/bin/qmake tests/terrain-render-benchmark.pro -o build/Desktop_Qt_6_12_0_Debug/Makefile.benchmark CONFIG+=debug CONFIG-=release
make -C build/Desktop_Qt_6_12_0_Debug -f Makefile.benchmark -j4
OMP_NUM_THREADS=4 ./build/Desktop_Qt_6_12_0_Debug/terrain-render-benchmark /tmp/terrain-render
```

Each case targets 750 ms of warm-up (with a minimum of three frames and a
256-frame cap), then reports the median and minimum of eleven complete framebuffer
captures. These are CPU/GPU frame latencies including readback, rather than GPU
shader-only timings. Vsync is disabled for the benchmark. Keep the Qt build,
window size, desktop session, GPU and OpenMP settings identical when comparing
runs. An optional filename prefix retains a screenshot for each case.

The cases cover heightmaps, voxel/layer/implicit Marching Cubes, heightmap and
voxel ambient occlusion, raw layer boxes, raw voxels, and soft shadows. Set
`TERRAIN_RENDER_TIMINGS=1` when launching the app or benchmark to print CPU
submission timings for interfaces, uniforms, picking and terrain texture work.

The active terrain pipelines are:

| View | Geometry stage | Fragment stage |
| --- | --- | --- |
| Heightmap | `grid.geom`, adaptive quarter-cell detail | `MarchingCubes.frag` |
| Voxels, implicit, smooth layers | `MarchingCubes.geom` | `MarchingCubes.frag` |
| Raw layers | `layer_based.geom`, cumulative layer heights | `MarchingCubes.frag` |

Heightmaps select 1, 2 or 4 subdivisions from projected bilinear curvature, using
0.35 pixels as the refinement threshold and retaining the previous maximum
detail. The shader's `adaptiveTessellation` uniform can force fixed maximum detail
for regression comparisons. Terrain normals use a six-read central-difference
stencil, and soft shadows use twenty stratified rays with bounded voxel traversal.
Raw layer ambient occlusion now skips all sampling when its strength is zero.

Shader programs cache link results, uniform locations and values. Texture storage
is reused for same-size edits; unchanged volume data is rebound without another
upload. Layer sampling and cumulative-height caches compare actual intervals,
including direct edits that do not update history counters. Heightmap texture
caching accounts for elevations, biome assignments and names, texture palettes,
OMP settings, and GPU rebuilds.

The desktop regression checks negative texture values, same-size updates,
resizing, uniform changes and program relinking, compute shader uniforms, linking
of terrain/mesh/rock/gizmo/particle shader families, adaptive versus fixed-detail
images, and context loss/recreation. The CPU suite also checks direct layer edits
and cumulative-height cache invalidation.

[Measured before/after results and validation scope](render-performance.md).
