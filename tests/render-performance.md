# Terrain rendering performance

Measured on 2026-11-09 with Qt 6.12.0, the Debug configuration, OpenGL 4.6,
an NVIDIA GeForce GTX 1650 Ti and four OpenMP threads. Both builds used the same
deterministic 64×64×32 scene, assets, 900×700 application window and benchmark source.
Each case targeted 750 ms of warm-up (minimum three, maximum 256 frames), then
captured eleven measured frames. Times include
CPU submission, GPU rendering and framebuffer readback; they are not GPU-only timings.

The baseline is commit `4cc64b0`, built separately from an archive of that revision.
Its image dictionaries, textures and scene assets were linked from the unchanged
working-tree assets. The optimized run used the current changes. Startup and shader
compilation are outside the measured samples. Results will vary with scene and hardware.

| Rendering path | Before median (ms) | After median (ms) | Speedup |
| --- | ---: | ---: | ---: |
| Heightmap | 9.04 | 1.83 | 4.9× |
| Voxels / Marching Cubes | 13.47 | 5.64 | 2.4× |
| Layers / Marching Cubes | 227.81 | 8.82 | 25.8× |
| Implicit terrain | 10.05 | 5.18 | 1.9× |
| Heightmap with AO | 39.49 | 5.62 | 7.0× |
| Voxels with AO | 21.33 | 16.77 | 1.3× |
| Raw layer boxes | 81.28 | 5.87 | 13.8× |
| Raw voxels | 11.61 | 9.15 | 1.3× |
| Voxels with soft shadows | 14.89 | 6.62 | 2.2× |

The largest improvement comes from removing repeated layer reconstruction and
material-set copies from rendering, and computing terrain bounds once per frame.
Other changes cache uniform state and texture uploads, reuse GPU storage, cache
heightmap/biome data and layer prefix heights, reduce volume corner/normal sampling,
adapt heightmap tessellation to projected curvature, and bound soft-shadow traversal.

Rendering correctness fixes initialize shader outputs and material values, isolate
the comparison mask from terrain samples, normalize lighting normals, and use the
existing albedo atlas for surface detail. Particle fragments outside the visible
circle are discarded instead of writing depth with zero opacity.

See [README.md](README.md) for the benchmark and regression commands. The benchmark
can save a screenshot for every rendering path; comparisons should include both
frame times and appearance. Experimental raymarching shaders are not part of the
active terrain renderer measured here.
