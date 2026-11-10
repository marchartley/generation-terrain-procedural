#ifndef TERRAINSHADOWMAP_H
#define TERRAINSHADOWMAP_H
#include "Graphics/Mesh.h"
#include <array>

// Six positional lights, using the terrain's existing geometry program for depth.
// Cached independently of the viewer camera. Resource leases survive context resets.
class TerrainShadowMap {
public:
    void update(Mesh& mesh, const Vector3& renderedExtent, size_t terrainRevision);
    static GridF lightMatrix(const Vector3& light, const Vector3& extent);
    static constexpr int resolution = 1024;
    size_t renderCount = 0;
    std::shared_ptr<GlobalsGL::Resources> resources;
    GLuint texture = 0;
private:
    GLuint framebuffer = 0;
    std::vector<unsigned char> key;
    std::array<GridF, 6> matrices;
};
#endif
