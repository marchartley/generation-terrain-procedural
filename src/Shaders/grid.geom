#version 430
//New G80 extensions
//#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

layout ( points ) in;
layout ( triangle_strip, max_vertices = 40) out;

//height data field texture
uniform sampler2D heightmapFieldTex;
uniform float maxHeight;
uniform vec2 viewportSize = vec2(0);
uniform bool adaptiveTessellation = true;
//uniform sampler2D allBiomesColorTextures;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec3 min_vertice_positions;
uniform vec3 max_vertice_positions;
uniform float ambiantOcclusionFactor = 0.0;
uniform float heightFactor = 0.1;
//Vertices position for fragment shader
in vec3 initialVertPos[];
in vec3 realNormal[];

out vec3 ginitialVertPos;
out vec3 grealNormal;
out vec4 gcolor;
out float gambiantOcclusion;
out float gdensity;


float linearInterpolate(sampler2D tex, vec2 uv, vec2 texSize, float maxHeight) {
    vec2 texelSize = 1.0 / texSize;
    uv = uv * texSize; // Transform to texel space

    vec2 iuv = floor(uv);
    vec2 fuv = fract(uv);

    vec2 offset = vec2(1.0, 1.0);

    float s00 = texture(tex, (iuv + vec2(0.0, 0.0)) * texelSize).a * maxHeight;
    float s10 = texture(tex, (iuv + vec2(1.0, 0.0)) * texelSize).a * maxHeight;
    float s01 = texture(tex, (iuv + vec2(0.0, 1.0)) * texelSize).a * maxHeight;
    float s11 = texture(tex, (iuv + vec2(1.0, 1.0)) * texelSize).a * maxHeight;

    float s0 = mix(s00, s10, fuv.x);
    float s1 = mix(s01, s11, fuv.x);

    return mix(s0, s1, fuv.y);
}

float cubicInterpolate(float p0, float p1, float p2, float p3, float t) {
    return p1 + 0.5 * t*(p2 - p0 + t*(2.0*p0 - 5.0*p1 + 4.0*p2 - p3 + t*(3.0*(p1 - p2) + p3 - p0)));
}

float bicubicInterpolate(sampler2D tex, vec2 uv, vec2 texSize, float maxHeight) {
    float cubicScale = 1.0;
    vec2 texelSize = cubicScale / texSize;
    uv = uv * texSize - 0.5 * cubicScale; // Transform to texel space

    vec2 iuv = floor(uv);
    vec2 fuv = fract(uv);

    float samples[4][4];

    float radius = 1.0;
    for(int y = 0; y < 4; y++) {
        for(int x = 0; x < 4; x++) {
            vec2 samplePos = (iuv + vec2(x - 1.0, y - 1.0) * radius) * texelSize;
            samplePos.x = clamp(samplePos.x, 0.0, 1.0);
            samplePos.y = clamp(samplePos.y, 0.0, 1.0);
            // samplePos *= texelSize;
            samples[y][x] = texture(tex, samplePos / cubicScale).a * maxHeight;
        }
    }

    float col0 = cubicInterpolate(samples[0][0], samples[0][1], samples[0][2], samples[0][3], fuv.x);
    float col1 = cubicInterpolate(samples[1][0], samples[1][1], samples[1][2], samples[1][3], fuv.x);
    float col2 = cubicInterpolate(samples[2][0], samples[2][1], samples[2][2], samples[2][3], fuv.x);
    float col3 = cubicInterpolate(samples[3][0], samples[3][1], samples[3][2], samples[3][3], fuv.x);

    return cubicInterpolate(col0, col1, col2, col3, fuv.y);
}

float getHeight(vec2 pos) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    // return linearInterpolate(heightmapFieldTex, pos / texSize, texSize, maxHeight);
    return texture(heightmapFieldTex, (pos + vec2(.5)) / texSize).a * maxHeight;
}
float getHeightLinear(vec2 pos) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    // return linearInterpolate(heightmapFieldTex, pos / texSize, texSize, maxHeight);
    return texture(heightmapFieldTex, (pos + vec2(.5)) / texSize).a * maxHeight;
}
/*float getHeight(vec2 pos) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    return (pos.x >= texSize.x || pos.y >= texSize.y) ? 0 : texture(heightmapFieldTex, pos/texSize).a * maxHeight;
}*/
float getDisplacementLength(vec2 pos) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    return (pos.x >= texSize.x || pos.y >= texSize.y) ? 0 : texture(heightmapFieldTex, pos/texSize).b;
}

vec3 getNormal(vec2 pos) {
    float h00 = getHeight(pos + vec2(-1,-1));
    float h10 = getHeight(pos + vec2( 0,-1));
    float h20 = getHeight(pos + vec2( 1,-1));
    float h01 = getHeight(pos + vec2(-1, 0));
    float h21 = getHeight(pos + vec2( 1, 0));
    float h02 = getHeight(pos + vec2(-1, 1));
    float h12 = getHeight(pos + vec2( 0, 1));
    float h22 = getHeight(pos + vec2( 1, 1));
    return normalize(vec3(h00 + 2.0*h01 + h02 - h20 - 2.0*h21 - h22,
                          h00 + 2.0*h10 + h20 - h02 - 2.0*h12 - h22, 8.0));
}
float getAmbiantOcclusion(vec3 pos, vec3 normal) {
    if (ambiantOcclusionFactor == 0.0) return 1.0;
    int nX = 5;
    int nY = 5;
    float fnX = float(nX);
    float fnY = float(nY);
    float pi = 3.141592;

    float maxDistance = 10.0;

    float occlusion = 0.0;
    float total = 0.0;
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nY; j++) {
            float theta = (i / (fnX)) * 2.0 * pi;
            float phi = (j / (fnY - 1)) * pi - (pi * 0.5);
            vec3 ray = vec3(cos(theta) * cos(phi), sin(theta) * cos(phi), sin(phi));
            ray *= (dot(ray, normal) <= 0 ? -1.0 : 1.0);

            bool good = true;
            for (int _r = 1; _r <= int(maxDistance) && good; _r++) {
                float r = _r;
                good = (pos.z + ray.z * r < getHeightLinear((pos + ray * r).xy) ? false : true);
            }
            occlusion += (good ? 1.0 : 0.0);
            total += 1;
        }
    }
    return occlusion / total; // pow(clamp(1.0 - smoothstep(0.0, 1.0, occlusion / total) + 0.5, 0.0, 1.0), 2.0);
}

void emitTerrainVertex(vec2 xy) {
    vec3 position = vec3(xy, getHeight(xy) * heightFactor);
    ginitialVertPos = position;
    gdensity = .65;
    grealNormal = getNormal(xy);
    gcolor = vec4(texture(heightmapFieldTex, xy / vec2(textureSize(heightmapFieldTex, 0))).rgb, 1.0);
    gambiantOcclusion = getAmbiantOcclusion(position, grealNormal);
    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);
    EmitVertex();
}

void main(void) {
    vec2 texSize = textureSize(heightmapFieldTex, 0);
    vec2 base = initialVertPos[0].xy;
    if (any(greaterThanEqual(base, texSize - vec2(1)))) return;
    if (any(lessThan(base, min_vertice_positions.xy)) || any(greaterThan(base, max_vertice_positions.xy))) return;
    int subdivisions = 4;
    if (adaptiveTessellation && all(greaterThan(viewportSize, vec2(0)))) {
        float h00 = getHeight(base), h10 = getHeight(base + vec2(1,0));
        float h01 = getHeight(base + vec2(0,1)), h11 = getHeight(base + vec2(1,1));
        float centerHeight = (h00 + h10 + h01 + h11) * .25;
        float triangleHeight = (h10 + h01) * .5;
        vec4 a = proj_matrix * mv_matrix * vec4(base + .5, centerHeight * heightFactor, 1.0);
        vec4 b = proj_matrix * mv_matrix * vec4(base + .5, triangleHeight * heightFactor, 1.0);
        if (a.w > 0.0 && b.w > 0.0) {
            float errorPixels = length((a.xy / a.w - b.xy / b.w) * viewportSize * .5);
            subdivisions = errorPixels <= .35 ? 1 : (errorPixels <= 1.4 ? 2 : 4);
        }
    }
    // The bilinear saddle error drops quadratically under subdivision. Smooth
    // cells emit one quad; sharp saddles keep the previous quarter-cell detail.
    float stepSize = 1.0 / float(subdivisions);
    for (int row = 0; row < subdivisions; ++row) {
        for (int col = 0; col <= subdivisions; ++col) {
            emitTerrainVertex(base + vec2(row, col) * stepSize);
            emitTerrainVertex(base + vec2(row + 1, col) * stepSize);
        }
        EndPrimitive();
    }
}
