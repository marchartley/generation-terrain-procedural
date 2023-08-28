#version 430
float random (in vec2 st) {
    return fract(sin(dot(st.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}

// Based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float hash(float p) { p = fract(p * 0.011); p *= p + 7.5; p *= p + p; return fract(p); }
float hash(vec2 p) {vec3 p3 = fract(vec3(p.x, p.y, p.x) * 0.13); p3 += dot(p3, p3.yzx + 3.333); return fract((p3.x + p3.y) * p3.z); }

float noise(float x) {
    float i = floor(x);
    float f = fract(x);
    float u = f * f * (3.0 - 2.0 * f);
    return mix(hash(i), hash(i + 1.0), u);
}
float noise (in vec2 st) {
    vec2 i = floor(st);
    vec2 f = fract(st);

    // Four corners in 2D of a tile
    float a = random(i);
    float b = random(i + vec2(1.0, 0.0));
    float c = random(i + vec2(0.0, 1.0));
    float d = random(i + vec2(1.0, 1.0));

    vec2 u = f * f * (3.0 - 2.0 * f);

    return mix(a, b, u.x) +
            (c - a)* u.y * (1.0 - u.x) +
            (d - b) * u.x * u.y;
}
float noise3(vec3 x) {
    const vec3 step = vec3(110, 241, 171);

    vec3 i = floor(x);
    vec3 f = fract(x);

    // For performance, compute the base input to a 1D hash from the integer part of the argument and the
    // incremental change to the 1D based on the 3D -> 1D wrapping
    float n = dot(i, step);

    vec3 u = f * f * (3.0 - 2.0 * f);
    return mix(mix(mix( hash(n + dot(step, vec3(0, 0, 0))), hash(n + dot(step, vec3(1, 0, 0))), u.x),
                   mix( hash(n + dot(step, vec3(0, 1, 0))), hash(n + dot(step, vec3(1, 1, 0))), u.x), u.y),
               mix(mix( hash(n + dot(step, vec3(0, 0, 1))), hash(n + dot(step, vec3(1, 0, 1))), u.x),
                   mix( hash(n + dot(step, vec3(0, 1, 1))), hash(n + dot(step, vec3(1, 1, 1))), u.x), u.y), u.z);
}

float fbm (in vec2 st) {
    // Initial values
    float value = 0.0;
    float amplitud = .5;
    float frequency = 0.;
    //
    // Loop of octaves
    for (int i = 0; i < 6; i++) {
        value += amplitud * noise(st);
        st *= 2.;
        amplitud *= .5;
    }
    return value;
}
float fbm3 (in vec3 st) {
    // Initial values
    float value = 0.0;
    float amplitud = .5;
    float frequency = 0.;
    //
    // Loop of octaves
    for (int i = 0; i < 6; i++) {
        value += amplitud * noise3(st);
        st *= 2.;
        amplitud *= .5;
    }
    return value;
}

vec2 fbmToVec2(in vec2 st) {
    return vec2(fbm(st), fbm(st + 100.0));
}
vec2 fbm3ToVec2(in vec3 st) {
    return vec2(fbm3(st), fbm3(st + 100.0));
}
vec3 fbmToVec3(in vec2 st) {
    return vec3(fbm(st), fbm(st + 18.0), fbm(st - 10.0));
}
vec3 fbm3ToVec3(in vec3 st) {
    return vec3(fbm3(st), fbm3(st + 18.0), fbm3(st - 10.0));
}

uniform vec3 min_vertice_positions;
uniform vec3 max_vertice_positions;

uniform sampler3D dataFieldTex;
uniform sampler3D dataChangesFieldTex;

uniform float fogNear;
uniform float fogFar;

uniform bool clipPlaneActive;
uniform vec3 clipPlanePosition;
uniform vec3 clipPlaneDirection;

uniform float waterRelativeHeight = 0.0;

uniform vec3 subterrainOffset = vec3(0, 0, 0);
uniform float subterrainScale = 1.0;

out vec4 fragColor;

in vec3 grealNormal;
in vec3 ginitialVertPos;
in vec4 gcolor;
in float gdensity;

struct PositionalLight {
    vec4 ambiant;
    vec4 diffuse;
    vec4 specular;
    vec3 position;
};

struct Material {
    vec4 ambiant;
    vec4 diffuse;
    vec4 specular;
    float shininness;
};

uniform vec4 globalAmbiant;
uniform PositionalLight light;
uniform bool isSpotlight;
uniform Material ground_material;
uniform Material grass_material;

uniform bool displayingIgnoredVoxels = false;
uniform bool wireframeMode = false;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;


uniform sampler2D biomeFieldTex;
uniform sampler2D heightmapFieldTex;
uniform sampler2D biomesAndDensitiesTex;
uniform sampler2D allBiomesColorTextures;
uniform int maxBiomesColorTextures;
uniform sampler2D allBiomesNormalTextures;
uniform int maxBiomesNormalTextures;

//vec4 sandColor = vec4(.761, .698, .502, 1.0);
//vec4 mudColor = vec4(.430, .320, .250, 1.0);
//vec4 rockColor = vec4(.755, .754, .748, 1.0);

//float densities[5] = float[5](0.3, 0.5, 0.8, 0.9, 1.2);
//vec4 colors[5] = vec4[5](vec4(.761, .698, .502, 1.0), vec4(.661, .598, .402, 1.0), vec4(.630, .320, .250, 1.0), vec4(.430, .320, .250, 1.0), vec4(.755, .754, .748, 1.0));

//float densities[4] = float[4](0.3, 0.5, 0.6, 0.7);
//vec4 colors[4] = vec4[4](vec4(.761, .698, .502, 1.0), vec4(.661, .598, .402, 1.0), vec4(.630, .320, .250, 1.0), vec4(.755, .754, .748, 1.0));
//int texIndices[4] = int[4](6, 5, 4, 2);

float densities[5] = float[5](0.1, 0.2, 0.6, 1.0, 2.0);
vec4 colors[5] = vec4[5](vec4(.761, .698, .502, 1.0), vec4(.761, .698, .502, 1.0), vec4(.661, .598, .402, 1.0), vec4(.630, .320, .250, 1.0), vec4(.755, .754, .748, 1.0));
int texIndices[5] = int[5](2, 6, 1, 5, 4);

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

float sigmoid(float x, float _null) {
    // Approximaton
    return 0.5 + x * (1.0 - abs(x) * 0.5);
}

vec3 getBiomeColor(vec2 pos) {
    vec2 texSize = textureSize(biomeFieldTex, 0);
    float val = max(texture(biomeFieldTex, pos/texSize).r, max(texture(biomeFieldTex, pos/texSize).g, texture(biomeFieldTex, pos/texSize).b));
//    return hsv2rgb(vec3(val, 1.0, 1.0));
    return vec3(val, 1.0, 1.0);
}
float getDensity(vec3 pos, float resolution) {
    if (gdensity <= 0) {
        vec3 texSize = vec3(textureSize(dataFieldTex, 0));
        vec3 offsets = vec3(0.0); // 0.5 / texSize;
        float density = texture(dataFieldTex, pos).a;
        float surrounding = 2.f;
        for (float x = 0; x < surrounding + 1.0; x += 1.0) {
            for (float y = 0; y < surrounding + 1.0; y += 1.0) {
                for (float z = 0; z < surrounding + 1.0; z += 1.0) {
                    float divisor = 0.1; //float(int(surrounding*.5f));
                    vec3 newPos = pos + vec3(x - surrounding*.5f, y - surrounding*.5f, z - surrounding*.5f) / texSize; // * (resolution / divisor);
    //                vec3 newPos = pos + vec3(resolution * (x/surrounding - surrounding * .5f), resolution * (y/surrounding - surrounding * .5f), resolution * (z/surrounding - surrounding * .5f));
                    float val = texture(dataFieldTex, newPos - offsets).a;
                    density = max(density, val);
                }
            }
        }
        return (density - 0.5);
    } else {
        return gdensity;
    }
}
int getDensityIndex(vec3 pos, vec3 terrainSize, float depth) {
//    return 1;
    float density = getDensity((pos /*- vec3(1)*/) / (terrainSize /*+ vec3(1)*/), 0.0); // 0.1 / terrainSize.x); //1.0 / terrainSize.x);
    float power = 1.0;
    int texIndex = 0;
    if (density <= densities[0]) {
        texIndex = texIndices[0];
    } else if (density >= densities[densities.length() - 1]) {
        texIndex = texIndices[densities.length() - 1];
    } else {
        for (int i = 0; i < densities.length() - 1; i++)
            texIndex = (densities[i] <= density && density < densities[i + 1] ? (density < (densities[i] + densities[i + 1])/2.0 ? texIndices[i] : texIndices[i+1]) : texIndex);
    }
    return texIndex;
}
/*
vec4 getDensityColor(vec3 pos, vec3 terrainSize) {
    float density = getDensity(pos/ terrainSize, 1.0 / terrainSize.x);
    float power = 1.0;
    vec4 color = vec4(1.0, 1.0, 1.0, 1.0);
    if (density <= densities[0]) {
        color *= colors[0];
    } else if (density >= densities[densities.length() - 1]) {
        color *= colors[colors.length() - 1];
    } else {
        for (int i = 0; i < densities.length() - 1; i++)
            color *= (densities[i] <= density && density < densities[i + 1] ? mix(colors[i], colors[i + 1], (density - densities[i])/(densities[i + 1] - densities[i])) : vec4(1.0, 1.0, 1.0, 1.0));
    }
//    if (density <= sandDensity) {
//        color *= sandColor;
//    } else if (density > sandDensity && density <= mudDensity) {
//        color *= mix(sandColor, mudColor, pow((density - sandDensity)/(mudDensity - sandDensity), power));
//    } else if (density > mudDensity && density <= rockDensity) {
//        color *= mix(mudColor, rockColor, pow((density - mudDensity)/(rockDensity - mudDensity), power));
//    } else  {
//        color *= rockColor;
//    }
    return color;
}
*/

vec3 getTriPlanarBlend(vec3 _wNorm){
    // in wNorm is the world-space normal of the fragment
    vec3 blending = abs( _wNorm );
    blending = normalize(max(blending, 0.00001)); // Force weights to sum to 1.0
    float b = (blending.x + blending.y + blending.z);
    blending /= vec3(b, b, b);
    return blending;
}

float wyvill(float x) {
    return pow((1 - pow(x, 2)), 3);
}
void main(void)
{
    if (clipPlaneActive) {
        if (dot((ginitialVertPos.xyz - clipPlanePosition), clipPlaneDirection) > 0) {
            discard;
        }
    }

    if (wireframeMode && !gl_FrontFacing)
        discard;

    Material material = ground_material;
    vec3 position = ginitialVertPos.xyz;
//    varyingColor = vec4(1.0, 1.0, 1.0, 1.0);
    vec3 light_position = vec4(mv_matrix * vec4(light.position, 1.0)).xyz;
    vec3 varyingVertPos = vec4(mv_matrix * vec4(position, 1.0)).xyz;
    vec3 varyingLightDir = light_position - varyingVertPos;
    vec3 varyingNormal = vec4(transpose(inverse(mv_matrix)) * vec4(grealNormal, 1.0)).xyz;
    vec3 varyingHalfH = (varyingLightDir - varyingVertPos).xyz;
    vec3 vertEyeSpacePos = vec4(mv_matrix * vec4(position, 1.0)).xyz;
    vec3 N = normalize(varyingNormal /*+ fbm3ToVec3(ginitialVertPos)*0.5*/);
    vec3 L = normalize(varyingLightDir);
    vec3 V = normalize(-varyingVertPos);
    vec3 R = reflect(-L, N);
    vec3 H = normalize(varyingHalfH);

    vec3 realFragmentPosition = (ginitialVertPos.xyz / vec3(subterrainScale, subterrainScale, 1.0)) + subterrainOffset;
    vec3 dataTexSize = textureSize(dataFieldTex, 0);

    vec2 texSize = textureSize(heightmapFieldTex, 0);
    vec2 fbmWrap = clamp(ginitialVertPos.xy + 2.0 * (fbm3ToVec2(realFragmentPosition) - vec2(.5, .5)) * 5.0 * subterrainScale, vec2(0.51, 0.51), texSize - vec2(0.51, 0.51)) - ginitialVertPos.xy;
    float biomeColorValue = texture(heightmapFieldTex, (realFragmentPosition.xy + fbmWrap)/texSize).r;
    float realBiomeColorValue = texture(heightmapFieldTex, (realFragmentPosition.xy)/texSize).r;
    float biomeNormalValue = texture(heightmapFieldTex, (realFragmentPosition.xy + fbmWrap)/texSize).r;

    float depth = (waterRelativeHeight * textureSize(dataFieldTex, 0).z - realFragmentPosition.z); // Depth from surface

//    fragColor = vec4(maxBiomesColorTextures / 7.0, 0.0, 0.0, 1.0);
//    return;
    realBiomeColorValue = float(getDensityIndex(realFragmentPosition, dataTexSize, depth)) / maxBiomesColorTextures;
    biomeColorValue = realBiomeColorValue; //biomeColorValue = float(getDensityIndex(realFragmentPosition + fbm3ToVec3(realFragmentPosition), dataTexSize, depth)) / maxBiomesColorTextures;
    biomeNormalValue = biomeColorValue;
//    fragColor = vec4((ginitialVertPos.xy/texSize).x, 0.0, (ginitialVertPos.xy/texSize).y, 1.0);
    float scale = 10.0;
    vec3 blending = getTriPlanarBlend(grealNormal);
    if (biomeNormalValue < 1.0) {

        vec3 Nx = vec3(1, 0, 0);
        vec3 Ny = vec3(0, 1, 0);
        vec3 Nz = vec3(0, 0, 1);

        vec3 Tx = vec3(0, 0, 1);
        vec3 Ty = vec3(1, 0, 0);
        vec3 Tz = vec3(0, 1, 0);

        vec3 Bx = cross(Nx, Tx); // vec3(0, 1, 0)
        vec3 By = cross(Ny, Ty); // vec3(0, 0, 1)
        vec3 Bz = cross(Nz, Tz); // vec3(1, 0, 0)

        mat3 TBNx = transpose(mat3(Tx, Bx, Nx));
        mat3 TBNy = transpose(mat3(Ty, By, Ny));
        mat3 TBNz = transpose(mat3(Tz, Bz, Nz));

        vec2 normalTextureOffset = vec2(biomeNormalValue, 0);
        vec3 xaxis = texture(allBiomesNormalTextures, normalTextureOffset + (fract(realFragmentPosition.yz / scale) * vec2(1.0/maxBiomesNormalTextures, 1.0)) * 0.99).rgb;
        vec3 yaxis = texture(allBiomesNormalTextures, normalTextureOffset + (fract(realFragmentPosition.xz / scale) * vec2(1.0/maxBiomesNormalTextures, 1.0)) * 0.99).rgb;
        vec3 zaxis = texture(allBiomesNormalTextures, normalTextureOffset + (fract(realFragmentPosition.xy / scale) * vec2(1.0/maxBiomesNormalTextures, 1.0)) * 0.99).rgb;
        vec3 normalMap = (TBNx * xaxis * blending.x + TBNy * yaxis * blending.y + TBNz * zaxis * blending.z);
        N = vec4(transpose(inverse(mv_matrix)) * vec4(normalMap, 1.0)).xyz + N; // = normalize(N + normalMap * sign(N));
    }
    float cosTheta = dot(L, N);
    float cosPhi = dot(H, N);

    float upValue = 1.0; // clamp(dot(normalize(grealNormal), normalize(vec3(0.0, 0.0, 1.0) /*+ fbm3ToVec3(ginitialVertPos) * 0.5*/)), 0.0, 1.0);
    vec4 material_ambiant = (ground_material.ambiant * (1 - upValue) + grass_material.ambiant * upValue) * .8;
    vec4 material_diffuse = (ground_material.diffuse * (1 - upValue) + grass_material.diffuse * upValue) * .8;
    vec4 material_specular = (ground_material.specular * (1 - upValue) + grass_material.specular * upValue) * 1.;
    float material_shininness = ground_material.shininness * (1 - upValue) + grass_material.shininness * upValue;

    vec4 ambiant = ((globalAmbiant * material_ambiant) + (light.ambiant * material_ambiant));
    vec4 diffuse = light.diffuse * material_diffuse * max(cosTheta, 0.0);
    vec4 specular = light.specular * material_specular * pow(max(cosPhi, 0.0), material_shininness * 32.0);

    vec4 material_color = vec4((ambiant + diffuse + specular).xyz*2, 1.0);

    vec3 cam_pos = vec4(mv_matrix * vec4(0.0, 0.0, 0.0, 1.0)).xyz;
    vec3 light_pos = vec4(mv_matrix * vec4(light.position, 1.0)).xyz;
    float dist_source = length(light.position - ginitialVertPos) / 50.0;
    float dist_cam = length(cam_pos - vertEyeSpacePos) / 50.0;
//    float albedo = 3.14;
    float albedo = 1.0;

    vec4 attenuation_coef = vec4(vec3(0.000245, 0.0027, 0.0020) * depth, 1.0);
    vec4 absorbtion_coef = vec4(vec3(0.210, 0.075, 0.095), 1.0);
//    vec4 attenuation_coef = vec4(vec3(0.245, 0.027, 0.020) * depth, 1.0);
//    vec4 absorbtion_coef = vec4(vec3(0.210, 0.0075, 0.0005) * depth, 1.0); // Source : http://web.pdx.edu/~sytsmam/limno/Limno09.7.Light.pdf
    vec4 water_light_attenuation = exp(-attenuation_coef * dist_cam)*(albedo/3.1415)*cosTheta*(1-absorbtion_coef)*exp(-attenuation_coef*dist_source);
    vec4 fogColor = vec4((vec4(water_light_attenuation.xyz, 1.0) * material_color).xyz, 1.0);

    float fogStart = fogNear;
    float fogEnd = fogFar;
    float dist = length(vertEyeSpacePos) / 100000;
    float fogFactor = 1.0; // depth < 0.0 ? clamp(((1 - (depth / 10.0)) * (fogEnd - dist) / (fogEnd - fogStart)), 0.0, 1.0) : 0.0;

    float lumin = 1.0;
    if (isSpotlight && false) {
        float cone_angle = radians(20.0);
        float cone_cos_angle = cos(cone_angle);
        float epsilon = radians(10.0);

        lumin = clamp(1 - (acos(dot(normalize(varyingLightDir), vec3(0.0, 0.0, 1.0))) - (cone_angle + epsilon))/(epsilon), 0.0, 1.0);
    }
    if (!gl_FrontFacing)
        lumin *= 0.6;
    fragColor = vec4(mix(fogColor, material_color, fogFactor).xyz * lumin, 1.0) * gcolor;

//    fragColor = material_color * gcolor;
//    if (depth > 0)
//        fragColor *= vec4(0.8, 1.0, 0.9, 1.0);
//    fragColor = vec4(texture(dataFieldTex, (realFragmentPosition + vec3(-1, -1, -1)) / dataTexSize, 0).a,
//                     texture(dataFieldTex, (realFragmentPosition + vec3(0, 0, 0)) / dataTexSize, 0).a,
//                     texture(dataFieldTex, (realFragmentPosition + vec3(1, 1, 1)) / dataTexSize, 0).a,
//                     1.0);
//    return;
//    fragColor *= getDensityColor(realFragmentPosition, dataTexSize);
//    fragColor *= (depth < 0.0 ? vec4(0.5, 0.8, 0.9, 1.0) : vec4(1.0, 1.0, 1.0, 1.0));

    if (displayingIgnoredVoxels) {
        fragColor = vec4(0, 0, 0, 0.1);
        return;
    }

    vec4 colorErod = vec4(158, 42, 43, 0) / 255.0;
    vec4 colorNull = vec4(242, 184, 128, 0) / 255.0;
    vec4 colorDepo = vec4(183, 195, 243, 0) / 255.0;

    vec3 changeSize = vec3(textureSize(dataChangesFieldTex, 0));
    vec3 evalPos = (changeSize.z > 1 ? ginitialVertPos / changeSize : vec3(ginitialVertPos.xy, 0.5) / changeSize);
    float changeVal = texture(dataChangesFieldTex, evalPos).a - 2.0;

    vec4 col = mix(colorErod, colorDepo, (changeVal + 1.0) / 2.0);
    col = mix(col, colorNull, 1.0 - abs(changeVal));
    fragColor = vec4((col * (ambiant + diffuse + specular)).xyz * 3.0, 1.0);
    fragColor = vec4(fragColor.xyz * (realFragmentPosition.z > waterRelativeHeight * dataTexSize.z ? vec3(1.0) : vec3(0.8, 1.1, 1.5)), 1.0);
    return;

    if (biomeColorValue < 1.0) {
        vec2 colorTextureOffset     = vec2(biomeColorValue, 0);
        vec2 realColorTextureOffset = vec2(realBiomeColorValue, 0);
        // biome val == real biome val => no need to have a mix (alpha = 0)
        float alpha = 1 - clamp(wyvill(length(fbmWrap)/5.0), 0.0, 1.0) * (biomeColorValue == realBiomeColorValue ? 0.0 : 1.0);
        vec3 xaxis = mix(texture(allBiomesColorTextures, colorTextureOffset     + (fract(realFragmentPosition.yz / scale) * vec2(1.0/maxBiomesColorTextures, 1.0)) * 0.99),
                         texture(allBiomesColorTextures, realColorTextureOffset + (fract(realFragmentPosition.yz / scale) * vec2(1.0/maxBiomesColorTextures, 1.0)) * 0.99),
                         alpha).rgb;
        vec3 yaxis = mix(texture(allBiomesColorTextures, colorTextureOffset     + (fract(realFragmentPosition.xz / scale) * vec2(1.0/maxBiomesColorTextures, 1.0)) * 0.99),
                         texture(allBiomesColorTextures, realColorTextureOffset + (fract(realFragmentPosition.xz / scale) * vec2(1.0/maxBiomesColorTextures, 1.0)) * 0.99),
                         alpha).rgb;
        vec3 zaxis = mix(texture(allBiomesColorTextures, colorTextureOffset     + (fract(realFragmentPosition.xy / scale) * vec2(1.0/maxBiomesColorTextures, 1.0)) * 0.99),
                         texture(allBiomesColorTextures, realColorTextureOffset + (fract(realFragmentPosition.xy / scale) * vec2(1.0/maxBiomesColorTextures, 1.0)) * 0.99),
                         alpha).rgb;
        fragColor = vec4((vec4(blending.x * xaxis + blending.y * yaxis + blending.z * zaxis, 1.0) * (ambiant + diffuse + specular)).xyz * 3.0, 1.0);

//        fragColor = vec4(fragColor.xyz * waterRelativeHeight, 1.0);
        fragColor = vec4(fragColor.xyz * (realFragmentPosition.z > waterRelativeHeight * dataTexSize.z ? vec3(1.0) : vec3(0.8, 1.1, 1.5)), 1.0);
    }
}
