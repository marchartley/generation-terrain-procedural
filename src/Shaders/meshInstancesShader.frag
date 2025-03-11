#version 430

out vec4 fragColor;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec4 color = vec4(1.0, 1.0, 1.0, 1.0);
uniform bool cullFace = true;


struct PositionalLight {
    vec4 ambiant;
    vec4 diffuse;
    vec4 specular;
    vec3 position;
};

in vec3 grealNormal;
in vec3 ginitialVertPos;
in vec4 gcolor;
in float gdensity;
in float gambiantOcclusion;

in vec4 varyingColor;
in vec3 varyingHalfH;

in vec3 varyingNormal;
in vec3 varyingLightDir;
in vec3 varyingVertPos;
in vec3 vertEyeSpacePos;
in vec3 initialVertPos;

uniform vec4 globalAmbiant;
const int nbLights = 6;
uniform PositionalLight lights[nbLights];

void main(void)
{
    if (!gl_FrontFacing && cullFace)
        discard;


    vec3 N = normalize(varyingNormal);
    for (int iLight = 0; iLight < nbLights; iLight++) {

        vec3 light_position = vec4(mv_matrix * vec4(lights[iLight].position, 1.0)).xyz;
        vec3 varyingLightDir = light_position - varyingVertPos;
        vec3 varyingHalfH = (varyingLightDir - varyingVertPos).xyz;
        vec3 L = normalize(varyingLightDir);
        vec3 R = reflect(-L, N);
        vec3 H = normalize(varyingHalfH);
        float cosTheta = dot(L, N);
        float cosPhi = dot(H, N);

        float lumin = 1.0;

        vec4 ambiant = ((globalAmbiant * color) + (lights[iLight].ambiant * color));
        vec4 diffuse = lights[iLight].diffuse * color * max(cosTheta, 0.0);
        vec4 specular = vec4(0.0, 0.0, 0.0, 1.0); //lights[iLight].specular * material_specular * pow(max(cosPhi, 0.0), material_shininness * 32.0);
        vec4 material_color = vec4((ambiant + diffuse + specular).xyz*2, 1.0);
        vec3 light_pos = vec4(mv_matrix * vec4(lights[iLight].position, 1.0)).xyz;
        float dist_source = length(lights[iLight].position - ginitialVertPos) / 50.0;

        fragColor += vec4((ambiant + diffuse + specular).xyz, 1.0);
    }
    fragColor /= float(nbLights);

//    fragColor = color;
}
