#version 430

out vec4 fragColor;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec4 color = vec4(1.0, 1.0, 1.0, 1.0);
uniform bool cullFace = true;

in vec3 vColor;

void main(void)
{
    if (!gl_FrontFacing && cullFace)
        discard;
    fragColor = (length(vColor) > 0 ? vec4(vColor, color.a) : color); // Use vertex color if possible, but keep alpha value
//    fragColor = color;
}
