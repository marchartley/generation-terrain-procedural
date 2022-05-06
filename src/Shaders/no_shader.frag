#version 120

//out vec4 fragColor;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform mat4 norm_matrix;

uniform vec4 color = vec4(1.0, 1.0, 1.0, 1.0);

void main(void)
{
//    if (!gl_FrontFacing)
//        discard;
    gl_FragColor = color;
}
