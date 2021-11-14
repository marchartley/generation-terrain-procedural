#version 430
layout (location=0) in vec3 position;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;

out vec4 varyingColor;

void main(void)
{
    float waterIndex = 20.0;
//    vec3 position = vec3(position.x, position.y, max(position.z, waterIndex-0.01));
    gl_Position = proj_matrix * mv_matrix * vec4(position, 1.0);
    if (position.z < waterIndex) {
        varyingColor = vec4(0.0, 0.0, 1.0, 1.0);
//        gl_Position.z = waterIndex;
    } else {
        varyingColor = vec4(0.0, 1.0, 0.0, 1.0);
    }
}
