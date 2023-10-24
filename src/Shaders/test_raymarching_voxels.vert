#version 430

layout (location = 0) in vec2 aPos;

out vec2 v_texCoord;
out vec4 near_4;
out vec4 far_4;
out vec3 rayDirection;

struct rt_scene {
//	vec4 quat_camera_rotation;
    vec4 camera_pos;
    mat4 proj_matrix;
    mat4 mv_matrix;
    vec4 bg_color;

    int canvas_width;
    int canvas_height;

    int reflect_depth;
};

uniform rt_scene scene;
uniform sampler3D dataFieldTex;

void main()
{
    vec4 reverseVec = vec4(aPos, 0.0, 1.0);
    mat4 inverseProjection = inverse(scene.proj_matrix);
    mat4 inverseModelView = inverse(scene.mv_matrix);
    reverseVec = inverseProjection * reverseVec;

    reverseVec.w = 0.0;
    reverseVec = inverseModelView * reverseVec;

    rayDirection = /*vec4(scene.mv_matrix * scene.proj_matrix * vec4(aPos, 0.0, 1.0)).xyz; //*/vec3(reverseVec);

//    near_4 = vec4(rayDirection, .5);
//    far_4 = vec4(rayDirection, 1.0);

    mat4 inverse_matrix = inverse(scene.proj_matrix * scene.mv_matrix);
    near_4 = inverse_matrix * (vec4(aPos, -1.0, 1.0));
    far_4 = inverse_matrix * (vec4(aPos, +1.0, 1.0));
//    near_4 = scene.proj_matrix * scene.mv_matrix * (vec4(aPos, -1.0, 1.0));
//    far_4 = scene.proj_matrix * scene.mv_matrix * (vec4(aPos, +1.0, 1.0));
    v_texCoord = aPos;
    gl_Position = vec4(aPos, 0.0, 1.0);
}
