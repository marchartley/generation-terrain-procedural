#ifndef SHADERELEMENT_H
#define SHADERELEMENT_H

class ShaderElement;
class ColoredElement;
class LightSource;
class PositionalLight;
class Material;

class rt_material;
class rt_sphere;
class rt_plane;
class rt_box;
class rt_ring;
class rt_surface;
class rt_torus;
class rt_light_direct;
class rt_light_point;
class rt_scene;

#include "DataStructure/Vector3.h"
#include "Graphics/Shader.h"
#include <glm/glm.hpp>

class ShaderElement
{
public:
    ShaderElement();
    ShaderElement(std::string name, std::shared_ptr<Shader> shader = nullptr);
    ~ShaderElement();

    virtual void update();
    void affectShader(std::shared_ptr<Shader> shader);

    std::string name;

    std::shared_ptr<Shader> shader;
};

class ShaderUBO
{
public:
    ShaderUBO();
    ShaderUBO(std::string name, int binding, size_t size);
    ~ShaderUBO();

    void initBuffer(size_t size, void* data);
    void updateBuffer(size_t size, void* data);
    void affectShader(std::shared_ptr<Shader> shader);

    GLuint UBO;
    int binding;
    std::string name;

    std::shared_ptr<Shader> shader;

    size_t data_size;
};

class ColoredElement : public ShaderElement
{
public:
    ColoredElement();
    ColoredElement(float ambiant[4], float diffuse[4], float specular[4]);
    virtual void update();

    float ambiant[4];
    float diffuse[4];
    float specular[4];
};

class LightSource : public ColoredElement
{
public:
    LightSource();
    LightSource(float ambiant[4], float diffuse[4], float specular[4]);
    virtual void update();
};
class PositionalLight : public LightSource
{
public:
    PositionalLight();
    PositionalLight(float ambiant[4], float diffuse[4], float specular[4], const Vector3& pos);
    virtual void update();

    Vector3 position;
};

class Material : public ColoredElement
{
public:
    Material();
    Material(float ambiant[4], float diffuse[4], float specular[4], float shininess);
    virtual void update();

    float shininess;
};


using glm::vec4, glm::mat4, glm::vec3;

class rt_material : public ShaderElement {
public:
    rt_material();
    rt_material(std::string name, std::shared_ptr<Shader> shader = nullptr);
    vec4 color;
    vec4 absorb;

    float diffuse;
    float reflection;
    float refraction;
    int specular;
    float kd;
    float ks;

    virtual void update();
};

class rt_sphere : public ShaderElement {
public:
    rt_sphere();
    rt_sphere(std::string name, std::shared_ptr<Shader> shader = nullptr);
    rt_material mat;
    vec4 obj;
    vec4 quat_rotation; // rotate normal
    int textureNum;
    bool hollow;

    virtual void update();
};

class rt_plane : public ShaderElement {
public:
    rt_plane();
    rt_plane(std::string name, std::shared_ptr<Shader> shader = nullptr);

    rt_material mat;
    vec4 pos;
    vec4 normal;

    virtual void update();
};

class rt_box : public ShaderElement {
public:
    rt_box();
    rt_box(std::string name, std::shared_ptr<Shader> shader = nullptr);

    rt_material mat;
    vec4 quat_rotation;
    vec4 pos;
    vec4 form;
    int textureNum;

    virtual void update();
};

class rt_ring : public ShaderElement {
public:
    rt_ring();
    rt_ring(std::string name, std::shared_ptr<Shader> shader = nullptr);

    rt_material mat;
    vec4 quat_rotation;
    vec4 pos;
    int textureNum;
    float r1; // square of min radius
    float r2; // square of max radius

    virtual void update();
};

class rt_surface : public ShaderElement {
public:
    rt_surface();
    rt_surface(std::string name, std::shared_ptr<Shader> shader = nullptr);

    rt_material mat;
    vec4 quat_rotation;
    vec4 v_min;
    vec4 v_max;
    vec4 pos;
    float a; // x2
    float b; // y2
    float c; // z2
    float d; // z
    float e; // y
    float f; // const

    virtual void update();
};

class rt_torus : public ShaderElement {
public:
    rt_torus();
    rt_torus(std::string name, std::shared_ptr<Shader> shader = nullptr);

    rt_material mat;
    vec4 quat_rotation;
    vec4 pos;
    vec4 form; // x - radius, y - ring thickness

    virtual void update();
};

class rt_light_direct : public ShaderElement {
public:
    rt_light_direct();
    rt_light_direct(std::string name, std::shared_ptr<Shader> shader = nullptr);

    vec4 direction;
    vec4 color;

    float intensity;

    virtual void update();
};

class rt_light_point : public ShaderElement {
public:
    rt_light_point();
    rt_light_point(std::string name, std::shared_ptr<Shader> shader = nullptr);

    vec4 pos; //pos + radius
    vec4 color;
    float intensity;

    float linear_k;
    float quadratic_k;

    virtual void update();
};

class rt_scene : public ShaderElement {
public:
    rt_scene();
    rt_scene(std::string name, std::shared_ptr<Shader> shader = nullptr);
    vec4 camera_pos;
    mat4 proj_matrix;
    mat4 mv_matrix;
    vec4 bg_color;

    int canvas_width;
    int canvas_height;

    int reflect_depth;

    virtual void update();
};

#endif // SHADERELEMENT_H
