#include "Graphics/ShaderElement.h"

ShaderElement::ShaderElement()
{

}

ShaderElement::ShaderElement(std::string name, std::shared_ptr<Shader> shader)
    : name(name)
{
    this->affectShader(shader);
}

ShaderElement::~ShaderElement()
{
}

void ShaderElement::update()
{
}

void ShaderElement::affectShader(std::shared_ptr<Shader> shader)
{
    this->shader = shader;
}

ShaderUBO::ShaderUBO()
{

}

ShaderUBO::ShaderUBO(std::string name, int binding, size_t size)
    : binding(binding), name(name), data_size(size)
{
}

ShaderUBO::~ShaderUBO()
{
    GlobalsGL::f()->glDeleteBuffers(1, &UBO);
}

void ShaderUBO::initBuffer(size_t size, void* data)
{
    this->shader->use();
    GlobalsGL::f()->glGenBuffers(1, &UBO);
    GlobalsGL::f()->glBindBuffer(GL_UNIFORM_BUFFER, UBO);
    GlobalsGL::f()->glBufferData(GL_UNIFORM_BUFFER, size, data, GL_STATIC_DRAW);
    GLuint blockIndex = GlobalsGL::f()->glGetUniformBlockIndex(shader->programID, name.c_str());
    if (blockIndex == 0xffffffff)
    {
        std::cout << "Invalid ubo block name '%s'" << name << std::endl;
//        exit(1);
    }
    GlobalsGL::f()->glUniformBlockBinding(shader->programID, blockIndex, binding);
    GlobalsGL::f()->glBindBufferBase(GL_UNIFORM_BUFFER, binding, UBO);
    GlobalsGL::f()->glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void ShaderUBO::updateBuffer(size_t size, void* data)
{
    this->shader->use();
    GlobalsGL::f()->glBindBuffer(GL_UNIFORM_BUFFER, UBO);
    GlobalsGL::f()->glBufferSubData(GL_UNIFORM_BUFFER, 0, size, data);
    GlobalsGL::f()->glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void ShaderUBO::affectShader(std::shared_ptr<Shader> shader)
{
    this->shader = shader;
    this->initBuffer(data_size, NULL);
}

ColoredElement::ColoredElement()
{
    Vector4 defaultAmbiant = {.8, .8, .8, 1.0};
    Vector4 defaultDiffuse = {.5, .5, .5, 1.0};
    Vector4 defaultSpecular = {.5, .5, .5, 1.0};
    ColoredElement(defaultAmbiant, defaultDiffuse, defaultSpecular);
}
ColoredElement::ColoredElement(const Vector4& ambiant, const Vector4& diffuse, const Vector4& specular)
{
    for (int i = 0; i < 4; i++)
    {
        this->ambiant[i] = ambiant[i];
        this->diffuse[i] = diffuse[i];
        this->specular[i] = specular[i];
    }
}

ColoredElement::~ColoredElement()
{
}

void ColoredElement::update()
{
    this->shader->setVector((name + ".ambiant").c_str(), ambiant);
    this->shader->setVector((name + ".diffuse").c_str(), diffuse);
    this->shader->setVector((name + ".specular").c_str(), specular);
}
LightSource::LightSource()
{

}
LightSource::LightSource(const Vector4& ambiant, const Vector4& diffuse, const Vector4& specular)
    : ColoredElement(ambiant, diffuse, specular)
{

}

void LightSource::update()
{
    ColoredElement::update();
}
PositionalLight::PositionalLight()
{

}
PositionalLight::PositionalLight(const Vector4& ambiant, const Vector4& diffuse, const Vector4& specular, const Vector3& pos)
    : LightSource(ambiant, diffuse, specular), position(pos)
{

}

void PositionalLight::update()
{
    LightSource::update();
    this->shader->setVector((name + ".position").c_str(), position);
}
Material::Material()
{

}
Material::Material(const Vector4& ambiant, const Vector4& diffuse, const Vector4& specular, float shininess)
    : ColoredElement(ambiant, diffuse, specular), shininess(shininess)
{

}

void Material::update()
{
    ColoredElement::update();
    this->shader->setFloat((name + ".shininness").c_str(), shininess);
}

rt_material::rt_material() : ShaderElement()
{}
rt_material::rt_material(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}

void rt_material::update()
{
    this->shader->setVector((name + ".color").c_str(), color);
    this->shader->setVector((name + ".absorb").c_str(), absorb);

    this->shader->setFloat((name + ".diffuse").c_str(), diffuse);
    this->shader->setFloat((name + ".reflection").c_str(), reflection);
    this->shader->setFloat((name + ".refraction").c_str(), refraction);
    this->shader->setInt((name + ".specular").c_str(), specular);
    this->shader->setFloat((name + ".kd").c_str(), kd);
    this->shader->setFloat((name + ".ks").c_str(), ks);
}


rt_sphere::rt_sphere() : ShaderElement()
{}
rt_sphere::rt_sphere(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_sphere::update()
{
    mat.name = name + ".mat";
    mat.affectShader(this->shader);
    mat.update();

    this->shader->setVector((name + ".obj").c_str(), obj);
    this->shader->setVector((name + ".quat_rotation").c_str(), quat_rotation);
    this->shader->setInt((name + ".textureNum").c_str(), textureNum);
    this->shader->setBool((name + ".hollow").c_str(), hollow);
}


rt_plane::rt_plane() : ShaderElement()
{}
rt_plane::rt_plane(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_plane::update()
{
    mat.name = name + ".mat";
    mat.affectShader(this->shader);
    mat.update();

    this->shader->setVector((name + ".pos").c_str(), pos);
    this->shader->setVector((name + ".normal").c_str(), normal);
}


rt_box::rt_box() : ShaderElement()
{}
rt_box::rt_box(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_box::update()
{
    mat.name = name + ".mat";
    mat.affectShader(this->shader);
    mat.update();

    this->shader->setVector((name + ".quat_rotation").c_str(), quat_rotation);
    this->shader->setVector((name + ".pos").c_str(), pos);
    this->shader->setVector((name + ".form").c_str(), form);
    this->shader->setInt((name + ".textureNum").c_str(), textureNum);
}


rt_ring::rt_ring() : ShaderElement()
{}
rt_ring::rt_ring(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_ring::update()
{
    mat.name = name + ".mat";
    mat.affectShader(this->shader);
    mat.update();

    this->shader->setVector((name + ".quat_rotation").c_str(), quat_rotation);
    this->shader->setVector((name + ".pos").c_str(), pos);
    this->shader->setInt((name + ".textureNum").c_str(), textureNum);
    this->shader->setFloat((name + ".r1").c_str(), r1);
    this->shader->setFloat((name + ".r2").c_str(), r2);
}


rt_surface::rt_surface() : ShaderElement()
{}
rt_surface::rt_surface(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_surface::update()
{
    mat.name = name + ".mat";
    mat.affectShader(this->shader);
    mat.update();

    this->shader->setVector((name + ".quat_rotation").c_str(), quat_rotation);
    this->shader->setVector((name + ".pos").c_str(), pos);
    this->shader->setVector((name + ".v_min").c_str(), v_min);
    this->shader->setVector((name + ".v_max").c_str(), v_max);
    this->shader->setFloat((name + ".a").c_str(), a);
    this->shader->setFloat((name + ".b").c_str(), b);
    this->shader->setFloat((name + ".c").c_str(), c);
    this->shader->setFloat((name + ".d").c_str(), d);
    this->shader->setFloat((name + ".e").c_str(), e);
    this->shader->setFloat((name + ".f").c_str(), f);
}


rt_torus::rt_torus() : ShaderElement()
{}
rt_torus::rt_torus(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_torus::update()
{
    mat.name = name + ".mat";
    mat.affectShader(this->shader);
    mat.update();

    this->shader->setVector((name + ".quat_rotation").c_str(), quat_rotation);
    this->shader->setVector((name + ".pos").c_str(), pos);
    this->shader->setVector((name + ".form").c_str(), form);
}


rt_light_direct::rt_light_direct() : ShaderElement()
{}
rt_light_direct::rt_light_direct(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_light_direct::update()
{
    this->shader->setVector((name + ".direction").c_str(), direction);
    this->shader->setVector((name + ".color").c_str(), color);
    this->shader->setFloat((name + ".intensity").c_str(), intensity);

}


rt_light_point::rt_light_point() : ShaderElement()
{}
rt_light_point::rt_light_point(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_light_point::update()
{
    this->shader->setVector((name + ".pos").c_str(), pos);
    this->shader->setVector((name + ".color").c_str(), color);
    this->shader->setFloat((name + ".intensity").c_str(), intensity);
    this->shader->setFloat((name + ".linear_k").c_str(), linear_k);
    this->shader->setFloat((name + ".quadratic_k").c_str(), quadratic_k);
}


rt_scene::rt_scene() : ShaderElement()
{}
rt_scene::rt_scene(std::string name, std::shared_ptr<Shader> shader) : ShaderElement(name, shader)
{}
void rt_scene::update()
{
    this->shader->setVector((name + ".camera_pos").c_str(), camera_pos);
    this->shader->setMat4((name + ".proj_matrix").c_str(), proj_matrix);
    this->shader->setMat4((name + ".mv_matrix").c_str(), mv_matrix);
    this->shader->setVector((name + ".bg_color").c_str(), bg_color);
    this->shader->setInt((name + ".canvas_width").c_str(), canvas_width);
    this->shader->setInt((name + ".canvas_height").c_str(), canvas_height);
    this->shader->setInt((name + ".reflect_depth").c_str(), reflect_depth);
}
