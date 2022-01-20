#include "Graphics/ShaderElement.h"

ShaderElement::ShaderElement()
{

}

ColoredElement::ColoredElement()
{
    float defaultAmbiant[4] = {.8, .8, .8, 1.0};
    float defaultDiffuse[4] = {.5, .5, .5, 1.0};
    float defaultSpecular[4] = {.5, .5, .5, 1.0};
    ColoredElement(defaultAmbiant, defaultDiffuse, defaultSpecular);
}
ColoredElement::ColoredElement(float ambiant[4], float diffuse[4], float specular[4])
{
    for (int i = 0; i < 4; i++)
    {
        this->ambiant[i] = ambiant[i];
        this->diffuse[i] = diffuse[i];
        this->specular[i] = specular[i];
    }
}
LightSource::LightSource()
{

}
LightSource::LightSource(float ambiant[4], float diffuse[4], float specular[4])
    : ColoredElement(ambiant, diffuse, specular)
{

}
PositionalLight::PositionalLight()
{

}
PositionalLight::PositionalLight(float ambiant[4], float diffuse[4], float specular[4], Vector3 pos)
    : LightSource(ambiant, diffuse, specular), position(pos)
{

}
Material::Material()
{

}
Material::Material(float ambiant[4], float diffuse[4], float specular[4], float shininess)
    : ColoredElement(ambiant, diffuse, specular), shininess(shininess)
{

}
