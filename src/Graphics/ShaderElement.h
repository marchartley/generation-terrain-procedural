#ifndef SHADERELEMENT_H
#define SHADERELEMENT_H

#include "DataStructure/Vector3.h"

class ShaderElement
{
public:
    ShaderElement();
};

class ColoredElement : public ShaderElement
{
public:
    ColoredElement();
    ColoredElement(float ambiant[4], float diffuse[4], float specular[4]);

    float ambiant[4];
    float diffuse[4];
    float specular[4];
};

class LightSource : public ColoredElement
{
public:
    LightSource();
    LightSource(float ambiant[4], float diffuse[4], float specular[4]);
};
class PositionalLight : public LightSource
{
public:
    PositionalLight();
    PositionalLight(float ambiant[4], float diffuse[4], float specular[4], Vector3 pos);

    Vector3 position;
};

class Material : public ColoredElement
{
public:
    Material();
    Material(float ambiant[4], float diffuse[4], float specular[4], float shininess);

    float shininess;
};

#endif // SHADERELEMENT_H
