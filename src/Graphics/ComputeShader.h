#ifndef COMPUTESHADER_H
#define COMPUTESHADER_H

#include "Shader.h"

class ComputeShader : public Shader
{
public:
    ComputeShader();
    ComputeShader(std::string shaderPath);

    void apply(int nbWorkersX = -1, int nbWorkersY = -1, int nbWorkersZ = -1);

    void compileShadersFromSource();

    std::string computeShaderFilename = "";
    int cShader = -1;
};

#endif // COMPUTESHADER_H
