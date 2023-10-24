#ifndef DEBUGSHADER_H
#define DEBUGSHADER_H

#include "Graphics/Shader.h"

class DebugShader : public Shader
{
public:
    DebugShader();
    DebugShader(const char* vertexShaderFilename);
    DebugShader(const char* vertexShaderFilename, const char* fragmentShaderFilename);
    DebugShader(const char* vertexShaderFilename, const char* fragmentShaderFilename,
                const char* geometryShaderFilename);
};

#endif // DEBUGSHADER_H
