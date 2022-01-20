#include "DebugShader.h"

DebugShader::DebugShader()
    : DebugShader(nullptr, nullptr, nullptr)
{

}
DebugShader::DebugShader(const char* vertexShaderFilename)
    : DebugShader(vertexShaderFilename, nullptr, nullptr)
{

}
DebugShader::DebugShader(const char* vertexShaderFilename, const char* fragmentShaderFilename)
    : DebugShader(vertexShaderFilename, fragmentShaderFilename, nullptr)
{

}
DebugShader::DebugShader(const char* vertexShaderFilename, const char* fragmentShaderFilename,
       const char* geometryShaderFilename)
    : Shader(vertexShaderFilename, fragmentShaderFilename, geometryShaderFilename)
{
}
