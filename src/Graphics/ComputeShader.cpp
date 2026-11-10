#include "ComputeShader.h"

ComputeShader::ComputeShader()
    : ComputeShader("")
{

}

ComputeShader::ComputeShader(const std::string& shaderPath)
    : Shader("", "", ""), computeShaderFilename(shaderPath)
{

}


void ComputeShader::compileShadersFromSource()
{
    GlobalsGL::ContextScope contextScope;
    ensureResources();
    if (programID) GlobalsGL::f()->glDeleteProgram(programID);
    gpuResources->programs.clear();
#if useModernOpenGL || !useModernOpenGL
    this->programID = GlobalsGL::f()->glCreateProgram();
    gpuResources->programs.push_back(programID);
    if (computeShaderFilename != "")
    {
        std::string content = Shader::readShaderSource(computeShaderFilename);
        if (!content.empty()) {
            this->cShader = GlobalsGL::f()->glCreateShader(GL_COMPUTE_SHADER);
            const char* src = content.c_str();
            GlobalsGL::f()->glShaderSource(this->cShader, 1, &src, NULL);
            GlobalsGL::f()->glCompileShader(this->cShader);
            GlobalsGL::f()->glAttachShader(this->programID, this->cShader);
            GlobalsGL::printShaderErrors(this->cShader);
        } else {
            computeShaderFilename = "";
        }
    }

    GlobalsGL::f()->glLinkProgram(this->programID);
    invalidateProgramCache();
    GLint linked = GL_FALSE;
    GlobalsGL::f()->glGetProgramiv(programID, GL_LINK_STATUS, &linked);
    programLinked = linked == GL_TRUE;
    GlobalsGL::printProgramErrors(programID);
    if (cShader >= 0) GlobalsGL::f()->glDeleteShader(cShader);
    cShader = -1;
#endif
}

void ComputeShader::apply(int nbWorkersX, int nbWorkersY, int nbWorkersZ) {
    GlobalsGL::ContextScope contextScope;
    if (this->use()) {
        if (nbWorkersX < 0)
            GlobalsGL::f()->glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &nbWorkersX);
        if (nbWorkersY < 0)
            GlobalsGL::f()->glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &nbWorkersY);
        if (nbWorkersZ < 0)
            GlobalsGL::f()->glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &nbWorkersZ);
        GlobalsGL::f()->glDispatchCompute(nbWorkersX, nbWorkersY, nbWorkersZ);
    }
}
