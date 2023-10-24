#include "ComputeShader.h"

ComputeShader::ComputeShader()
    : ComputeShader("")
{

}

ComputeShader::ComputeShader(std::string shaderPath)
    : Shader("", "", ""), computeShaderFilename(shaderPath)
{

}


void ComputeShader::compileShadersFromSource()
{
#if useModernOpenGL || !useModernOpenGL
    this->programID = GlobalsGL::f()->glCreateProgram();
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
#endif
}

void ComputeShader::apply(int nbWorkersX, int nbWorkersY, int nbWorkersZ) {
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
