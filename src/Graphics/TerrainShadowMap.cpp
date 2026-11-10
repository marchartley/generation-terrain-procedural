#include "Graphics/TerrainShadowMap.h"
#include "Graphics/RenderMatrices.h"
#include <limits>
#include <stdexcept>

namespace {
// Every state changed by a depth pass is restored, including Qt's framebuffer.
struct ShadowState {
    GLint draw, read, viewport[4], active, depthFunc, polygon[2];
    GLboolean depthMask, depth, cull, blend, scissor, discard;
    GLdouble clearDepth;
    ShadowState() {
        glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &draw);
        glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &read);
        glGetIntegerv(GL_VIEWPORT, viewport);
        glGetIntegerv(GL_ACTIVE_TEXTURE, &active);
        glGetIntegerv(GL_DEPTH_FUNC, &depthFunc);
        glGetIntegerv(GL_POLYGON_MODE, polygon);
        glGetBooleanv(GL_DEPTH_WRITEMASK, &depthMask);
        glGetDoublev(GL_DEPTH_CLEAR_VALUE, &clearDepth);
        depth=glIsEnabled(GL_DEPTH_TEST); cull=glIsEnabled(GL_CULL_FACE);
        blend=glIsEnabled(GL_BLEND); scissor=glIsEnabled(GL_SCISSOR_TEST);
        discard=glIsEnabled(GL_RASTERIZER_DISCARD);
    }
    ~ShadowState() {
        auto* gl = GlobalsGL::f();
        gl->glBindFramebuffer(GL_DRAW_FRAMEBUFFER, draw);
        gl->glBindFramebuffer(GL_READ_FRAMEBUFFER, read);
        glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
        gl->glActiveTexture(active);
        glDepthFunc(depthFunc); glDepthMask(depthMask); glClearDepth(clearDepth);
        glPolygonMode(GL_FRONT, polygon[0]); glPolygonMode(GL_BACK, polygon[1]);
        for (auto state : {std::pair<GLenum,GLboolean>{GL_DEPTH_TEST,depth}, {GL_CULL_FACE,cull},
                          {GL_BLEND,blend}, {GL_SCISSOR_TEST,scissor}, {GL_RASTERIZER_DISCARD,discard}})
            if (state.second) glEnable(state.first); else glDisable(state.first);
    }
};
}

GridF TerrainShadowMap::lightMatrix(const Vector3& light, const Vector3& extent)
{
    // Fit the perspective frustum to the actual terrain, with double precision
    // while multiplying matrices: the application's lights are far from the mesh.
    const Vec3<double> size(std::max(double(extent.x()),.01), std::max(double(extent.y()),.01), std::max(double(extent.z()),.01));
    const Vec3<double> center = size * .5;
    Vec3<double> eye(light.x(),light.y(),light.z());
    if (RenderMatrices::length(eye-center) < .001) eye += Vec3<double>(0,0,.01);
    const auto direction = RenderMatrices::unit(center-eye);
    const auto up = std::abs(direction.z()) > .99 ? Vec3<double>(0,1,0) : Vec3<double>(0,0,1);
    const auto view = RenderMatrices::lookAt(eye, center, up);
    double nearZ=std::numeric_limits<double>::max(), farZ=0;
    double left=nearZ, bottom=nearZ, right=-nearZ, top=-nearZ;
    bool crossesEye = false;
    for (int i=0; i<8; ++i) {
        const Vec3<double> corner((i&1)?size.x():0, (i&2)?size.y():0, (i&4)?size.z():0);
        Vec3<double> p;
        for(int row=0;row<3;++row) {
            p.data()[row]=view(row,3);
            for(int col=0;col<3;++col) p.data()[row]+=view(row,col)*corner.data()[col];
        }
        const double distance = -p.z();
        if (distance <= .001) { crossesEye=true; continue; }
        nearZ=std::min(nearZ,distance); farZ=std::max(farZ,distance);
        left=std::min(left,p.x()/distance); right=std::max(right,p.x()/distance);
        bottom=std::min(bottom,p.y()/distance); top=std::max(top,p.y()/distance);
    }
    if (crossesEye) {
        // A positional light inside the terrain cannot cover all directions
        // with one perspective map. The application uses external lights.
        // Return a wide forward frustum; out-of-frustum receivers stay lit.
        nearZ=.001; farZ=RenderMatrices::length(size)*2+RenderMatrices::length(eye-center);
        left=bottom=-10; right=top=10;
    }
    const double margin = std::max(.01, RenderMatrices::length(size)*.002);
    nearZ=std::max(.001,nearZ-margin); farZ=std::max(nearZ+.01,farZ+margin);
    const double padX=std::max(1e-8,(right-left)*.01), padY=std::max(1e-8,(top-bottom)*.01);
    return RenderMatrices::floats(RenderMatrices::multiply(
        RenderMatrices::frustum((left-padX)*nearZ,(right+padX)*nearZ,
                               (bottom-padY)*nearZ,(top+padY)*nearZ,nearZ,farZ), view));
}

void TerrainShadowMap::update(Mesh& mesh, const Vector3& extent, size_t terrainRevision)
{
    GlobalsGL::ContextScope scope;
    auto shader=mesh.shader;
    if (!shader || !shader->use()) return;
    auto* gl=GlobalsGL::f();
    shader->setInt("shadowMapTex",9);
    ShadowState state;
    if (!resources || resources->released || resources->owner!=QOpenGLContext::currentContext()) {
        resources=GlobalsGL::Resources::create(); key.clear();
        gl->glGenTextures(1,&texture); resources->textures.push_back(texture);
        gl->glActiveTexture(GL_TEXTURE9);
        glBindTexture(GL_TEXTURE_2D_ARRAY,texture);
        gl->glTexImage3D(GL_TEXTURE_2D_ARRAY,0,GL_DEPTH_COMPONENT32F,resolution,resolution,6,0,GL_DEPTH_COMPONENT,GL_FLOAT,nullptr);
        glTexParameteri(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_BORDER);
        const float border[]={1,1,1,1}; glTexParameterfv(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_BORDER_COLOR,border);
        glTexParameteri(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_COMPARE_MODE,GL_NONE);
        glTexParameteri(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_COMPARE_FUNC,GL_LEQUAL);
        gl->glGenFramebuffers(1,&framebuffer); resources->framebuffers.push_back(framebuffer);
    }
    std::vector<unsigned char> nextKey;
    auto append=[&](const void* data,size_t size) {
        const auto* bytes=static_cast<const unsigned char*>(data); nextKey.insert(nextKey.end(),bytes,bytes+size);
    };
    append(&shader->programID,sizeof(shader->programID));
    append(&shader->volumeUploadCount,sizeof(shader->volumeUploadCount));
    append(&terrainRevision,sizeof(terrainRevision));
    const auto pointCount=mesh.vertexArray.size(); append(&pointCount,sizeof(pointCount));
    const float dims[]={extent.x(),extent.y(),extent.z()}; append(dims,sizeof(dims));
    // All uniforms affecting terrain position, topology or clipping. Camera and
    // material uniforms deliberately do not invalidate the maps.
    for (const auto& entry : std::vector<std::pair<std::string,int>>{
        {"min_vertice_positions",3},{"max_vertice_positions",3},{"scale",3},{"rotation",3},
        {"offsetX",1},{"offsetY",1},{"offsetZ",1},{"heightFactor",1},{"maxHeight",1},
        {"clipPlanePosition",3},{"clipPlaneDirection",3},{"min_isolevel",1},{"max_isolevel",1}}) {
        float value[4]={}; const GLint location=shader->uniformLocation(entry.first);
        if (location>=0) gl->glGetUniformfv(shader->programID,location,value);
        append(value,sizeof(float)*entry.second);
    }
    for (const char* name : {"clipPlaneActive","useMarchingCubes","voxels_displayed_on_borders"}) {
        GLint value=0, location=shader->uniformLocation(name);
        if (location>=0) gl->glGetUniformiv(shader->programID,location,&value);
        append(&value,sizeof(value));
    }
    std::array<Vector3,6> lights;
    std::array<bool,6> castsShadows;
    GLint activeCount=6;
    const auto countLocation=shader->uniformLocation("terrainLightCount");
    if(countLocation>=0) gl->glGetUniformiv(shader->programID,countLocation,&activeCount);
    append(&activeCount,sizeof(activeCount));
    for (int i=0;i<6;++i) {
        lights[i]=Vector3(0,0,0);
        const auto location=shader->uniformLocation("lights["+std::to_string(i)+"].position");
        if (location>=0) gl->glGetUniformfv(shader->programID,location,lights[i].data());
        append(lights[i].data(),sizeof(float)*3);
        GLint enabled=0;
        const auto castLocation=shader->uniformLocation("lights["+std::to_string(i)+"].castsShadows");
        if(castLocation>=0) gl->glGetUniformiv(shader->programID,castLocation,&enabled);
        castsShadows[i]=enabled!=0 && i<activeCount;
        append(&castsShadows[i],sizeof(bool));
    }
    if (key!=nextKey) {
        // Unbind the depth array while it is attached: never sample a render target.
        gl->glActiveTexture(GL_TEXTURE9); glBindTexture(GL_TEXTURE_2D_ARRAY,0);
        gl->glBindFramebuffer(GL_FRAMEBUFFER,framebuffer);
        glDrawBuffer(GL_NONE); glReadBuffer(GL_NONE);
        glViewport(0,0,resolution,resolution);
        glEnable(GL_DEPTH_TEST); glDepthFunc(GL_LESS); glDepthMask(GL_TRUE); glClearDepth(1);
        glDisable(GL_BLEND); glDisable(GL_CULL_FACE); glDisable(GL_SCISSOR_TEST);
        glDisable(GL_RASTERIZER_DISCARD); glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        shader->setBool("shadowPass",true);
        for (int i=0;i<6;++i) {
            matrices[i]=lightMatrix(lights[i],Vector3(dims[0],dims[1],dims[2]));
            shader->setMatrix("shadowDepthMatrix",matrices[i].data.data(),4,4);
            gl->glFramebufferTextureLayer(GL_FRAMEBUFFER,GL_DEPTH_ATTACHMENT,texture,0,i);
            if (gl->glCheckFramebufferStatus(GL_FRAMEBUFFER)!=GL_FRAMEBUFFER_COMPLETE) {
                shader->setBool("shadowPass",false);
                throw std::runtime_error("Incomplete terrain shadow framebuffer");
            }
            glClear(GL_DEPTH_BUFFER_BIT);
            if(castsShadows[i]) mesh.display(GL_POINTS);
        }
        shader->setBool("shadowPass",false);
        key=std::move(nextKey); ++renderCount;
    }
    for (int i=0;i<6;++i) shader->setMatrix("shadowMatrices["+std::to_string(i)+"]",matrices[i].data.data(),4,4);
    shader->setInt("shadowMapTex",9);
    gl->glActiveTexture(GL_TEXTURE9); glBindTexture(GL_TEXTURE_2D_ARRAY,texture);
}
