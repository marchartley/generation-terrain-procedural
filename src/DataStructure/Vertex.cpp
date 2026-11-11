#include "DataStructure/Vertex.h"
#include <QGLViewer/qglviewer.h>
#include <QOpenGLContext>
#include <stdexcept>

Vertex::Vertex() : Vertex(0.0, 0.0, 0.0, -1.0)
{

}
Vertex::Vertex(float x, float y, float z, float isosurface)
    : Vector3(x, y, z), isosurface(isosurface)
{

}
Vertex::Vertex(const Vector3& v, float isosurface)
    : Vector3(v), isosurface(isosurface) {

}

void Vertex::display() {
    if (!QOpenGLContext::currentContext())
        throw std::logic_error("Vertex display requires a current OpenGL context");
    glPushAttrib(GL_CURRENT_BIT);
    glPushMatrix();
    glTranslatef(x(), y(), z());
    float mini_change = 0.01; //(rand() % 100)/3000.0;
    if (isosurface > 0)
        glColor3f(1., 1., 1.);
    else
        glColor3f(0., 0., 0.);
    glBegin(GL_QUADS);
    glVertex3f(- mini_change, - mini_change, - mini_change);
    glVertex3f(- mini_change, - mini_change, + mini_change);
    glVertex3f(- mini_change, + mini_change, + mini_change);
    glVertex3f(- mini_change, + mini_change, - mini_change);

    glVertex3f(- mini_change, - mini_change, - mini_change);
    glVertex3f(+ mini_change, - mini_change, - mini_change);
    glVertex3f(+ mini_change, + mini_change, - mini_change);
    glVertex3f(- mini_change, + mini_change, - mini_change);

    glVertex3f(+ mini_change, - mini_change, - mini_change);
    glVertex3f(+ mini_change, - mini_change, + mini_change);
    glVertex3f(- mini_change, - mini_change, + mini_change);
    glVertex3f(- mini_change, - mini_change, - mini_change);
    glEnd();
    glPopMatrix();
    glPopAttrib();
}
