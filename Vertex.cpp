#include "Vertex.h"
#include <QGLViewer/qglviewer.h>

Vertex::Vertex() : Vertex(0.0, 0.0, 0.0, 0.0)
{

}
Vertex::Vertex(float x, float y, float z, float isosurface)
    : Vector3(x, y, z), isosurface(isosurface)
{

}
Vertex::Vertex(Vector3 v, float isosurface)
    : Vector3(v), isosurface(isosurface) {

}

void Vertex::display() {
    if (isosurface > 0)
        glColor3f(1., 1., 1.);
    else
        glColor3f(0., 0., 0.);
    glBegin(GL_QUADS);
    glVertex3f(- 0.01, - 0.01, - 0.01);
    glVertex3f(- 0.01, - 0.01, + 0.01);
    glVertex3f(- 0.01, + 0.01, + 0.01);
    glVertex3f(- 0.01, + 0.01, - 0.01);

    glVertex3f(- 0.01, - 0.01, - 0.01);
    glVertex3f(+ 0.01, - 0.01, - 0.01);
    glVertex3f(+ 0.01, + 0.01, - 0.01);
    glVertex3f(- 0.01, + 0.01, - 0.01);

    glVertex3f(+ 0.01, - 0.01, - 0.01);
    glVertex3f(+ 0.01, - 0.01, + 0.01);
    glVertex3f(- 0.01, - 0.01, + 0.01);
    glVertex3f(- 0.01, - 0.01, - 0.01);
    glEnd();
}
