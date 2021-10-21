#ifndef VIEWER_H
#define VIEWER_H

#include "Grid.h"
#include "VoxelGrid.h"
#include <QGLViewer/qglviewer.h>

enum ViewerMode {
    GRID_MODE = 0,
    VOXEL_MODE = 1
};

class Viewer : public QGLViewer {
public:
  Viewer(Grid g);
  Viewer(VoxelGrid g);

protected:
  virtual void init();
  virtual void draw();

  ViewerMode mode;

private:
  Grid grid;
  VoxelGrid voxelGrid;
};


#endif // VIEWER_H
