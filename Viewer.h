#ifndef VIEWER_H
#define VIEWER_H

#include "Grid.h"
#include "VoxelGrid.h"
#include <QGLViewer/qglviewer.h>

enum ViewerMode {
    GRID_MODE  = 0b0001,
    VOXEL_MODE = 0b0010,
    FILL_MODE  = 0b0100,
    WIRE_MODE  = 0b1000
};

class Viewer : public QGLViewer {
public:
    Viewer(Grid* grid, VoxelGrid* voxelGrid, int mode);
  Viewer(Grid* g);
  Viewer(VoxelGrid* g);

  void setMode(int newMode) { this->viewerMode = newMode; }

protected:
  virtual void init();
  virtual void draw();

  int viewerMode;

private:
  Grid* grid;
  VoxelGrid* voxelGrid;
};


#endif // VIEWER_H
