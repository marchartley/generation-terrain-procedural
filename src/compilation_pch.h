#ifndef TERRAIN_COMPILATION_PCH_H
#define TERRAIN_COMPILATION_PCH_H

// Stable, widely shared dependencies only. In particular, do not include STB
// or other implementation-macro headers: main.cpp enables their definitions.
#ifdef __cplusplus
#include <QGLViewer/vec.h>
#include <QGLViewer/quaternion.h>
#include "Utils/json.h"
#include "Utils/FastNoiseLit.h"
#endif

#endif
