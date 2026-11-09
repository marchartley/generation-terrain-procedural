# Reuse the application's sources, replacing its interactive main with a smoke test.
PROJECT_ROOT = $$clean_path($$PWD/..)
include($$PROJECT_ROOT/generation-terrain-procedural.pro)

SMOKE_SOURCES =
for(source, SOURCES): SMOKE_SOURCES += $$absolute_path($$source, $$PROJECT_ROOT)
SOURCES = $$SMOKE_SOURCES
SOURCES -= $$PROJECT_ROOT/src/main.cpp
SOURCES += $$PWD/opengl_lifecycle_smoke.cpp

SMOKE_HEADERS =
for(header, HEADERS): SMOKE_HEADERS += $$absolute_path($$header, $$PROJECT_ROOT)
HEADERS = $$SMOKE_HEADERS
SMOKE_RESOURCES =
for(resource, RESOURCES): SMOKE_RESOURCES += $$absolute_path($$resource, $$PROJECT_ROOT)
RESOURCES = $$SMOKE_RESOURCES
SMOKE_INCLUDES =
for(directory, INCLUDEPATH): SMOKE_INCLUDES += $$absolute_path($$directory, $$PROJECT_ROOT)
INCLUDEPATH = $$SMOKE_INCLUDES
TARGET = opengl-lifecycle-smoke
