#include "TerrainModel.h"

TerrainModel::TerrainModel()
{

}

TerrainModel::~TerrainModel()
{

}

float TerrainModel::getHeight(Vector3 pos)
{
    return this->getHeight(pos.x, pos.y);
}

bool TerrainModel::contains(Vector3 v)
{
    return Vector3::isInBox(v, Vector3(), this->getDimensions());
}

bool TerrainModel::contains(float x, float y, float z)
{
    return this->contains(Vector3(x, y, z));
}

size_t TerrainModel::getCurrentHistoryIndex() const
{
    return this->_cachedHistoryIndex; // Not good
}

