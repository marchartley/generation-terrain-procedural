#ifndef INSTANCESTREE_H
#define INSTANCESTREE_H

#include "Biomes/BiomeInstance.h"

class InstancesTree
{
public:
    InstancesTree();
    InstancesTree(std::shared_ptr<BiomeInstance> rootBiome);

    int getHeight();
    std::vector<std::shared_ptr<BiomeInstance>> getNodesAtHeight(int height, bool ignorePriorityOffsets = true);
    std::vector<std::shared_ptr<BiomeInstance>> extractAllNodes();

    std::shared_ptr<BiomeInstance> root;
};

#endif // INSTANCESTREE_H
