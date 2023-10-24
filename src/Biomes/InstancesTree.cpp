#include "InstancesTree.h"

InstancesTree::InstancesTree()
{

}

InstancesTree::InstancesTree(std::shared_ptr<BiomeInstance> rootBiome)
    : root(rootBiome)
{

}

int InstancesTree::getHeight()
{
}

std::vector<std::shared_ptr<BiomeInstance> > InstancesTree::getNodesAtHeight(int height, bool ignorePriorityOffsets)
{

}

std::vector<std::shared_ptr<BiomeInstance> > InstancesTree::extractAllNodes()
{

}
