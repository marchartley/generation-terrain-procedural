#ifndef GRAPHNODE_H
#define GRAPHNODE_H

#include "DataStructure/Vector3.h"
template <class T>
class GraphNode
{
public:
    GraphNode();
    GraphNode(T value);
    GraphNode(T value, Vector3 pos, int index = 0);

    Vector3 pos;
    T value;
    int index;

    Vector3 privateVector;
    int privateIndex;
};

template<class T>
GraphNode<T>::GraphNode() : GraphNode(T(), Vector3(), 0)
{

}
template<class T>
GraphNode<T>::GraphNode(T value) : GraphNode(value, Vector3(), 0)
{

}
template<class T>
GraphNode<T>::GraphNode(T value, Vector3 pos, int index)
    : value(value), pos(pos), index(index)
{

}

#endif // GRAPHNODE_H
