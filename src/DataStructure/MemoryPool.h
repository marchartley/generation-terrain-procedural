#ifndef MEMORYPOOL_H
#define MEMORYPOOL_H

#include <list>

template <typename T>
class NodeBlock {
public:
    NodeBlock();
    ~NodeBlock();

    static const size_t BLOCK_SIZE = 256; // or another suitable number
    T* nodes;
    size_t usedNodes = 0;
};

template <typename T>
NodeBlock<T>::NodeBlock() {
    this->nodes = new T[NodeBlock::BLOCK_SIZE];
}

template <typename T>
NodeBlock<T>::~NodeBlock() {
    delete[] nodes;
}

template <typename T>
class MemoryPool {
private:
    std::list<NodeBlock<T>*> blocks;  // Change to a list of pointers
    NodeBlock<T>* currentBlock = nullptr;

public:
    MemoryPool();
    ~MemoryPool();

    T* allocate();
    T* parallelAllocate();

private:
    void allocateBlock();
    void parallelAllocateBlock();
};

template <typename T>
MemoryPool<T>::MemoryPool() {
    allocateBlock();
}

template <typename T>
MemoryPool<T>::~MemoryPool() {
    for (NodeBlock<T>* block : blocks) {
        delete block;  // Delete each allocated block
    }
}

template <typename T>
T* MemoryPool<T>::allocate() {
    if (currentBlock->usedNodes == NodeBlock<T>::BLOCK_SIZE) {
        allocateBlock();
    }
    return &(currentBlock->nodes[currentBlock->usedNodes++]);
}

template <typename T>
T* MemoryPool<T>::parallelAllocate() {
    T* res;
    #pragma omp critical
    {
        if (currentBlock->usedNodes == NodeBlock<T>::BLOCK_SIZE) {
            allocateBlock();
        }
        res = &(currentBlock->nodes[currentBlock->usedNodes++]);
    }
    return res;
}

template <typename T>
void MemoryPool<T>::allocateBlock() {
    currentBlock = new NodeBlock<T>();  // Use new to allocate
    blocks.emplace_back(currentBlock);
}

template <typename T>
void MemoryPool<T>::parallelAllocateBlock() {
    #pragma omp critical
    {
        currentBlock = new NodeBlock<T>();  // Use new to allocate
        blocks.emplace_back(currentBlock);
    }
}

#endif // MEMORYPOOL_H
