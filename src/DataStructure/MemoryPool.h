#ifndef MEMORYPOOL_H
#define MEMORYPOOL_H

#include <algorithm>
#include <cstddef>
#include <list>
#include <memory>
#include <mutex>
#include <utility>

template <class T>
class NodeBlock {
public:
  static constexpr size_t BLOCK_SIZE = 256;
  T *nodes;
  size_t usedNodes = 0;

  NodeBlock() : nodes(new T[BLOCK_SIZE]) {}
  ~NodeBlock() { delete[] nodes; }
  NodeBlock(const NodeBlock &other) : nodes(nullptr) {
      // Only copy allocated slots; unused arithmetic slots are uninitialized.
      // Do not delegate construction: a throwing T assignment must have exactly
      // one owner responsible for destroying the temporary block.
      std::unique_ptr<T[]> buffer(new T[BLOCK_SIZE]);
      std::copy_n(other.nodes, other.usedNodes, buffer.get());
      usedNodes = other.usedNodes;
      nodes = buffer.release();
  }
  NodeBlock &operator=(const NodeBlock &other) {
      if (this != &other) {
          NodeBlock copy(other);
          std::swap(nodes, copy.nodes);
          std::swap(usedNodes, copy.usedNodes);
      }
      return *this;
  }
};

// Allocated addresses stay valid until clear/destruction. Concurrent allocations
// must use parallelAllocate; copying/clearing requires exclusive caller access.
template <class T>
class MemoryPool {
    std::list<std::unique_ptr<NodeBlock<T>>> blocks;
    NodeBlock<T>* currentBlock = nullptr;
    std::mutex allocationMutex;

  public:
    MemoryPool() = default;
    ~MemoryPool() = default;
    MemoryPool(const MemoryPool &other) {
        for (const auto &block : other.blocks)
            blocks.push_back(std::make_unique<NodeBlock<T>>(*block));
        if (!blocks.empty())
            currentBlock = blocks.back().get();
    }
    MemoryPool &operator=(const MemoryPool &other) {
        if (this != &other) {
            MemoryPool copy(other);
            blocks.swap(copy.blocks);
            std::swap(currentBlock, copy.currentBlock);
        }
        return *this;
    }
    MemoryPool(MemoryPool &&other) noexcept
        : blocks(std::move(other.blocks)), currentBlock(std::exchange(other.currentBlock, nullptr)) {}
    MemoryPool &operator=(MemoryPool &&other) noexcept {
        if (this != &other) {
            blocks = std::move(other.blocks);
            currentBlock = std::exchange(other.currentBlock, nullptr);
        }
        return *this;
    }
    void clear() noexcept {
        blocks.clear();
        currentBlock = nullptr;
    }
    T *allocate() {
        if (!currentBlock || currentBlock->usedNodes == NodeBlock<T>::BLOCK_SIZE) {
            auto block = std::make_unique<NodeBlock<T>>();
            blocks.push_back(std::move(block));
            currentBlock = blocks.back().get();
        }
        return &currentBlock->nodes[currentBlock->usedNodes++];
    }
    T *parallelAllocate() {
        std::lock_guard<std::mutex> lock(allocationMutex);
        return allocate();
    }
};
#endif
