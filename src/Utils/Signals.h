#ifndef SIGNALS_H
#define SIGNALS_H
#include <functional>
#include <vector>
#include <cstddef>
#include <utility>


// Here is an easy to use Callback framework. Example usage (in a class):
// DECALRE_EVENT(MouseMoved, (int pos_x, int pos_y), (pos_x, pos_y))
// Note: second parameter's names MUST match third parameter's names.
#define DECLARE_EVENT(EVENT_NAME, PARAMS, ARGS)                             \
private:                                                                    \
    std::vector<std::function<void PARAMS>> on##EVENT_NAME##Callbacks;          \
public:                                                                     \
    void setOn##EVENT_NAME(std::function<void PARAMS> func)         \
{                                                                       \
    on##EVENT_NAME##Callbacks.push_back(std::move(func));                   \
}                                                                       \
                                                                            \
    void emit##EVENT_NAME PARAMS                                           \
{                                                                       \
    for (auto& func : on##EVENT_NAME##Callbacks)                            \
        func ARGS;                                                      \
}


// More complex Signal/Slot system, but code completion may not recognize parameters needed, so I do not use it....
// Example of usage :
// class MyClass {
// public:
//   Signal<> onUpdate;
//   Signal<int, int> onMouseMoved;
//
//   void foo() {
//      onUpdate.emit();
//      onMouseMoved(5, 10);
//   }
// };
// int main() {
//   MyClass obj;
//   obj.onMouseMoved.connect([](int x, int y) { std::cout << x << " " << y << std::endl; });
// }
template<typename... Args>
class Signal
{
public:
    using Slot = std::function<void(Args...)>;
    using ConnectionId = std::size_t;

    ConnectionId connect(Slot slot)
    {
        slots_.push_back({nextId_, std::move(slot), true});
        return nextId_++;
    }

    void disconnect(ConnectionId id)
    {
        for (auto& s : slots_)
        {
            if (s.id == id)
            {
                s.active = false;
                break;
            }
        }
    }

    void emit(Args... args)
    {
        for (auto& s : slots_)
        {
            if (s.active)
                s.slot(args...);
        }
    }

private:
    struct Entry
    {
        ConnectionId id;
        Slot slot;
        bool active;
    };

    std::vector<Entry> slots_;
    ConnectionId nextId_ = 0;
};

#endif // SIGNALS_H
