#ifndef IMPLICITTREE_H
#define IMPLICITTREE_H

#include "TerrainGen/ImplicitPatch.h"

// Preserve the shared root's identity without slicing a binary, unary or 2D
// source root into an ordinary group. Assignment clones before clearing state.
inline void replaceImplicitRoot(ImplicitNaryOperator& root, const ImplicitPatch& source)
{
    if (typeid(root) != typeid(ImplicitNaryOperator))
        throw std::invalid_argument("Shared terrain replacement requires an ordinary implicit root");
    if (typeid(source) == typeid(ImplicitNaryOperator)) {
        root = static_cast<const ImplicitNaryOperator&>(source);
    } else {
        ImplicitNaryOperator wrapper;
        wrapper.addChild(source.copy());
        wrapper.scaling = source.scaling;
        wrapper.translation = source.translation;
        wrapper.name = source.name;
        root = wrapper;
    }
    root.used_json_filename.clear();
    root.referenceProperties = nlohmann::json::object();
    root.updateCache();
}

#endif
