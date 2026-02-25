#ifndef MATRIX3SERIALIZER_H
#define MATRIX3SERIALIZER_H

#include "DataStructure/Matrix3.h"
#include "Utils/json.h"

template <class Json>
void to_json(Json& json, const GridF& grid);
template <class Json>
void to_json(Json& json, const GridV3& grid);

template <class Json>
void from_json(const Json& json, GridF& grid);
template <class Json>
void from_json(const Json& json, GridV3& grid);




template <class Json>
void to_json(Json& json, const GridF &grid)
{
    json = stringifyGridF(grid, false);
}

template <class Json>
void to_json(Json& json, const GridV3 &grid)
{
    json = stringifyGridV3(grid, false);
}

template <class Json>
void from_json(const Json& json, GridF &grid)
{
    grid = loadGridF(json, false);
}

template <class Json>
void from_json(const Json& json, GridV3 &grid)
{
    grid = loadGridV3(json, false);
}

#endif // MATRIX3SERIALIZER_H
