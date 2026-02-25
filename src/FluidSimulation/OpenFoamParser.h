#ifndef OPENFOAMPARSER_H
#define OPENFOAMPARSER_H

#include "DataStructure/Matrix3.h"

class OpenFoamParser
{
public:
    OpenFoamParser();

    static GridV3 parseSimulation(const std::string& foldername);

    static std::string createSimulationFile(const std::string& foldername, const GridF &boundaries);
};

#endif // OPENFOAMPARSER_H
