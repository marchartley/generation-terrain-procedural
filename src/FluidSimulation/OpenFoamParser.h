#ifndef OPENFOAMPARSER_H
#define OPENFOAMPARSER_H

#include "DataStructure/Matrix3.h"

class OpenFoamParser
{
public:
    OpenFoamParser();

    static GridV3 parseSimulation(std::string foldername);
};

#endif // OPENFOAMPARSER_H
