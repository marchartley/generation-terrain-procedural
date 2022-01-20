#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include "DataStructure/Vector3.h"

std::vector<std::string> split(std::string str, char c = ' ');
bool makedir(std::string path);
Vector3 HSVtoRGB(float H, float S,float V);

#endif // UTILS_H
