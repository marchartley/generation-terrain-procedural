#ifndef HOTRELOADFILE_H
#define HOTRELOADFILE_H

#include <string>
#include <vector>
#include <filesystem>
#include <functional>
#include <fstream>
#include "Utils/Signals.h"

class HotreloadFile
{
public:
    HotreloadFile();
    HotreloadFile(const std::string& path);
    HotreloadFile(const std::string& path, const std::function<void(const std::string&)>& onChangeFunc);

    std::string path;

    bool check(bool verbose = true);
    // void onChange(const std::function<void(const std::string&)>& func);
    std::string read();

    DECLARE_EVENT(Change, (const std::string& content), (content))
    // std::vector<std::function<void(const std::string&)>> onChangeCallbacks;
    std::filesystem::file_time_type lastChange;
};

#endif // HOTRELOADFILE_H
