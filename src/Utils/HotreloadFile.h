#ifndef HOTRELOADFILE_H
#define HOTRELOADFILE_H

#include <string>
#include <vector>
#include <filesystem>
#include <functional>
#include <fstream>

class HotreloadFile
{
public:
    HotreloadFile();
    HotreloadFile(std::string path);
    HotreloadFile(std::string path, const std::function<void(std::string)>& onChangeFunc);

    std::string path;

    bool check(bool verbose = true);
    void onChange(const std::function<void(std::string)>& func);
    std::string read();

    std::vector<std::function<void(std::string)>> onChangeCallbacks;
    std::filesystem::file_time_type lastChange;
};

#endif // HOTRELOADFILE_H
