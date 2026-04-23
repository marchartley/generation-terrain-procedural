#include "HotreloadFile.h"

#include <iostream>
#include <chrono>

HotreloadFile::HotreloadFile()
{

}

HotreloadFile::HotreloadFile(const std::string& path)
    : path(path)
{

}

HotreloadFile::HotreloadFile(const std::string& path, const std::function<void (const std::string&)> &onChangeFunc)
    : HotreloadFile(path)
{
    setOnChange(onChangeFunc);
}

bool HotreloadFile::check(bool verbose)
{
    if (this->path.empty()) {
        std::cerr << "File still empty!" << std::endl;
        return false;
    }

    std::filesystem::file_time_type timeModif = std::filesystem::last_write_time(this->path);
    bool difference = timeModif != lastChange;
    this->lastChange = timeModif;

    if (difference) {
        if (verbose) {
            std::cout << "File " << this->path << " modified." << std::endl;
        }
        if (onChangeCallbacks.size() > 0) {
            std::string content = this->read();
            /*for (auto& callback : onChangeCallbacks) {
                callback(content);
            }*/
            emitChange(content);
        }
    }
    return difference;
}

// void HotreloadFile::onChange(const std::function<void(const std::string&)> &func)
// {
    // onChangeCallbacks.push_back(func);
// }

std::string HotreloadFile::read()
{
    if (this->path.empty()) return "";
    std::ifstream file(this->path);
    std::string fileContent((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    return fileContent;
}
