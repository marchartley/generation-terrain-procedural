#ifndef REWRITABLEFILE_H
#define REWRITABLEFILE_H

#include <fstream>
#include <string>
#include <iostream>

class RewritableFile
{
public:
    RewritableFile();
    RewritableFile(std::string filename);

    void clearFile();
    std::ostream& operator<<(std::string text) {
        stream.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
        stream.seekg(0);
        stream << text;
        stream.flush();
        std::cout << text << "\nFile '" << filename << "' is " << (stream.is_open() ? "open" : "closed") << std::endl;
        stream.close();
        return stream;
    }

    std::string filename;
protected:
    std::fstream stream;
};

//template <typename T>
//RewritableFile& operator<<(RewritableFile& file, const T& text) {
//    file.open(file.filename, std::ios_base::out | std::ios_base::ate);
//    file << text;
//    file.close();
//}
#endif // REWRITABLEFILE_H
