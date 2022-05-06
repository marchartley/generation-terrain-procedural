#include "RewritableFile.h"

RewritableFile::RewritableFile()
{
}

RewritableFile::RewritableFile(std::string filename)
    : filename(filename)
{
}

void RewritableFile::clearFile()
{
    stream.close();
    stream.open(this->filename, std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
}
