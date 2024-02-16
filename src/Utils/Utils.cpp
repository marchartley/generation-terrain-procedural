#include "Utils.h"

#include <sys/stat.h>
#include <thread>
#include <chrono>

#include "Utils/BSpline.h"

std::vector<std::string> split(std::string str, std::string c)
{
    std::vector<std::string> result;
    size_t pos = str.npos;
    do {
        pos = str.rfind(c, pos);
        if (pos != str.npos) {
            std::string sub = str.substr(pos + 1);
            if(sub != "")
                result.insert(result.begin(), sub);
            str = str.substr(0, pos);
        }
    } while(pos != str.npos);
    result.insert(result.begin(), str);
    return result;
}

std::vector<std::string> split(std::string str)
{
    std::vector<std::string> result;
    for (char c : str) {
        result.push_back(std::string(1, c));
    }
    return result;
}

bool checkPathExists(std::string path)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false;
    }
    return true;
}

bool makedir(std::string path)
{
    std::vector<std::string> splitted = split(path, "/");
    int result = 0;
    std::string currentPath = "";
    for (size_t i = 0; i < splitted.size(); i++) {
        currentPath += splitted[i] + '/';
        struct stat info;
        if(stat(currentPath.c_str(), &info) != 0) { // Folder doesn't exist
#ifdef linux
            mode_t prevMode = umask(0011);
            result = mkdir(currentPath.c_str(), 0777); // Create with full permission for linux
            chmod(currentPath.c_str(), 0777);
            umask(prevMode);
#elif _WIN32
            result = mkdir(currentPath.c_str()); // Create for windows
#endif
            if (result != 0)
                return false;
        }
    }
    return true;
}



Vector3 HSVtoRGB(float H, float S,float V){
    H *= 360;
    float s = S;
    float v = V;
    float C = s*v;
    float X = C*(1-abs(fmod(H/60.0, 2)-1));
    float m = v-C;
    float r,g,b;
    if(H >= 0 && H < 60){
        r = C,g = X,b = 0;
    }
    else if(H >= 60 && H < 120){
        r = X,g = C,b = 0;
    }
    else if(H >= 120 && H < 180){
        r = 0,g = C,b = X;
    }
    else if(H >= 180 && H < 240){
        r = 0,g = X,b = C;
    }
    else if(H >= 240 && H < 300){
        r = X,g = 0,b = C;
    }
    else{
        r = C,g = 0,b = X;
    }
    float R = (r+m);
    float G = (g+m);
    float B = (b+m);
    return Vector3(R, G, B);
}

namespace interpolation {
float linear(float x, float _min, float _max) {
    return (x - _min) / (_max - _min);
}
float inv_linear(float x, float _min, float _max) {
    return x * (_max - _min) + _min;
}
float sigmoid(float _x, float lambda, float offset, float _min, float _max) {
    float x = linear(_x, _min, _max);
    float s_0 = 1 / (1 + std::exp(-lambda * (0.f+offset)));
    float s_1 = 1 / (1 + std::exp(-lambda * (1.f+offset)));
    return inv_linear(linear(1 / (1 + std::exp(-lambda * (x+offset))), s_0, s_1), _min, _max);
}
float smooth(float _x, float _min, float _max) {
    float x = linear(_x, _min, _max);
    return inv_linear(3*x*x-2*x*x*x, _min, _max);
}
float quadratic(float _x, float _min, float _max) {
    float x = linear(_x, _min, _max);
    return inv_linear(x * x, _min, _max);
}
float cubic(float _x, float _min, float _max) {
    float x = linear(_x, _min, _max);
    return inv_linear(x * x * x, _min, _max);
}
float cosine(float _x, float _min, float _max) {
    float x = linear(_x, _min, _max);
    float pi = 3.141592;
    return inv_linear((std::cos(x * pi + pi) / 2.f) + 0.5, _min, _max);
}
float binary(float _x, float _min, float _max) {
    float x = linear(_x, _min, _max);
    return inv_linear((x < 0.5 ? 0.f : 1.f), _min, _max);
}

float wyvill(float _x, float _min, float _max)
{
    float x = linear(_x, _min, _max);
    return inv_linear(std::pow(1 - x*x, 3), _min, _max);
}

float fault_distance(float distance, float impactRadius)
{
//    float a = distance / impactRadius;
//    float b = std::pow(a, 2);
//    float c = 1 - b;
//    float d = std::pow(c, 2);
    return (distance < impactRadius ? std::pow(1 - std::pow(distance / impactRadius, 2), 2) : 0);
}

}

std::string toUpper(std::string s)
{
    std::string res = s;
    for (auto& c : res) {
        c = toupper(c);
    }
    return res;
}

std::string toLower(std::string s)
{
    std::string res = s;
    for (auto& c : res) {
        c = tolower(c);
    }
    return res;
}

std::string toCapitalize(std::string s)
{
    std::string res = s;
    bool needCapital = true;
    std::vector<char> replacedBySpace = {'_', '-'};
    for (auto& c : res) {
        if (isIn(c, replacedBySpace)) {
            c = ' ';
            needCapital = true;
        } else {
            c = (needCapital ? toupper(c) : tolower(c));
            needCapital = false;
        }
    }
    return res;
}

std::string getExtension(std::string file)
{
    std::string ext = file.substr(file.find_last_of('.') + 1);
    return ext;
}

float rad2deg(float rad)
{
    return (rad * 180.f) / PI;
}

float deg2rad(float deg)
{
    return (deg * PI) / 180.f;
}

float gaussian(float sigma, float sqrDist) {
    float oneOverSqrt2Pi = 1.f/(2 * M_PI * sigma * sigma);
    float sqrSigma = 2 * sigma * sigma;
    return std::exp(-sqrDist/sqrSigma) * oneOverSqrt2Pi;
}

float gaussian(const Vector3& size, const Vector3& position, float sigma)
{
//    float oneOverSqrt2Pi = 1.f/(2 * M_PI * sigma * sigma);
//    float sqrSigma = 2 * sigma * sigma;
    return gaussian(sigma, (position - (size * .5f)).norm2());
//    float gaussian = std::exp(-position.norm2()/sqrSigma) * oneOverSqrt2Pi;
//    return gaussian;
}

float normalizedGaussian(float sigma, float sqrDist)
{
    float maxValue = gaussian(sigma, 0.f);
    if (maxValue > 0.f)
        return gaussian(sigma, sqrDist) / maxValue;
    return 0.f;
}

float normalizedGaussian(const Vector3& size, const Vector3& position, float sigma)
{
    float maxValue = gaussian(size, size * .5f, sigma);
    if (maxValue > 0.f)
        return gaussian(size, position, sigma) / maxValue;
    return 0.f;
}

float normalDistribution(const Vector3& size, const Vector3& position, float sigma)
{
    float oneOverSqrt2Pi = 1.f/std::sqrt(2 * 3.141592);
    float sqrSigma = 2 * sigma * sigma;
//    position -= (size * .5f);
    float normal = std::exp(-(position - size * .5f).norm2()/(sqrSigma)) * oneOverSqrt2Pi;
    return normal;
}

float normalizedNormalDistribution(const Vector3& size, const Vector3& position, float sigma)
{
    float maxValue = normalDistribution(size, size * .5f, sigma);
    if (maxValue > 0.f)
        return normalDistribution(size, position, sigma) / maxValue;
    return 0.f;
}

bool startsWith(std::string text, std::string needle)
{
    return text.substr(0, needle.size()) == needle;
}

bool endsWith(std::string text, std::string needle)
{
    return text.substr(text.size() - needle.size()) == needle;
}


void displayProgress(float percent, bool displayPercent, int consoleWidth, std::string progressSign)
{
    std::string wait[4] = {"-", "\\", "|", "/"};
    std::cout << "\r[";
    int availableSpaceForSymbols = consoleWidth - 2;
    int nbSymbols = int(percent * float(availableSpaceForSymbols));
    int waitingState = int(((percent * float(availableSpaceForSymbols)) - nbSymbols) * 8) % 4;
    availableSpaceForSymbols -= 1;
    std::string percentValue = to_string_with_precision(percent * 100, 1);
    if (displayPercent)
        availableSpaceForSymbols -= (percentValue.size() + 1);
    for (int i = 0; i < std::min(nbSymbols, availableSpaceForSymbols); i++)
        std::cout << progressSign;
    std::cout << wait[waitingState];
    for (int i = std::min(nbSymbols, availableSpaceForSymbols); i < availableSpaceForSymbols; i++)
        std::cout << " ";
    if (displayPercent)
        std::cout << percentValue << "%";
    std::cout << "]" << std::flush;
}

std::string replaceInString(std::string initial, std::string toReplace, std::string replacing)
{
    std::string replaced = initial;
    while (replaced.find(toReplace.c_str()) != replaced.npos) {
        size_t pos = replaced.find(toReplace.c_str());
        replaced.replace(pos, toReplace.size(), replacing);
    }
    return replaced;
}

std::string getFilename(std::string path)
{
    std::vector<std::string> fullPath = split(path, "/");
    fullPath = split(fullPath.back(), "\\");
    return fullPath.back();
}

std::string simplify(std::string s)
{
    return toLower(replaceInString(replaceInString(s, "_", ""), "-", ""));
}



double timeIt(std::function<void ()> func)
{
    auto start = std::chrono::system_clock::now();
    func();
    auto end = std::chrono::system_clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

void sleep(int milliseconds)
{
    std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds) );
}

std::string showTime(double nanoseconds)
{
    try {
        std::ostringstream oss;
        nanoseconds /= 1000000.f;
        if (nanoseconds != nanoseconds) // NaN
            oss << "0ms";
        else if (nanoseconds < 10) { // < 10ms
            oss << std::setprecision(3) << nanoseconds << "ms";
        } else if (nanoseconds < 1000 * 10) { // < 10 sec
            oss << int(nanoseconds) << "ms";
        } else if (nanoseconds < 1000 * 60 * 10) { // < 10 min
            oss << int(nanoseconds / 1000) << "s";
        } else { // > 10 min
            oss << int(nanoseconds / (1000 * 60)) << "min" << int(nanoseconds) % (1000 * 60) << "s (" + int(nanoseconds / 1000) << "s)";
        }
        return oss.str();
    } catch (std::exception& e) {
        return e.what();
    }
}

/*
std::string displayTable(const Table &table, const std::vector<std::string> &colNames, const std::vector<std::string> &rowNames) {
    std::stringstream ss;

    const auto& data = table.rows;

    // Rest of the function remains the same...

    if (data.size() != rowNames.size() || (data.size() > 0 && data[0].size() != colNames.size())) {
        return "Invalid data or column/row names!";
    }

    std::vector<size_t> colWidths(colNames.size(), 0);
    for (size_t j = 0; j < colNames.size(); ++j) {
        colWidths[j] = colNames[j].size();
        for (size_t i = 0; i < data.size(); ++i) {
            colWidths[j] = std::max(colWidths[j], variantToStr(data[i][j]).size());
        }
    }

    ss << std::setw(8) << " ";
    for (size_t j = 0; j < colNames.size(); ++j) {
        ss << std::setw(colWidths[j] + 2) << colNames[j];
    }
    ss << '\n';

    for (size_t i = 0; i < data.size(); ++i) {
        ss << std::setw(8) << rowNames[i];
        for (size_t j = 0; j < data[i].size(); ++j) {
            ss << std::setw(colWidths[j] + 2) << variantToStr(data[i][j]);
        }
        ss << '\n';
    }

    return ss.str();
}

std::string toCSV(const Table &table, const std::vector<std::string> &colNames, const std::vector<std::string> &rowNames) {
    std::stringstream ss;

    // Check if data dimensions match provided column and row names
    if (table.rows.size() != rowNames.size() || (table.rows.size() > 0 && table.rows[0].size() != colNames.size())) {
        return "Invalid data or column/row names!";
    }

    // Start with column names
    ss << ",";  // The first column will be for row names
    for (size_t j = 0; j < colNames.size(); ++j) {
        ss << colNames[j];
        if (j < colNames.size() - 1) ss << ",";
    }
    ss << '\n';

    // Now, the data
    for (size_t i = 0; i < table.rows.size(); ++i) {
        ss << rowNames[i] << ",";  // Row name first
        for (size_t j = 0; j < table.rows[i].size(); ++j) {
            ss << variantToStr(table.rows[i][j]);
            if (j < table.rows[i].size() - 1) ss << ",";
        }
        ss << '\n';
    }

    return ss.str();
}

std::string variantToStr(const DataVariant &var) {
    if (std::holds_alternative<float>(var)) {
        return std::to_string(std::get<float>(var));
    } else {
        return std::get<std::string>(var);
    }
}
*/

std::vector<std::string> getAllFiles(std::string folderName)
{
    std::vector<std::string> filenames;
    QDirIterator it(QString::fromStdString(folderName), QDir::Files, QDirIterator::Subdirectories);
    while (it.hasNext()) {
        QString dir = it.next();
        filenames.push_back(dir.toStdString());
    }
    return filenames;
}

// Function to perform FFT
std::vector<std::complex<float>> fft(const std::vector<std::complex<float>>& x, bool inverse) {
    const size_t N = x.size();
    if (N <= 1) return x;

    std::vector<std::complex<float>> result(N);

    // Bit-reversal permutation (optional but enhances performance)
    std::vector<size_t> permutation(N);
    size_t logN = static_cast<size_t>(std::log2(N));
    for (size_t i = 0; i < N; ++i) {
        size_t j = 0;
        for (size_t bit = 0; bit < logN; ++bit) {
            if (i & (1 << bit)) {
                j |= (1 << (logN - 1 - bit));
            }
        }
        permutation[i] = j;
    }

    // Perform FFT
    for (size_t i = 0; i < N; ++i) {
        result[permutation[i]] = x[i];
    }

    for (size_t s = 1; s <= logN; ++s) {
        size_t m = 1 << s;
        std::complex<float> wm = std::polar(1.0f, (inverse ? -2.0f : 2.0f) * float(M_PI) / float(m));
        for (size_t k = 0; k < N; k += m) {
            std::complex<float> w = 1.0f;
            for (size_t j = 0; j < m / 2; ++j) {
                std::complex<float> t = w * result[k + j + m / 2];
                std::complex<float> u = result[k + j];
                result[k + j] = u + t;
                result[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }

    if (inverse) {
        for (size_t i = 0; i < N; ++i) {
            result[i] /= N; // Scaling for the inverse FFT
        }
    }

    return result;
}

// Function to perform the inverse FFT on a given FFT result
std::vector<std::complex<float>> inverseFFT(const std::vector<std::complex<float>>& fft_result) {
    return fft(fft_result, true);
    /*
    size_t size = fft_result.size();
    std::vector<std::complex<float>> conjugate(size);

    // Take the conjugate of the FFT result
    for (size_t i = 0; i < size; ++i) {
        conjugate[i] = std::conj(fft_result[i]);
    }

    // Perform another FFT (or IFFT) on the conjugate result
    std::vector<std::complex<float>> inverse_result = fft(conjugate); // Assuming fft() is your IFFT function

    // Normalize the result by dividing by the size
    for (size_t i = 0; i < size; ++i) {
        inverse_result[i] /= static_cast<float>(size);
    }

    return inverse_result;
    */
}

bool isPowerOf2(int n)
{
    return n > 0 &&  (n & (n-1)) == 0;
}

int runCommand(std::string command)
{
    command = "/bin/bash -c 'source ~/.bashrc && " + command + "'";
    return std::system(command.c_str());
}

Vector3 colorPalette(float t, const Vector3 &startColor, const Vector3 &endColor)
{
    return Vector3::slerp(t, startColor, endColor);
}

void displayProcessTime(std::string textToDisplay, std::function<void ()> func, bool print)
{
    if (print) std::cout << textToDisplay << std::flush;
    float time = timeIt(func);
    if (print) std::cout << " " << showTime(time) << std::endl;
}

Vector3 colorPalette(float t, const std::vector<Vector3> &colors)
{
    return BSpline(colors).getPoint(t);
}
