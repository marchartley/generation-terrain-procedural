#include "OpenFoamParser.h"

OpenFoamParser::OpenFoamParser()
{

}





std::vector<Vector3> parseOpenFoamUFile(const std::string& content) {
    std::vector<Vector3> velocities;
    std::istringstream iss(content);
    std::string line;

    // Search for the start of the internalField data
    while (std::getline(iss, line)) {
        if (line.find("internalField") != std::string::npos) {
            break;
        }
    }

    // Skip next lines until we find the starting parenthesis
    while (std::getline(iss, line)) {
        if (line.find("(") != std::string::npos) {
            break;
        }
    }

    // Read the vectors
    while (std::getline(iss, line)) {
        std::istringstream lineStream(line);
        char openPar, closePar;
        double x, y, z;

        if (lineStream >> openPar >> x >> y >> z >> closePar) {
            velocities.emplace_back(x, y, z);
        } else {
            // Assuming the vectors end when we can't read a proper format
            break;
        }
    }

    return velocities;
}

std::vector<Vector3> parseOpenFoamPointsFile(const std::string& content) {
    std::vector<Vector3> points;
    std::istringstream iss(content);
    std::string line;

    // Skip next lines until we find the starting parenthesis
    while (std::getline(iss, line)) {
        if (line.find("(") != std::string::npos) {
            break;
        }
    }

    // Read the vectors
    while (std::getline(iss, line)) {
        std::istringstream lineStream(line);
        char openPar, closePar;
        double x, y, z;

        if (lineStream >> openPar >> x >> y >> z >> closePar) {
            points.emplace_back(x, y, z);
        } else {
            // Assuming the vectors end when we can't read a proper format
            break;
        }
    }

    return points;
}

struct Face {
    std::vector<int> points;

    Face(int a, int b, int c, int d)
        : points({a, b, c, d})
    {
    }
    Face(int a, int b, int c)
        : points({a, b, c})
    {
    }
};

std::vector<Face> parseOpenFoamFacesFile(const std::string& content) {
    std::vector<Face> faces;
    std::istringstream iss(content);
    std::string line;

    // Skip next lines until we find the starting parenthesis
    while (std::getline(iss, line)) {
        if (line.find("(") != std::string::npos) {
            break;
        }
    }

    // Read the vectors
    while (std::getline(iss, line)) {
        std::istringstream lineStream(line);
        char nbPoints, openPar, closePar;
        int v1, v2, v3, v4;

        if (lineStream >> nbPoints >> openPar >> v1 >> v2 >> v3 >> v4 >> closePar) {
            faces.emplace_back(Face(v1, v2, v3, v4));
        } else {
            if (nbPoints == '3') {
                faces.emplace_back(Face(v1, v2, v3));
            } else {
                break;
            }
        }
    }

    return faces;
}

struct Cell {
    std::vector<int> faces;
    Vector3 velocity;

    Vector3 getMidpoint(std::vector<Face>& _faces, std::vector<Vector3>& points) {
        Vector3 sum;
        int divisor = 0;
        for (auto& fId : faces) {
            auto& f = _faces[fId];
            for (auto& p : f.points) {
                sum += points[p];
                divisor ++;
            }
        }
        return sum / float(divisor);
    }
};

std::map<int, Cell> parseOpenFoamCellsFile(const std::string& content) {
    std::map<int, Cell> cells;
    std::istringstream iss(content);
    std::string line;

    // Skip next lines until we find the starting parenthesis
    while (std::getline(iss, line)) {
        if (line.find("(") != std::string::npos) {
            break;
        }
    }

    // Read the vectors
    int faceNumber = 0;
    while (std::getline(iss, line)) {
        std::istringstream lineStream(line);
        int cellNumber;
        if (lineStream >> cellNumber) {
            if (cells.count(cellNumber) == 0)
                cells[cellNumber] = Cell();
            cells[cellNumber].faces.push_back(faceNumber);
            faceNumber++;
        } else {
            // Assuming the vectors end when we can't read a proper format
            break;
        }
    }
    return cells;
}

AABBox parseOpenFoamBoundariesFile(const std::string& content, std::vector<Face>& faces, std::vector<Vector3>& points) {
    std::istringstream iss(content);
    std::string line;

    int nFaces, startFace;

    // Skip next lines until we find the starting parenthesis
    while (std::getline(iss, line)) {
        if (line.find("(") != std::string::npos) {
            break;
        }
    }
    // Skip next lines until we find the starting parenthesis
    while (std::getline(iss, line)) {
        if (line.find("cylinder") != std::string::npos) {
            std::getline(iss, line); // Ignore the next line (just a "}" )
            break;
        }
    }

    // Read the vectors
    while (std::getline(iss, line)) {
        std::istringstream lineStream(line);
        std::string var_name;
        std::string value;

        if (lineStream >> var_name >> value) {
            if (var_name == "nFaces") nFaces = std::atoi(value.c_str());
            else if (var_name == "startFace") startFace = std::atoi(value.c_str());
        } else {
            break;
        }
    }

    AABBox limits = AABBox({points[faces[startFace].points.front()]});
//    Vector3 mini = Vector3::min(), maxi = Vector3::max();

    for (int iFaceOffset = 0; iFaceOffset < nFaces; iFaceOffset++) {
        for (const auto& point : faces[startFace + iFaceOffset].points) {
            limits.expand(points[point]);
        }
    }
    return limits;
}
void linkUandCells(std::vector<Vector3>& U, std::map<int, Cell>& cells) {
    for (int i = 0; i < U.size(); i++) {
        if (cells.count(i) != 0)
            cells[i].velocity = U[i];
    }
}

GridV3 transformToGrid(std::map<int, Cell>& cells, std::vector<Face>& faces, std::vector<Vector3>& points, const AABBox& meshBoundaries) {

    AABBox limits(points);
    Vector3 boundaryDimensions = meshBoundaries.dimensions();
    boundaryDimensions = boundaryDimensions.xzy();
    std::cout << "All points in range: " << limits << "\nUseful in range " << meshBoundaries << std::endl;
    auto transfo = [&](const Vector3& point) {
        Vector3 newPos = (point - meshBoundaries.mini);
//        newPos.y = meshBoundaries.dimensions().y - newPos.y;
//        Vector3 newPos = ((point - limits.mini) / (limits.maxi - limits.mini)) * meshBoundaries.dimensions();
        return Vector3(newPos.x, boundaryDimensions.y - newPos.z, newPos.y);
    };
    GridV3 velocities(boundaryDimensions); //(dimX, dimY, dimZ);
    std::cout << "Velocity is of size " << velocities.getDimensions() << std::endl;

    for (auto& [idx, cell] : cells) {
        Vector3 midPoint = cell.getMidpoint(faces, points);
        Vector3 gridPos = transfo(midPoint);
        Vector3 cellVelocity(cell.velocity.x, -cell.velocity.z, cell.velocity.y);
        velocities.addValueAt(cellVelocity, gridPos);
//        velocities.at(gridPos) += cell.velocity;
    }
    return velocities;
}



GridV3 OpenFoamParser::parseSimulation(std::string foldername)
{
    int highestIteration = -1;
    std::vector<std::string> filenames;
    QDirIterator it(QString::fromStdString(foldername), QDir::Dirs);
    while (it.hasNext()) {
        QString dir = it.next();
        int iteration = std::atoi(getFilename(dir.toStdString()).c_str());
        if (iteration > highestIteration)
            highestIteration = iteration;
    }
    std::ifstream filePoints(foldername + "/constant/polyMesh/points");
    std::stringstream bufferPoints;
    bufferPoints << filePoints.rdbuf();
    filePoints.close();
    auto points = parseOpenFoamPointsFile(bufferPoints.str());
    bufferPoints.clear();

    std::ifstream fileFaces(foldername + "/constant/polyMesh/faces");
    std::stringstream bufferFaces;
    bufferFaces << fileFaces.rdbuf();
    fileFaces.close();
    auto faces = parseOpenFoamFacesFile(bufferFaces.str());
    bufferFaces.clear();

    std::ifstream fileCells(foldername + "/constant/polyMesh/owner");
    std::stringstream bufferCells;
    bufferCells << fileCells.rdbuf();
    fileCells.close();
    auto cells = parseOpenFoamCellsFile(bufferCells.str());
    bufferCells.clear();


    std::ifstream fileU(foldername + "/" + std::to_string(highestIteration) + "/U");
    std::stringstream bufferU;
    bufferU << fileU.rdbuf();
    fileU.close();
    auto U = parseOpenFoamUFile(bufferU.str());
    bufferU.clear();

    linkUandCells(U, cells);


    std::ifstream fileBoundaries(foldername + "/constant/polyMesh/boundary");
    std::stringstream bufferBoundaries;
    bufferBoundaries << fileBoundaries.rdbuf();
    fileBoundaries.close();
    auto boundaries = parseOpenFoamBoundariesFile(bufferBoundaries.str(), faces, points);
    bufferCells.clear();


    auto velocities = transformToGrid(cells, faces, points, boundaries);

    return velocities;
}
