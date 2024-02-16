#include "OpenFoamParser.h"
#include "TerrainGen/VoxelGrid.h"
#include "Graphics/CubeMesh.h"

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

float parseOpenFoamScale(const std::string& content) {
    std::istringstream iss(content);
    std::string line;

    while (std::getline(iss, line)) {
        if (line.find("convertToMeters") != std::string::npos) {
            break;
        }
    }

    line = line.substr(line.find("convertToMeters") + std::string("convertToMeters").size());
    line = line.substr(0, line.find(";"));
    float scale = std::atof(line.c_str());
    return scale;
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
            points.push_back(Vector3(x, y, z) * 20.f);
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
    return AABBox(points);
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

GridV3 transformToGrid(std::map<int, Cell>& cells, std::vector<Face>& faces, std::vector<Vector3>& points, const AABBox& meshBoundaries, const Vector3& scaling) {
    AABBox limits(points);
    Vector3 boundaryDimensions = meshBoundaries.dimensions();
    boundaryDimensions = boundaryDimensions.xzy();
//    std::cout << "All points in range: " << limits << "\nUseful in range " << meshBoundaries << std::endl;
    auto transfo = [&](const Vector3& point) {
        Vector3 newPos = (point - meshBoundaries.mini);
//        newPos.y = meshBoundaries.dimensions().y - newPos.y;
//        Vector3 newPos = ((point - limits.mini) / (limits.maxi - limits.mini)) * meshBoundaries.dimensions();
        return Vector3(newPos.x, boundaryDimensions.y - newPos.z, newPos.y) / scaling;
    };
    GridV3 velocities(boundaryDimensions / scaling); //(dimX, dimY, dimZ);
//    std::cout << "Velocity is of size " << velocities.getDimensions() << std::endl;

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
    bool regularGrid = true;
    int highestIteration = -1;
    std::vector<std::string> filenames;
    QDirIterator it(QString::fromStdString(foldername), QDir::Dirs);
    while (it.hasNext()) {
        QString dir = it.next();
        int iteration = std::atoi(getFilename(dir.toStdString()).c_str());
        if (iteration > highestIteration)
            highestIteration = iteration;
    }
    GridV3 velocities;

    displayProcessTime("Reading velocities from " + foldername + "/" + std::to_string(highestIteration) + "... ", [&]() {
        std::ifstream blockMeshFile(foldername + "/system/blockMeshDict");
        std::stringstream blockMeshData;
        blockMeshData << blockMeshFile.rdbuf();
        blockMeshFile.close();
        float scaleFactor = parseOpenFoamScale(blockMeshData.str());
        Vector3 scaling = Vector3(scaleFactor, scaleFactor, scaleFactor);
        blockMeshData.clear();

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


        velocities = transformToGrid(cells, faces, points, boundaries, scaling);
        float offsetX = .25f;
        float offsetY = .25f;
        AABBox meshDimensions(Vector3(offsetX, offsetY, .0f) * velocities.getDimensions(), Vector3(1.f - offsetX, 1.f - offsetY, 1.f) * velocities.getDimensions());
    //    AABBox meshDimensions = AABBox(Vector3(), velocities.getDimensions());
        velocities = velocities.subset(meshDimensions.min(), meshDimensions.max());
    });
    return velocities;
}





std::string OpenFoamParser::createSimulationFile(std::string foldername, const GridF& _boundaries)
{
    float scaleFactor = 1000.f;
    Mesh m = Mesh::applyMarchingCubes(_boundaries);
    std::string file_name = foldername + "/constant/triSurface/cylinder.stl";
    Vector3 targetSize = Vector3(1.f, .5f, 3.f);
    Vector3 rescale = (Vector3(1.f, 1.f, 1.f) * scaleFactor) / Vector3(_boundaries.sizeX, _boundaries.sizeY, _boundaries.sizeZ - 1); // Take into account that we will reduce Z by 1 unit
    m.translate(-(_boundaries.getDimensions().xy() * .5f + Vector3(0, 0, 1)));
    m.scale(rescale);
    m.rotate(deg2rad(90), 0, 0);

    displayProcessTime("Saving as STL at " + file_name + "... ", [&]() {
        std::ofstream off(file_name);
        off << m.toSTL();
        off.close();
    });

    displayProcessTime("Update of blockMeshDict and snappyHexMeshDict...", [&]() {
        std::string line;
        std::vector<std::string> fileContent;

        std::ifstream blockFileIn(foldername + "/system/blockMeshDict");
        while (std::getline(blockFileIn, line)) {
            size_t patternFoundAt = line.find("convertToMeters");
            if (patternFoundAt != std::string::npos)
                line = "convertToMeters " + replaceInString(std::to_string(scaleFactor), ",", ".") + ";";
            fileContent.push_back(line);
        }
        blockFileIn.close();
        std::ofstream blockFileOut(foldername + "/system/blockMeshDict");
        blockFileOut << join(fileContent, "\n");
        blockFileOut.close();
        fileContent.clear();

        std::ifstream snappyHexMeshIn(foldername + "/system/snappyHexMeshDict");
        while (std::getline(snappyHexMeshIn, line)) {
            size_t patternFoundAt = line.find("locationInMesh");
            if (patternFoundAt != std::string::npos)
                line = "\tlocationInMesh (0 " + replaceInString(std::to_string(.9999f * scaleFactor), ",", ".") + " 0);";
            fileContent.push_back(line);
        }
        snappyHexMeshIn.close();
        std::ofstream snappyHexMeshOut(foldername + "/system/snappyHexMeshDict");
        snappyHexMeshOut << join(fileContent, "\n");
        snappyHexMeshOut.close();
        fileContent.clear();
    });

    std::vector<std::string> foldersToRemove;
    QDirIterator it(QString::fromStdString(foldername), QDir::Dirs);
    while (it.hasNext()) {
        QString dir = it.next();
        int iteration = std::atoi(getFilename(dir.toStdString()).c_str());
        if (iteration > 0) {
            std::string folderWithGuards = "\"" + foldername + "/" + getFilename(dir.toStdString()) + "\"";
            foldersToRemove.push_back(folderWithGuards);
        }
    }

    std::vector<std::string> commands = {
        "blockMesh -case \"" + foldername + "\"",
        "surfaceFeatures -case \"" + foldername + "\"",
        "snappyHexMesh -case \"" + foldername + "\"",
        "rm -rf \"" + foldername + "/constant/polyMesh\" && cp -rf \"" + foldername + "/0.5/polyMesh/\" \"" + foldername + "/constant/\" && rm -rf \"" + foldername + "/0.5/\"",
        "rm -rf " + join(foldersToRemove, " "),
        "simpleFoam -case \"" + foldername + "\"",
    };

//    std::string command;
    int result;

    for (auto command : commands) {
        displayProcessTime("Running " + command + "...", [&]() {
            command += " > /dev/null";
            result = std::system(command.c_str());
        });
        if (result != 0) {
            std::cerr << "Oups, the command `" << command << "` didn't finished as expected... Maybe OpenFOAM is missing?" << std::endl;
            break;
        }
    }

    return "";
}
