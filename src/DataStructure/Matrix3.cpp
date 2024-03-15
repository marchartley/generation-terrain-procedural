#include "DataStructure/Matrix3.h"
//#include "Utils/stb_image.h"
#include "Utils/Skeletonize.h"
#include <queue>

template<>
Matrix3<Vector3> Matrix3<Vector3>::curl(float radius) const {
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
//    float radius = 1;
    iterateParallel([&] (int x, int y, int z) {
        // Central finite difference
        Vector3 dX = at(x + radius, y, z) - at(x - radius, y, z);
        Vector3 dY = at(x, y + radius, z) - at(x, y - radius, z);
        Vector3 dZ = at(x, y, z + radius) - at(x, y, z - radius);

        Vector3 curl(dY.z - dZ.y, dZ.x - dX.z, dX.y - dY.x);
        returningGrid(x, y, z) = curl / (radius * 2.f);
//        Vector3 dF = Vector3(at(x + radius, y, z) - at(x - radius, y, z), at(x, y + radius, z) - at(x, y - radius, z), at(x, y, z + radius) - at(x, y, z - radius));
//        returningGrid(x, y, z) = Vector3(dF.z - dF.y, dF.x - dF.z, dF.y - dF.x) / (2 * radius);
//        const Vector3& vec = this->at(x, y, z);
//        returningGrid.at(x, y, z) = Vector3(vec.z - vec.y, vec.x - vec.z, vec.y - vec.x);
    });
    /*#pragma omp parallel for collapse(3)
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                Vector3& vec = this->at(x, y, z);
                returningGrid.at(x, y, z) = Vector3(vec.z - vec.y, vec.x - vec.z, vec.y - vec.x);
            }
        }
    }*/
    return returningGrid;
}
template<>
Matrix3<Vector3> Matrix3<Vector3>::rot() const {
    return this->curl();
}

template<>
Matrix3<int> Matrix3<int>::skeletonize() const
{
    Matrix3<int> self = ((Matrix3<float>)*this).binarize(0.5f);
    Matrix3<int> initial = *this;
    skeleton_tracer_t* skel = new skeleton_tracer_t();
    skel->W = self.sizeX; // width of image
    skel->H = self.sizeY; // height of image

    // allocate the input image
    unsigned char* data = (unsigned char*)malloc(sizeof(unsigned char)*skel->W*skel->H); //new uchar(self.sizeX * self.sizeY);
    self.iterateParallel([&] (size_t i) {
        data[i] = (unsigned char)(self[i]);
    });
    /*for (size_t i = 0; i < self.size(); i++)
        data[i] = (unsigned char)(self[i]);*/
    skel->im = data;

    skel->thinning_zs(); // perform raster thinning

    // run the algorithm
    skeleton_tracer_t::polyline_t* p = (skeleton_tracer_t::polyline_t*)skel->trace_skeleton(0, 0, skel->W, skel->H, 0);

    self.reset(); // Back to all 0's
    // print out points in every polyline
    skeleton_tracer_t::polyline_t* it = p; //iterator
    while(it){
      skeleton_tracer_t::point_t* jt = it->head;
      Vector3 prevPos(false);
      while(jt){
          Vector3 newPos(jt->x, jt->y);
          // Clearly not the best way to do this!
          if (prevPos.isValid()/* && initial.at(newPos) != 0*/) {
              Vector3 move = newPos - prevPos;
              float length = move.length();
              for (int i = 0; i < length; i++)
                  self.at(prevPos + move * clamp(i / length, 0.f, 1.f)) = 1;
          }
          prevPos = newPos;
        jt = jt->next;
      }
      it = it->next;
    }
    free(skel->im);
    skel->destroy_polylines(p);
    skel->destroy_rects();
    delete skel;
    return self;
}

template<>
std::vector<BSpline> Matrix3<int>::skeletonizeToBSplines() const
{
    Matrix3<int> initial = ((Matrix3<float>)*this).binarize(0.5f);
    skeleton_tracer_t* skel = new skeleton_tracer_t();
    skel->W = this->sizeX; // width of image
    skel->H = this->sizeY; // height of image

    // allocate the input image
    unsigned char* data = (unsigned char*)malloc(sizeof(unsigned char)*skel->W*skel->H); //new uchar(self.sizeX * self.sizeY);
    iterateParallel([&] (size_t i) {
        data[i] = (unsigned char)(initial(i));
    });
    /*for (size_t i = 0; i < self.size(); i++)
        data[i] = (unsigned char)(self[i]);*/
    skel->im = data;

    skel->thinning_zs(); // perform raster thinning

    // run the algorithm
    skeleton_tracer_t::polyline_t* p = (skeleton_tracer_t::polyline_t*)skel->trace_skeleton(0, 0, skel->W, skel->H, 0);


    std::vector<BSpline> splines;
    skeleton_tracer_t::polyline_t* it = p; //iterator
    while(it){
      skeleton_tracer_t::point_t* jt = it->head;
      BSpline spline;
      while(jt){
          ;
          spline.points.push_back(Vector3(jt->x, jt->y));
          jt = jt->next;
      }
      it = it->next;
      splines.push_back(spline);
    }
    free(skel->im);
    skel->destroy_polylines(p);
    skel->destroy_rects();
    delete skel;
//    return splines;

    // Try at best to merge the curves
    float limitDistance = 20.f;
    float sqrLim = limitDistance * limitDistance;
    bool atLeastOneMerging = true;
    while (atLeastOneMerging) {
        atLeastOneMerging = false;
        std::vector<BSpline> merged;
        for (int i = int(splines.size()) - 1; i >= 0; i--) {
            auto& spline = splines[i];
            bool isMerged = false;
            for (auto& merge : merged) {
                // Try back to front, front to back, back to back and front to front
                auto frontSpline = spline.points.front();
                auto backSpline = spline.points.back();
                auto frontMerge = merge.points.front();
                auto backMerge = merge.points.back();
                // back to front:
                if ((frontMerge - backSpline).norm2() < sqrLim) {
                    merge.points.insert(merge.points.begin(), spline.points.begin(), spline.points.end());
                    isMerged = true;
                    break;
                }
                // front to back:
                else if ((backMerge - frontSpline).norm2() < sqrLim) {
                    merge.points.insert(merge.points.end(), spline.points.begin(), spline.points.end());
                    isMerged = true;
                    break;
                }
                // back to back:
                else if ((backMerge - backSpline).norm2() < sqrLim) {
                    std::reverse(spline.points.begin(), spline.points.end());
                    merge.points.insert(merge.points.end(), spline.points.begin(), spline.points.end());
                    isMerged = true;
                    break;
                }
                // front to front:
                else if ((frontMerge - frontSpline).norm2() < sqrLim) {
                    std::reverse(spline.points.begin(), spline.points.end());
                    merge.points.insert(merge.points.begin(), spline.points.begin(), spline.points.end());
                    isMerged = true;
                    break;
                }
            }
            if (!isMerged)
                merged.push_back(spline);
            else
                atLeastOneMerging = true;
        }
        splines = merged;
    }
    return splines;
}

template<>
Matrix3<int> Matrix3<int>::computeConnectedComponents(bool use4Connect) const
{
    int currentLabel = 1;
    GridI labelMap = (getDimensions());  // Initialize labelMap with zeros
    std::vector<std::vector<int>> equivalences;

//    iterateParallel([&] (size_t i) {
//        labelMap[i] = (this->at(i) == 0 ? 0 : i + 1);
//    });

    auto findForegroundNeighbors = [&](int x, int y, int z) {
        std::vector<int> neighbors;
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (use4Connect && !(dx == 0 || dy == 0 || dz == 0)) continue;
                    int _x = x + dx, _y = y + dy, _z = z + dz;
                    if (checkCoord(_x, _y, _z) && this->at(_x, _y, _z) != 0 && labelMap.at(_x, _y, _z) != 0)
                        neighbors.push_back(labelMap(_x, _y, _z));
                }
            }
        }
        return neighbors;
    };

    auto recordEquivalence = [&](const std::vector<int>& neighbors, int minNeighborLabel) {
        std::vector<int> newEquivalenceGroup;
        newEquivalenceGroup.push_back(minNeighborLabel);

        for (int neighborLabel : neighbors) {
            if (neighborLabel != minNeighborLabel) {
                newEquivalenceGroup.push_back(neighborLabel);
            }
        }

        // Check if any of the labels in the new group already appear in existing groups.
        for (auto& existingGroup : equivalences) {
            bool found = false;
            for (int label : newEquivalenceGroup) {
                if (std::find(existingGroup.begin(), existingGroup.end(), label) != existingGroup.end()) {
                    found = true;
                    break;
                }
            }

            // If a label from the new group appears in an existing group, merge the groups.
            if (found) {
                existingGroup.insert(existingGroup.end(), newEquivalenceGroup.begin(), newEquivalenceGroup.end());
                // Remove duplicates and sort.
                std::sort(existingGroup.begin(), existingGroup.end());
                auto last = std::unique(existingGroup.begin(), existingGroup.end());
                existingGroup.erase(last, existingGroup.end());
                return;
            }
        }

        // If the new group doesn't overlap with any existing groups, add it as a new group.
        equivalences.push_back(newEquivalenceGroup);
    };

    auto findRootLabel = [&](int label) {
        // Iterate through each equivalence group
        for (const auto& group : equivalences) {
            // Check if the label is in the current group
            if (std::find(group.begin(), group.end(), label) != group.end()) {
                // Return the smallest label in the group
                return *std::min_element(group.begin(), group.end());
            }
        }
        // If the label does not appear in any equivalence group, it's its own root
        return label;
    };

    // First pass
    iterate([&] (int i, int j, int k) {
        if (this->at(i, j, k) != 0) {  // Assuming foreground is represented by non-zero values
            std::vector<int> neighbors = findForegroundNeighbors(i, j, k);
            if (neighbors.empty()) {
                labelMap(i, j, k) = currentLabel++;
            } else {
                int minNeighborLabel = *std::min_element(neighbors.begin(), neighbors.end());
                labelMap(i, j, k) = minNeighborLabel;
                recordEquivalence(neighbors, minNeighborLabel);
            }
        }
    });

    // Reduce labels values
    std::map<int, int> rootToNewLabel;
    int newLabel = 1;
    for (const auto& group : equivalences) {
        int rootLabel = *std::min_element(group.begin(), group.end());
        if (rootToNewLabel.find(rootLabel) == rootToNewLabel.end()) {
            rootToNewLabel[rootLabel] = newLabel++;
        }
    }

    // Second pass
    labelMap.iterateParallel([&] (size_t i) {
        if (labelMap[i] != 0) {
            int rootLabel = findRootLabel(labelMap[i]);
            labelMap[i] = rootToNewLabel[rootLabel];
        }
    });

    return labelMap;
}

template<>
Matrix3<int> Matrix3<int>::findContour(bool use2D) const
{
//    return this->dilate(true) - *this;
    return *this - this->erode(use2D);
}

template<>
std::vector<ShapeCurve> Matrix3<int>::findContoursAsCurves() const
{
    std::vector<ShapeCurve> curves;
    auto grid = this->resize(sizeX * 2.f, sizeY * 2.f, RESIZE_MODE::MAX_VAL).findContour(true);
//    std::cout << "Grid : " << grid << "\n" << this->resize(sizeX * .5f, sizeY * .5f, RESIZE_MODE::MAX_VAL).findContour(true).displayAsPlot() << std::endl;
    grid.raiseErrorOnBadCoord = false;

    while (grid.max() != 0) {
        std::vector<Vector3> contour;
        // Find the starting point (the first 1 encountered in the grid)
        Vector3 start(false);

        grid.iterate([&] (size_t i) {
            if (!start.isValid() && grid[i] == 1) {
                start = grid.getCoordAsVector3(i);
                return;
            }
        });

        // If no point was found, return the empty contour
        if (!start.isValid()) return curves;

        // Start from the initial point and do a simple DFS to find the contour
        std::vector<Vector3> q;
        q.push_back(start);
        grid(start) = 0;  // Mark as visited

        while (!q.empty()) {
            Vector3 current = q.back();
            q.pop_back();
            // Add to contour
            contour.push_back(current);

            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    if (!(dx == 0 || dy == 0)) continue;
                    Vector3 nextPos(current.x + dx, current.y + dy);
                    if (grid(nextPos)) {
                        q.push_back(nextPos);
                        grid(nextPos) = 0;
                    }
                }
            }
        }

        for (auto& p : contour)
            p *= .5f;
        curves.push_back(BSpline(contour)/*.resamplePoints().simplifyByRamerDouglasPeucker(5.f).resamplePoints()*/);
    }
    return curves;
}

template <>
Matrix3<float> Matrix3<Vector3>::divergence() const
{
    auto self = *this;
    self.raiseErrorOnBadCoord = false;
    Matrix3<float> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    iterateParallel([&] (int x, int y, int z) {
        returningGrid.at(x, y, z) = ((self.at(x + 1, y, z) - self.at(x - 1, y, z)).x +
                                     (self.at(x, y + 1, z) - self.at(x, y - 1, z)).y +
                                     (self.at(x, y, z + 1) - self.at(x, y, z - 1)).z) * .5f;
    });
    return returningGrid;
}

template<>
Vector3 Matrix3<Vector3>::gradient(const Vector3& position) const
{
    auto self = *this;
    self.raiseErrorOnBadCoord = false;
    Vector3 flooredPos = position.floor();
    Vector3 offset = position - flooredPos;
    return Vector3(
                self.at(flooredPos + Vector3(1, 0, 0)).x * (1 - offset.x) + self.at(flooredPos).x * offset.x,
                self.at(flooredPos + Vector3(0, 1, 0)).y * (1 - offset.y) + self.at(flooredPos).y * offset.y,
                self.at(flooredPos + Vector3(0, 0, 1)).z * (1 - offset.z) + self.at(flooredPos).z * offset.z
                );
}

template<>
Vector3 Matrix3<Vector3>::gradient(float posX, float posY, float posZ) const
{
    return gradient(Vector3(posX, posY, posZ));
}

template<>
Matrix3<Vector3> Matrix3<Vector3>::gradient() const
{
    auto self = *this;
    self.raiseErrorOnBadCoord = false;
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    iterateParallel([&] (int x, int y, int z) {
        returningGrid.at(x, y, z) = Vector3((self.at(x + 1, y, z) - self.at(x - 1, y, z)).x * .5f,
                                            (self.at(x, y + 1, z) - self.at(x, y - 1, z)).y * .5f,
                                            (self.at(x, y, z + 1) - self.at(x, y, z - 1)).z * .5f);
    });
    /*#pragma omp parallel for collapse(3)
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
//                returningGrid.at(x, y, z) = gradient(x, y, z);
//                continue;
                // Need to change the divergence function...
//                returningGrid.at(x, y, z) = this->at(x, y, z).divergence();
                returningGrid.at(x, y, z) = Vector3((this->at(x + 1, y, z) - this->at(x - 1, y, z)).x * .5f,
                                                    (this->at(x, y + 1, z) - this->at(x, y - 1, z)).y * .5f,
                                                    (this->at(x, y, z + 1) - this->at(x, y, z - 1)).z * .5f);
            }
        }
    }*/
//    this->raiseErrorOnBadCoord = true;
    return returningGrid;
}


template<>
Matrix3<Vector3> Matrix3<Vector3>::random(size_t sizeX, size_t sizeY, size_t sizeZ)
{
    Matrix3<Vector3> mat(sizeX, sizeY, sizeZ);
    mat.iterateParallel([&] (size_t i) {
        mat[i] = Vector3::random();
    });
    /*for (Vector3& val : mat)
        val = Vector3::random();*/
    return mat;
}

template<>
Matrix3<Vector3> Matrix3<Vector3>::fromImageRGB(std::string filename)
{
    if (!checkPathExists(filename)) {
        throw std::runtime_error("Error: Impossible to load '" + filename + "', file not found");
    }
    int imgW, imgH, nbChannels;
    unsigned char *c_data = stbi_load(filename.c_str(), &imgW, &imgH, &nbChannels, STBI_rgb); // Load image, force 3 channel
    if (c_data == NULL)
    {

        throw std::runtime_error("Error : Impossible to load RGB image at " + filename + "\n" +
                                "Either file is not found, or type is incorrect. Available file types are : \n" +
                                "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC");
    }
    float *data = new float[imgW * imgH * 3];
    for (int i = 0; i < imgW * imgH * 3; i++)
        data[i] = c_data[i];
    stbi_image_free(c_data);

    Matrix3<Vector3> map(imgW, imgH);
    map.iterateParallel([&] (int x, int y, int _z) {
        int index = x + y * imgW;
        map(x, y) = Vector3(data[3 * index + 0],
                data[3 * index + 1],
                data[3 * index + 1]) / 255.f;
    });
    /*for (int x = 0; x < imgW; x++) {
        for (int y = 0; y < imgH; y++) {
            int index = x + y * imgW;
            map.at(x, y) = Vector3(data[3 * index + 0],
                    data[3 * index + 1],
                    data[3 * index + 1]) / 255.f;
        }
    }*/
    if (data != nullptr)
        delete[] data;//stbi_image_free(data);

    return map;
}

template<>
Matrix3<float> Matrix3<float>::fromImageBW(std::string filename)
{
    if (!checkPathExists(filename)) {
        throw std::runtime_error("Error: Impossible to load '" + filename + "', file not found");
    }
    int imgW, imgH, nbChannels;
    unsigned char *c_data = stbi_load(filename.c_str(), &imgW, &imgH, &nbChannels, STBI_grey); // Load image, force 1 channel
    if (c_data == NULL)
    {

        throw std::runtime_error("Error : Impossible to load BW image at " + filename + "\n" +
                                "Either file is not found, or type is incorrect. Available file types are : \n" +
                                "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC");
    }
    float *data = new float[imgW * imgH];
    for (int i = 0; i < imgW * imgH; i++)
        data[i] = c_data[i];
    stbi_image_free(c_data);

    Matrix3<float> map(imgW, imgH);
    map.iterateParallel([&] (int x, int y, int _z) {
        int index = x + y * imgW;
        map(x, y) = data[index] / 255.f;
    });
    /*for (int x = 0; x < imgW; x++) {
        for (int y = 0; y < imgH; y++) {
            int index = x + y * imgW;
            map.at(x, y) = data[index] / 255.f;
        }
    }*/
    if (data != nullptr)
        delete[] data;//stbi_image_free(data);

    return map;
}
//template<class T>
Matrix3<float> operator-(const float a, Matrix3<float> b) {
    Matrix3<float> res = b;
    res.iterateParallel([&] (size_t i) {
        res[i] = a - res[i];
    });
    /*for (size_t i = 0; i < res.size(); i++)
        res[i] = a - res[i];*/
    return res;
}
//template<class T>
Matrix3<float> operator+(const float a, Matrix3<float> b) {
    return b + a;
    /*Matrix3<float> res = b;
    for (size_t i = 0; i < res.size(); i++)
        res[i] = a + res[i];
    return res;*/
}
