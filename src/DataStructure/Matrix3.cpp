#include "DataStructure/Matrix3.h"
//#include "Utils/stb_image.h"
#include "Utils/Skeletonize.h"

template<>
Matrix3<Vector3> Matrix3<Vector3>::curl() {
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                Vector3& vec = this->at(x, y, z);
                returningGrid.at(x, y, z) = Vector3(vec.z - vec.y, vec.x - vec.z, vec.y - vec.x);
            }
        }
    }
    return returningGrid;
}
template<>
Matrix3<Vector3> Matrix3<Vector3>::rot() { return this->curl(); }

template<>
Matrix3<int> Matrix3<int>::skeletonize()
{
    Matrix3<int> self = ((Matrix3<float>)*this).binarize(0.5f);
    Matrix3<int> initial = *this;
    skeleton_tracer_t* skel = new skeleton_tracer_t();
    skel->W = self.sizeX; // width of image
    skel->H = self.sizeY; // height of image

    // allocate the input image
    unsigned char* data = (unsigned char*)malloc(sizeof(unsigned char)*skel->W*skel->H); //new uchar(self.sizeX * self.sizeY);
    for (size_t i = 0; i < self.size(); i++)
        data[i] = (unsigned char)(self[i]);
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
Matrix3<float> Matrix3<Vector3>::divergence()
{
    this->raiseErrorOnBadCoord = false;
    Matrix3<float> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                // Need to change the divergence function...
//                returningGrid.at(x, y, z) = this->at(x, y, z).divergence();
                returningGrid.at(x, y, z) = ((this->at(x + 1, y, z) - this->at(x - 1, y, z)).x +
                                             (this->at(x, y + 1, z) - this->at(x, y - 1, z)).y +
                                             (this->at(x, y, z + 1) - this->at(x, y, z - 1)).z) * .5f;
            }
        }
    }
    this->raiseErrorOnBadCoord = true;
    return returningGrid;
}

template<>
Vector3 Matrix3<Vector3>::gradient(Vector3 position)
{
    this->raiseErrorOnBadCoord = false;
    Vector3 flooredPos = position.floor();
    Vector3 offset = position - flooredPos;
    return Vector3(
                at(flooredPos + Vector3(1, 0, 0)).x * (1 - offset.x) + at(flooredPos).x * offset.x,
                at(flooredPos + Vector3(0, 1, 0)).y * (1 - offset.y) + at(flooredPos).y * offset.y,
                at(flooredPos + Vector3(0, 0, 1)).z * (1 - offset.z) + at(flooredPos).z * offset.z
                );
}

template<>
Vector3 Matrix3<Vector3>::gradient(float posX, float posY, float posZ)
{
    return gradient(Vector3(posX, posY, posZ));
}

template<>
Matrix3<Vector3> Matrix3<Vector3>::gradient()
{
    this->raiseErrorOnBadCoord = false;
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    #pragma omp parallel for collapse(3)
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
    }
    this->raiseErrorOnBadCoord = true;
    return returningGrid;
}


template<>
Matrix3<Vector3> Matrix3<Vector3>::random(size_t sizeX, size_t sizeY, size_t sizeZ)
{
    Matrix3<Vector3> mat(sizeX, sizeY, sizeZ);
    for (Vector3& val : mat.data)
        val = Vector3::random();
    return mat;
}

template<>
Matrix3<Vector3> Matrix3<Vector3>::fromImageRGB(std::string filename)
{
    int imgW, imgH, nbChannels;
    unsigned char *c_data = stbi_load(filename.c_str(), &imgW, &imgH, &nbChannels, STBI_rgb); // Load image, force 3 channel
    if (c_data == NULL)
    {
        std::cerr << "Error : Impossible to load " << filename << "\n";
        std::cerr << "Either file is not found, or type is incorrect. Available file types are : \n";
        std::cerr << "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC";
        exit (-1);
        return Matrix3<Vector3>();
    }
    float *data = new float[imgW * imgH * 3];
    for (int i = 0; i < imgW * imgH * 3; i++)
        data[i] = c_data[i];
    stbi_image_free(c_data);

    Matrix3<Vector3> map(imgW, imgH);
    for (int x = 0; x < imgW; x++) {
        for (int y = 0; y < imgH; y++) {
            int index = x + y * imgW;
            map.at(x, y) = Vector3(data[3 * index + 0],
                    data[3 * index + 1],
                    data[3 * index + 1]) / 255.f;
        }
    }
    if (data != nullptr)
        delete[] data;//stbi_image_free(data);

    return map;
}

template<>
Matrix3<float> Matrix3<float>::fromImageBW(std::string filename)
{
    int imgW, imgH, nbChannels;
    unsigned char *c_data = stbi_load(filename.c_str(), &imgW, &imgH, &nbChannels, STBI_grey); // Load image, force 1 channel
    if (c_data == NULL)
    {
        std::cerr << "Error : Impossible to load " << filename << "\n";
        std::cerr << "Either file is not found, or type is incorrect. Available file types are : \n";
        std::cerr << "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC";
        exit (-1);
        return Matrix3<float>();
    }
    float *data = new float[imgW * imgH];
    for (int i = 0; i < imgW * imgH; i++)
        data[i] = c_data[i];
    stbi_image_free(c_data);

    Matrix3<float> map(imgW, imgH);
    for (int x = 0; x < imgW; x++) {
        for (int y = 0; y < imgH; y++) {
            int index = x + y * imgW;
            map.at(x, y) = data[index] / 255.f;
        }
    }
    if (data != nullptr)
        delete[] data;//stbi_image_free(data);

    return map;
}
//template<class T>
Matrix3<float> operator-(const float a, Matrix3<float> b) {
    Matrix3<float> res = b;
    for (size_t i = 0; i < res.size(); i++)
        res[i] = a - res[i];
    return res;
}
//template<class T>
Matrix3<float> operator+(const float a, Matrix3<float> b) {
    Matrix3<float> res = b;
    for (size_t i = 0; i < res.size(); i++)
        res[i] = a + res[i];
    return res;
}
