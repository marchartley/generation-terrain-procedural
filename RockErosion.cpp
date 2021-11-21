#include "RockErosion.h"

RockErosion::RockErosion()
{

}
RockErosion::RockErosion(int size, float maxStrength)
    : size(size), maxStrength(maxStrength)
{
    std::cout.precision(2);
    attackMask = new float**[size];
    float maxVal = sqrt(3);
    for (int i = 0; i < size; i++) {
        attackMask[i] = new float*[size];
        for (int j = 0; j < size; j++) {
            attackMask[i][j] = new float[size];
            for (int k = 0; k < size; k++) {
                float i_i = (i - (size-1)/2.0) / ((size-1)/2.0);
                float j_i = (j - (size-1)/2.0) / ((size-1)/2.0);
                float k_i = (k - (size-1)/2.0) / ((size-1)/2.0);
                attackMask[i][j][k] = -(maxVal - sqrt(i_i*i_i + j_i*j_i + k_i*k_i))/(maxVal) * maxStrength;
//                std::cout << attackMask[i][j][k] << "\t";
            }
//            std::cout << std::endl;
        }
//        std::cout << std::endl;
    }
}

void RockErosion::Apply(Voxel* main_v, bool addingMatterMode, bool applyRemeshing) {
    for (int x = -size/2; x < size/2; x++) {
        for (int y = -size/2; y < size/2; y++) {
            for (int z = -size/2; z < size/2; z++){
                Voxel* v = main_v->parent->parent->getVoxel(main_v->globalPos() + Vector3(x, y, z) + Vector3(.5, .5, .5));
                if(v != nullptr) {
                    v->manual_isosurface += this->attackMask[x+size/2][y+size/2][z+size/2] * (addingMatterMode ? -1 : 1);
                    v->manual_isosurface = std::max(v->manual_isosurface, -2.0f);
                    v->manual_isosurface = std::min(v->manual_isosurface, 2.0f);
                    v->parent->needRemeshing = true;
                }
            }
        }
    }
    if (applyRemeshing)
        main_v->parent->parent->remeshAll();
}
