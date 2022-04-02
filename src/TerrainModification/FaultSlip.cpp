#include "FaultSlip.h"

#include "DataStructure/Matrix3.h"
#include "Utils/Utils.h"

FaultSlip::FaultSlip()
{

}

void FaultSlip::Apply(std::shared_ptr<VoxelGrid> grid, bool applyRemeshing)
{
    Matrix3<float> voxelValues = grid->getVoxelValues();
    Matrix3<float> transportMatrix(voxelValues.sizeX, voxelValues.sizeY, voxelValues.sizeZ, 0);
    transportMatrix.raiseErrorOnBadCoord = false;
    voxelValues.raiseErrorOnBadCoord = false;
    for (int x = 0; x < voxelValues.sizeX; x++) {
        for (int y = 0; y < voxelValues.sizeY; y++) {
            for (int z = 0; z < voxelValues.sizeZ; z++) {
                Vector3 pos(x, y, z);
                bool voxelIsConcerned =
                        tetrahedronSignedVolume(this->firstPointFault, this->secondPointFault,
                                                this->firstPointFault+this->slippingDirection, pos) > 0 == this->positiveSideFalling;
                if (voxelIsConcerned) {
                    transportMatrix.at(pos) -= voxelValues.at(pos);
                    transportMatrix.at(pos + this->slippingDirection * this->slippingDistance) += voxelValues.at(pos) /*- voxelValues.at(pos + this->slippingDirection * this->slippingDistance)*/;
                }
            }
        }
    }
    grid->applyModification(transportMatrix);
    if (applyRemeshing)
        grid->remeshAll();
}
