#include "TerrainComparatorInterface.h"
#include "Interface/FancySlider.h"
#include "Interface/InterfaceUtils.h"
#include "Interface/RangeSlider.h"

#include "DataStructure/Image.h"

TerrainComparatorInterface::TerrainComparatorInterface(QWidget *parent)
    : ActionInterface("terraincomparator", "Terrain Comparator", "digging", "Compare different terrains", "need_to_find_logo.png", parent)
{

}

void TerrainComparatorInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

    this->terrainMesh.displayWithOutlines(std::vector<float> {.1f, .9f, .3f, .5f});

    return ActionInterface::display(camPos);
}

void TerrainComparatorInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);

    const char* vParticleShader = "src/Shaders/particle.vert";
    const char* gParticleShader = "src/Shaders/particle.geom";
    const char* fParticleShader = "src/Shaders/particle.frag";
    const char* vNoShader = "src/Shaders/no_shader.vert";
    const char* fNoShader = "src/Shaders/no_shader.frag";
    const char* vMCShader = "src/Shaders/MarchingCubes.vert";
    const char* fMCShader = "src/Shaders/no_shader.frag";
    const char* gMCShader = "src/Shaders/MarchingCubes.geom";

    this->terrainMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->terrainMesh.useIndices = false;
//    this->terrainMesh = Mesh(std::make_shared<Shader>(vMCShader, fMCShader, gMCShader));
//    this->terrainMesh.useIndices = false;
//    this->otherMeshToDisplay = Mesh(std::make_shared<Shader>(vNoShader, fNoShader)); //(vMCShader, fMCShader, gMCShader));
//    this->otherMeshToDisplay.useIndices = false;
    this->terrainMesh.shader->setVector("color", std::vector<float> {.1f, .9f, .3f, .5f});
}

void TerrainComparatorInterface::replay(nlohmann::json action)
{
//     return ActionInterface::replay(action);
}

void TerrainComparatorInterface::mouseMoveEvent(QMouseEvent *event)
{
    return ActionInterface::mouseMoveEvent(event);
}

void TerrainComparatorInterface::keyPressEvent(QKeyEvent *event)
{
    return ActionInterface::keyPressEvent(event);
}

void TerrainComparatorInterface::keyReleaseEvent(QKeyEvent *event)
{
    return ActionInterface::keyReleaseEvent(event);
}

void TerrainComparatorInterface::wheelEvent(QWheelEvent *event)
{
    return ActionInterface::wheelEvent(event);
}

void TerrainComparatorInterface::mousePressEvent(QMouseEvent *event)
{
    return ActionInterface::mousePressEvent(event);
}

QLayout *TerrainComparatorInterface::createGUI()
{
    QLayout* layout = new QVBoxLayout;

    auto updateButton = new ButtonElement("Update", [&]() { this->updateStuff(); });
    auto unionCheckbox = new RadioButtonElement("Union", displayUnion);
    auto intersectionCheckbox = new RadioButtonElement("Intersection", displayIntersection);
    auto substractionCheckbox = new RadioButtonElement("Substract A-B", displaySubstractionAB);
    auto substraction2Checkbox = new RadioButtonElement("Substract B-A", displaySubstractionBA);
    auto interpolationSlider = new SliderElement("Interpolate", 0.f, 1.f, 0.01f);
    interpolationSlider->setOnValueChanged([&](float newValue) { interpolate(newValue); });

    unionCheckbox->setOnChecked([&](bool) { updateStuff(); });
    intersectionCheckbox->setOnChecked([&](bool) { updateStuff(); });
    substractionCheckbox->setOnChecked([&](bool) { updateStuff(); });
    substraction2Checkbox->setOnChecked([&](bool) { updateStuff(); });

    auto ui = new InterfaceUI(new QVBoxLayout, "Form");
    ui->add({
                unionCheckbox,
                intersectionCheckbox,
                substractionCheckbox,
                substraction2Checkbox,
                interpolationSlider
//                updateButton,
                 });
    layout->addWidget(ui->get());
    return layout;
}

void TerrainComparatorInterface::hide()
{
    ActionInterface::hide();
}

void TerrainComparatorInterface::show()
{
    this->updateStuff();
    return ActionInterface::show();
}

void TerrainComparatorInterface::mouseClickedOnMapEvent(const Vector3& mousePosInMap, bool mouseInMap, QMouseEvent *event, TerrainModel *model)
{
    return ActionInterface::mouseClickedOnMapEvent(mousePosInMap, mouseInMap, event, model);
}

void TerrainComparatorInterface::afterTerrainUpdated()
{

}

void TerrainComparatorInterface::afterWaterLevelChanged()
{

}

void TerrainComparatorInterface::interpolate(float t)
{
    const GridF& A = displayedVoxels;
    const GridF& B = voxelsFromHeightmap._cachedVoxelValues;

    this->terrainMesh.fromArray(Mesh::applyMarchingCubes(A * t + B * (1.f - t)).vertexArray);
    Q_EMIT updated();
}

void TerrainComparatorInterface::updateStuff()
{
    this->voxelsFromHeightmap._cachedVoxelValues = this->voxelGrid->_cachedVoxelValues; // To get same dimensions
    this->voxelsFromHeightmap.from2DGrid(*this->heightmap);
    displayedVoxels = voxelsFromHeightmap._cachedVoxelValues.binarize();
    GridF initial = voxelGrid->_cachedVoxelValues.binarize();

    if (displayUnion) {
        displayedVoxels.iterateParallel([&](size_t i) {
            displayedVoxels[i] = std::max(displayedVoxels[i], initial[i]);
        });
    } else if (displayIntersection) {
        displayedVoxels.iterateParallel([&](size_t i) {
            displayedVoxels[i] = (displayedVoxels[i] > 0 && initial[i] > 0 ? std::max(displayedVoxels[i], initial[i]) : std::min(displayedVoxels[i], initial[i]));
        });
    } else if (displaySubstractionAB) {
        displayedVoxels.iterateParallel([&](size_t i) {
            displayedVoxels[i] = (displayedVoxels[i] > 0 && initial[i] > 0 ? -std::max(std::abs(displayedVoxels[i] - initial[i]), 0.01f) : displayedVoxels[i]);
        });
    } else if (displaySubstractionBA) {
        displayedVoxels.iterateParallel([&](size_t i) {
            displayedVoxels[i] = (displayedVoxels[i] > 0 && initial[i] > 0 ? -std::max(std::abs(displayedVoxels[i] - initial[i]), 0.01f) : initial[i]);
        });
    }
//    displayedVoxels -= .5f;
    this->terrainMesh.fromArray(Mesh::applyMarchingCubes(displayedVoxels.meanSmooth().meanSmooth()).vertexArray);
    Q_EMIT updated();
}
