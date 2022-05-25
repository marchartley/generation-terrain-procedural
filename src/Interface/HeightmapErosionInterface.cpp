#include "HeightmapErosionInterface.h"

#include "Interface/InterfaceUtils.h"

HeightmapErosionInterface::HeightmapErosionInterface(QWidget *parent)
    : ActionInterface("heightmap-erosion", parent)
{
    hydraulicMesh = Mesh({}, true, GL_LINES);
    hydraulicMesh.useIndices = false;
}

void HeightmapErosionInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}

void HeightmapErosionInterface::display()
{
    if (hydraulicMesh.shader != nullptr)
        hydraulicMesh.shader->setVector("color", std::vector<float>({1.f, .0f, .0f, 1.f}));
    hydraulicMesh.display(GL_LINES);
}

void HeightmapErosionInterface::replay(nlohmann::json action)
{
//    if (this->isConcerned(action)) {
//        auto& parameters = action.at("parameters");
//        Vector3 pos = json_to_vec3(parameters.at("position"));
//        Vector3 dir = json_to_vec3(parameters.at("direction"));
//        float size = parameters.at("size").get<float>();
//        int qtt = parameters.at("quantity").get<int>();
//        float strength = parameters.at("strength").get<float>();
//        float randomness = parameters.at("randomness").get<float>();
//        UnderwaterErosion erod(this->voxelGrid, size, strength, qtt);
//        erod.Apply(pos, dir, 0.f, 0.f, randomness, false);
//    }
}

void HeightmapErosionInterface::hydraulicErosion()
{
    if (this->heightmap != nullptr) {
        std::vector<std::vector<Vector3>> traces = heightmap->hydraulicErosion();
        std::vector<Vector3> segments;
        for (const auto& trace : traces) {
            if (trace.size() < 2) continue;
            for (size_t i = 0; i < trace.size() - 1; i++) {
                segments.push_back(trace[i]);
                segments.push_back(trace[i + 1]);
            }
        }
        hydraulicMesh.fromArray(segments);
        Q_EMIT updated();
    }
}

void HeightmapErosionInterface::windErosion()
{
    if (this->heightmap != nullptr) {
        heightmap->windErosion();
        Q_EMIT updated();
    }
}

QLayout *HeightmapErosionInterface::createGUI()
{
    this->erosionLayout = new QHBoxLayout();

    QPushButton* hydraulicErosionButton = new QPushButton("Erosion hydraulique");
    QPushButton* thermicErosionButton = new QPushButton("Erosion thermique");
    QPushButton* windErosionButton = new QPushButton("Erosion de vent");

    erosionLayout->addWidget(createVerticalGroup({
                                                     hydraulicErosionButton,
                                                     thermicErosionButton,
                                                     windErosionButton
                                                 }));

    QObject::connect(hydraulicErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::hydraulicErosion);
//    QObject::connect(thermicErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::thermicErosion);
    QObject::connect(windErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::windErosion);

    return erosionLayout;
}

void HeightmapErosionInterface::show()
{
    this->hydraulicMesh.show();
    ActionInterface::show();
}

void HeightmapErosionInterface::hide()
{
    this->hydraulicMesh.hide();
    ActionInterface::hide();
}
