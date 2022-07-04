#include "BiomeInterface.h"
#include "TerrainModification/RockErosion.h"
#include "Utils/Voronoi.h"

BiomeInterface::BiomeInterface(QWidget* parent)
    : ActionInterface("biome-generation", parent)
{

}

void BiomeInterface::display()
{
    for (size_t i = 0; i < selectionPlanes.size(); i++) {
        auto& selectionPlane = selectionPlanes[i];
        if (selectionPlane.shader != nullptr) {
            Vector3 color = HSVtoRGB(i / (float)(selectionPlanes.size() - 1), 1.f, 1.f);
            selectionPlane.shader->setVector("color", std::vector<float>{color.x, color.y, color.z, .5f});
        }
        selectionPlane.display();
    }
}


Vector3 getSurfacePosition(std::shared_ptr<VoxelGrid> grid, Vector3 pos) {
    pos.x = clamp(pos.x, 0.f, grid->getDimensions().x - 1);
    pos.y = clamp(pos.y, 0.f, grid->getDimensions().y - 1);
    pos.z = std::max(pos.z, 0.f); // In case of small imprecision
    while (grid->getVoxelValue(pos) > 0) {
        pos += Vector3(0, 0, 1); // Move the position one voxel heigher
        if (!grid->contains(pos)) { // If it gets too high, the whole column must be filled, I guess we should cancel it...
            pos.z = 0;
            break;
        }
    }
    return pos;
}
std::shared_ptr<BiomeInstance> recursivelyCreateBiome(nlohmann::json json_content, Vector3 biomePosition, ShapeCurve area) {
    std::string biomeClass = json_content.at("class").get<std::string>();
    // Should be able to retrieve the parameters of the biome...
    std::shared_ptr<BiomeInstance> instance = std::make_shared<BiomeInstance>(BiomeInstance::fromClass(biomeClass));
    instance->position = biomePosition;
    instance->area = area;
    auto children = json_content.at("children");
    Voronoi diagram(children.size(), area);
    std::vector<BSpline> subarea_borders = diagram.solve();
    for (size_t i = 0; i < children.size(); i++) {
        std::shared_ptr<BiomeInstance> childBiome = recursivelyCreateBiome(children[i], diagram.pointset[i], subarea_borders[i]);
        childBiome->parent = instance;
        instance->instances.push_back(childBiome);
    }
    return instance;
}
Matrix3<float> archeTunnel(BSpline path, float size, float strength, bool addingMatter, std::shared_ptr<VoxelGrid> grid) {
    Matrix3<float> erosionMatrix(grid->getDimensions());
    float nb_points_on_path = path.length() / (size/5.f);
    RockErosion rock(size, strength);
    for (const auto& pos : path.getPath(nb_points_on_path)) {
        erosionMatrix = rock.computeErosionMatrix(erosionMatrix, pos);
    }
    erosionMatrix = erosionMatrix.abs();
    erosionMatrix.toDistanceMap();
    erosionMatrix.normalize();
    for (float& m : erosionMatrix) {
        m = interpolation::linear(m, 0.f, 1.0) * strength * (addingMatter ? 1.f : -1.f);
//        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
    }
    return erosionMatrix;
}
void BiomeInterface::generateBiomes(std::shared_ptr<BiomeInstance> predefinedBiomeInstance)
{
    this->heightmap->heights = Matrix3<float>(124, 124, 1, 20.f);
    this->heightmap->maxHeight = 40.f;

    this->voxelGrid->from2DGrid(*this->heightmap);
    this->voxelGrid->fromIsoData();

    ShapeCurve terrainArea({
                   heightmap->heights.getDimensions() * Vector3(0, 0, 0),
                   heightmap->heights.getDimensions() * Vector3(1, 0, 0),
                   heightmap->heights.getDimensions() * Vector3(1, 1, 0),
                   heightmap->heights.getDimensions() * Vector3(0, 1, 0)
               });
    Vector3 initialSpawn = heightmap->heights.getDimensions() / 2.f + Vector3(10, 0, 0);

    // All the biome hierarchy is created from the json
    // Each biome has a class name, a position and his biome children
    if (predefinedBiomeInstance != nullptr)
        rootBiome = predefinedBiomeInstance;
    else
        rootBiome = std::make_shared<BiomeInstance>(biomeModel.createInstance(initialSpawn, terrainArea));

    heightmap->biomeIndices = Matrix3<int>(heightmap->heights.getDimensions(), 0);
    std::vector<std::shared_ptr<BiomeInstance>> biomeQueue;
    std::vector<std::shared_ptr<BiomeInstance>> sortedBiomes;

    biomeQueue.push_back(rootBiome);

    int biomeID = 0;
    while (!biomeQueue.empty()) {
        std::shared_ptr<BiomeInstance> current = biomeQueue.front();
        biomeQueue.erase(biomeQueue.begin());
        if (!terrainArea.inside(current->position * Vector3(1, 1, 0), true)) {
            std::cout << "Biome #" << current->instanceID << " (" << current->classname << ") is outside the terrain" << std::endl;
            continue;
        }
        if (!current->valid) {
            std::cout << "Biome #" << current->instanceID << " (" << current->classname << ") is not valid" << std::endl;
            continue;
        }
        if (current->area.computeArea() == 0) {
            std::cout << "Biome #" << current->instanceID << " (" << current->classname << ") has no space : " << current->area << std::endl;
            continue;
        }
        sortedBiomes.push_back(current);
        biomeQueue.insert(biomeQueue.end(), current->instances.begin(), current->instances.end());
    }
    biomeQueue.push_back(rootBiome);

    // Sort the biomes by level, so that modifications are done top-to-bottom
    std::sort(sortedBiomes.begin(),
              sortedBiomes.end(),
              [](const auto& a, const auto& b) -> bool {
        return a->getLevel() < b->getLevel();
    });

    possibleBiomeInstances.clear();
    selectedBiomeIDs.clear();
    for (auto& current : sortedBiomes) {
        ShapeCurve area = current->area;
        int level = current->getLevel();
//        area = area.grow(-1); // Shrink the area to be able to see all layers
        Vector3 AABBoxMin, AABBoxMax;
        std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
        std::ostringstream out;
        out << "Checking for " << current->classname << " (#" << biomeID << ") :\n";
        bool atLeastOne = false;
        if (biomeID == 34) {
            int a = 0;
//            area = current->area.grow(-1);
//            std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
        }
        for (int x = AABBoxMin.x; x < AABBoxMax.x; x++) {
            for (int y = AABBoxMin.y; y < AABBoxMax.y; y++) {
                if (area.inside(Vector3(x, y, area.points[0].z)) && heightmap->biomeIndices.checkCoord(Vector3(x, y))) {
                    heightmap->biomeIndices.at(x, y).push_back(biomeID);
                    out << "Found at (" << x << ", " << y << ")\n";
                    atLeastOne = true;
                }
            }
        }
        if (!atLeastOne) {
            out << "Never seen...\n";
        }
        out << "Area was\n";
        for (auto& p : current->area.points)
            out << "- " << p << "\n";
        current->instanceID = biomeID;
        BiomeInstance::registerBiome(current);
        possibleBiomeInstances.push_back(*current);
        biomeID++;
//        if (current->getLevel() > 1)
//            std::cout << out.str() << std::endl;
    }
    // First, level the terrain as the biomes design it
    Matrix3<float> heightChange(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1);
    for (auto& current : sortedBiomes) {
        if (!current->depthShape.points.empty()) {
            ShapeCurve area = current->area;
            Vector3 AABBoxMin, AABBoxMax;
            std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
            Matrix3<float> falloff2D(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1, 0.f);
            for (size_t i = 0; i < falloff2D.size(); i++) {
                Vector3 pos = falloff2D.getCoordAsVector3(i);
                if(AABBoxMin.x <= pos.x && pos.x <= AABBoxMax.x && AABBoxMin.y <= pos.y && pos.y <= AABBoxMax.y) {
                    float dist = area.estimateDistanceFrom(pos);
                    falloff2D[i] = dist < 0 ? 1.f : 0.f;
                }
            }
            float desiredDepth = current->depthShape.getPoint(.5f).y;
            Matrix3<float> newHeight = falloff2D.binarize() * -desiredDepth/2.f;
            heightChange += newHeight;
        }
    }
    heightChange = heightChange.meanSmooth(15, 15, 1);
    voxelGrid->add2DHeightModification(heightChange, 1.5f);

    // Now add the primitives on top
    for (auto& current : sortedBiomes) {
        Vector3 pos = current->position;
//        Vector3 surfacePos = getSurfacePosition(voxelGrid, pos);
        ShapeCurve area = current->area;
        Vector3 areaSize = area.containingBoxSize() + Vector3(0, 0, heightmap->maxHeight);
        Vector3 AABBoxMin, AABBoxMax;
        std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
        Matrix3<float> falloff2D(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1, 1.f);
        Matrix3<float> falloff3D(voxelGrid->getDimensions(), 1.f);
        for (size_t i = 0; i < falloff2D.size(); i++) {
            Vector3 pos = falloff2D.getCoordAsVector3(i);
            if(AABBoxMin.x <= pos.x && pos.x <= AABBoxMax.x && AABBoxMin.y <= pos.y && pos.y <= AABBoxMax.y) {
                float dist = area.estimateDistanceFrom(pos);
                falloff2D[i] = dist < 0 ? 1.f - interpolation::fault_distance(-dist, 0.1f) : 0.f;
            }
            for (int z = 0; z < falloff3D.sizeZ; z++) {
                falloff3D.at(pos.x, pos.y, pos.z) = falloff2D[i];
            }
        }
//        Matrix3<float> subFalloff2D

        Matrix3<float> modifications(voxelGrid->getDimensions(), 0.f);
        Matrix3<float> heightmapModifier(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1, 0.f);

        bool modif3D = false;
        bool modif2D = false;

        if (current->classname == "arche") {
            // Create an arch from the two "point" children
            Vector3 point1 = getSurfacePosition(voxelGrid, current->instances[0]->position); // Start
            Vector3 point4 = getSurfacePosition(voxelGrid, current->instances[1]->position); // End
            Vector3 point2 = point1 + (point4 - point1) * .3f + Vector3(0, 0, (point4 - point1).norm() / 2.f); // Midpoint1
            Vector3 point3 = point1 + (point4 - point1) * .6f + Vector3(0, 0, (point4 - point1).norm() / 2.f); // Midpoint2
            modifications = archeTunnel(BSpline({point1, point2, point3, point4}), 5.f, 10.f, true, voxelGrid);
            modif3D = true;
        } else if (current->classname == "patate-corail") {
            float radius = current->area.containingBoxSize().norm() / 2.f;
            Vector3 patatePosition = getSurfacePosition(voxelGrid, current->position) + Vector3(0, 0, radius/10.f);
            modifications = archeTunnel(BSpline({patatePosition, patatePosition + Vector3(0, 0, 0.001f)}), radius, 2.f, true, voxelGrid);
            modif3D = true;
        } else if (current->classname == "tranchee" || current->classname == "passe-corail") {
            // Dig a tunnel on the surface
            Vector3 start = getSurfacePosition(voxelGrid, current->instances[0]->position); // Start
            Vector3 end = getSurfacePosition(voxelGrid, current->instances[1]->position); // End
            Vector3 midpoint1 = getSurfacePosition(voxelGrid, start + (end - start) * .4f + (end - start).cross(Vector3(0, 0, 1)) * .1f);
            Vector3 midpoint2 = getSurfacePosition(voxelGrid, start + (end - start) * .6f + (end - start).cross(Vector3(0, 0, 1)) * -.1f);
            modifications = archeTunnel(BSpline({start, midpoint1, midpoint2, end}), 8.f, 3.f, false, voxelGrid);
            modif3D = true;
        } else if (current->classname == "mur-corail") {
            // Add a small wall depending on the "point" children
            Vector3 start = getSurfacePosition(voxelGrid, current->instances[0]->position); // Start
            Vector3 end = getSurfacePosition(voxelGrid, current->instances[1]->position); // End
            Vector3 gradient = voxelGrid->getVoxelValues().gradient(getSurfacePosition(voxelGrid, (start + end) * .5f));
            Vector3 midpoint = (start + end) * .5f + gradient * (end - start).norm() * .2f;
            modifications = archeTunnel(BSpline({start, midpoint, end}), 5.f, 3.f, true, voxelGrid);
            modif3D = true;
        } else if (current->classname == "point") {
            // Nothing to do, this is the smallest primitive
        } else {
//            std::cout << "How the fuck did I get here? Class was " << current->classname << std::endl;
        }

        if (modif2D)
            voxelGrid->add2DHeightModification(heightmapModifier, 10.f);
        if (modif3D)
            voxelGrid->applyModification(modifications * falloff3D);
    }

    // Yeah whatever... Just do the translation once again
//    voxelGrid->saveState();
    this->heightmap->fromVoxelGrid(*voxelGrid);

/*
    for (auto& tuple : BiomeInstance::instancedBiomes) {
        auto& biome = tuple.second;
        if (biome->instanceID == 34) {
            std::cout << biome->area.inside(Vector3(15, 42, 0) * Vector3(1, 1, 0)) << std::endl;
            std::cout << biome->area.inside(Vector3(15, 42, 0) * Vector3(1, 1, 0)) << std::endl;
        }
//            std::cout << "Distance to biome #" << biome->instanceID << " (" << biome->classname << ") : ";
//            std::cout << biome->area.estimateDistanceFrom(mousePosInMap * Vector3(1, 1, 0)) << " ";
//            std::cout << "(" << (biome->area.inside(mousePosInMap * Vector3(1, 1, 0)) ? "contained" : "outside") << ")" << std::endl;
    }*/
}

void BiomeInterface::replaceBiome(std::shared_ptr<BiomeInstance> biomeToReplace, BiomeInstance newBiome)
{
    auto previousParent = biomeToReplace->parent;
    Vector3 previousPos = biomeToReplace->position;
    ShapeCurve previousArea = biomeToReplace->area;
    int previousID = biomeToReplace->instanceID;

    *biomeToReplace = *(newBiome.clone(previousArea));
    for (auto& child : biomeToReplace->instances) {
        child->parent = biomeToReplace;
    }

    biomeToReplace->parent = previousParent;
    biomeToReplace->position = previousPos;
//    biomeToReplace->area = previousArea;
    biomeToReplace->instanceID = previousID;

    biomeToReplace->completeIfNeeded();
}

void BiomeInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}
void BiomeInterface::affectHeightmap(std::shared_ptr<Grid> heightmap)
{
    this->heightmap = heightmap;
}

void BiomeInterface::replay(nlohmann::json action)
{

}

void BiomeInterface::mouseMoveEvent(QMouseEvent* event)
{

}
void BiomeInterface::keyPressEvent(QKeyEvent* event)
{

}
void BiomeInterface::keyReleaseEvent(QKeyEvent* event)
{

}
void BiomeInterface::wheelEvent(QWheelEvent* event)
{

}

void BiomeInterface::mousePressEvent(QMouseEvent *event)
{
//    if (event->)
}
void BiomeInterface::mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent* event)
{
    if (this->isVisible()) {
        if (!mouseInMap) return;

        this->selectedBiomeIDs = this->heightmap->biomeIndices.at(mousePosInMap * Vector3(1, 1, 0));
        updateBiomeSelectionGui();
    }
}

void BiomeInterface::updateSelectionPlaneToFitBiome(int biomeID, int planeIndex)
{
    if (biomeID == -1) {
        // Remove all triangles
//        selectionPlane.fromArray(std::vector<float>{});
    } else {
//        std::cout << "Displaing surface of " << BiomeInstance::instancedBiomes[biomeID]->classname << std::endl;
        int level = BiomeInstance::instancedBiomes[biomeID]->getLevel();
        ShapeCurve biomeArea = BiomeInstance::instancedBiomes[biomeID]->area;
        std::vector<Vector3> upperPoints, lowerPoints;
        float maxHeight = std::numeric_limits<float>::min(),
                minHeight = std::numeric_limits<float>::max();
        for (auto& point : biomeArea.points) {
            Vector3 surfacePoint = getSurfacePosition(voxelGrid, point);
            upperPoints.push_back(surfacePoint);
            lowerPoints.push_back(surfacePoint);
            maxHeight = std::max(maxHeight, surfacePoint.z);
            minHeight = std::min(maxHeight, surfacePoint.z);
        }
        for (size_t i = 0; i < upperPoints.size(); i++) {
            upperPoints[i].z = maxHeight + 5;
            lowerPoints[i].z = minHeight - 5;
        }
        std::vector<Vector3> vertices;
        for (size_t i = 0; i < upperPoints.size(); i++) {
            size_t next_i = (i + 1) % upperPoints.size();
            vertices.insert(vertices.end(), {
                                upperPoints[i],
                                lowerPoints[i],
                                upperPoints[next_i],

                                lowerPoints[i],
                                lowerPoints[next_i],
                                upperPoints[next_i]
                            });
        }
        selectionPlanes[planeIndex].fromArray(vertices);
        selectionPlanes[planeIndex].cullFace = false;
    }
}

void BiomeInterface::displayUniqueSelection(int selectionIndex)
{
    if (this->selectedBiomeIDs.size() > selectionIndex) {
        this->selectionPlanes.resize(1);
        this->updateSelectionPlaneToFitBiome(this->selectedBiomeIDs[selectionIndex], 0);
    } else {
        this->selectionPlanes.clear();
    }
}

QLayout* BiomeInterface::createGUI()
{
//    if (this->layout != nullptr) return this->layout;

    layout = new QVBoxLayout();
    biomeSelectionGui = new QListWidget;
//    biomeSelectionGui = new QGroupBox("Selection");
//    biomeSelectionGuiLayout = new QVBoxLayout();

    QPushButton* regenerationButton = new QPushButton("Regenerer");
    QPushButton* interchangeBiomeButton = new QPushButton("Changer un biome...");

//    biomeSelectionGui->setLayout(this->biomeSelectionGuiLayout);
    layout->addWidget(biomeSelectionGui);
    layout->addWidget(interchangeBiomeButton);
    layout->addWidget(regenerationButton);
    updateBiomeSelectionGui();

    QObject::connect(regenerationButton, &QPushButton::pressed, this, [&]() -> void { this->generateBiomes(rootBiome); });
    QObject::connect(biomeSelectionGui, &QListWidget::currentRowChanged, this, &BiomeInterface::displayUniqueSelection);
    QObject::connect(interchangeBiomeButton, &QPushButton::pressed, this, [&]() -> void {
        BiomeReplacementDialog dialog(this);
//        this->replaceDialog->allAvailableBiomes->clear();
        for (auto& biome : possibleBiomeInstances) {
//            std::ostringstream oss;
//            oss << biome.classname << (biome.depthShape.length() > 0 ? " - profondeur : " + std::to_string(biome.depthShape.getPoint(.5f).y) + "m" : "") << " (ID : " << biome.instanceID << ")";
            dialog.allAvailableBiomes->addItem(QString::fromStdString(biome.classname));
        }
        int selectedIndex = this->biomeSelectionGui->currentRow();
        dialog.show();
        int replacementIndex = dialog.exec();
        replacementIndex = tempIndex;

//        std::cout << "Biome #" << selectedIndex << " will be replaced by #" << replacementIndex << std::endl;
//        std::cout << "(" << BiomeInstance::instancedBiomes[selectedIndex]->classname << ") by (" << possibleBiomeInstances[replacementIndex].classname << ")" << std::endl;

        if (selectedIndex > -1 && replacementIndex > -1)
            this->replaceBiome(BiomeInstance::instancedBiomes[this->selectedBiomeIDs[selectedIndex]], possibleBiomeInstances[replacementIndex]);
    });

    this->replaceDialog = new BiomeReplacementDialog(this);
    return layout;
}

void BiomeInterface::hide()
{
    for (auto& selectionPlane : selectionPlanes)
        selectionPlane.hide();
    CustomInteractiveObject::hide();
}

void BiomeInterface::show()
{
    for (auto& selectionPlane : selectionPlanes)
        selectionPlane.show();
    CustomInteractiveObject::show();
}

void BiomeInterface::setBindings()
{

}

void BiomeInterface::updateBiomeSelectionGui()
{
//    qDeleteAll(this->biomeSelectionGui->findChildren<QWidget *>(QString(), Qt::FindDirectChildrenOnly));
    biomeSelectionGui->clear();

    for (size_t i = 0; i < this->selectedBiomeIDs.size(); i++) {
        int biomeID = this->selectedBiomeIDs[i];
        if (BiomeInstance::instancedBiomes.find(biomeID) == BiomeInstance::instancedBiomes.end())
            continue; // Biome not registered... This shouldn't occur.
        auto biome = BiomeInstance::instancedBiomes[biomeID];
//        std::ostringstream oss;
//        oss << biome->classname << (biome->depthShape.length() > 0 ? " - profondeur : " + std::to_string(biome->depthShape.getPoint(.5f).y) + "m" : "") << " (ID : " << biome->instanceID << ")";
        biomeSelectionGui->addItem(QString::fromStdString(biome->classname));
//        QLabel* nameLabel = new QLabel(QString::fromStdString(biome->classname));
//        Vector3 textColor = HSVtoRGB(i / (float)(selectedBiomeIDs.size() - 1), 1, 1) * 255.f;
//        std::ostringstream oss;
//        oss << "QLabel{color: rgb(" << int(textColor.x) << "," << int(textColor.y) << "," << int(textColor.z) << ");}";
//        nameLabel->setStyleSheet(oss.str().c_str());
//        this->biomeSelectionGuiLayout->addWidget(nameLabel);
    }

//    if (this->selectedBiomeIDs.size() > 3) {
    this->selectionPlanes.resize(this->selectedBiomeIDs.size());
    for (size_t i = 0; i < this->selectedBiomeIDs.size(); i++)
        this->updateSelectionPlaneToFitBiome(this->selectedBiomeIDs[i], i);
//    }
}

BiomeReplacementDialog::BiomeReplacementDialog(BiomeInterface* caller)
    : QDialog(), caller(caller)
{
    QVBoxLayout * vBoxLayout = new QVBoxLayout(this);
    allAvailableBiomes = new QListWidget(this);
    cancelButton = new QPushButton("Annuler", this);
    validButton = new QPushButton("Confirmer", this);

    vBoxLayout->addWidget(allAvailableBiomes);
    vBoxLayout->addWidget(cancelButton);
    vBoxLayout->addWidget(validButton);

    QObject::connect(cancelButton, &QPushButton::pressed, this, &BiomeReplacementDialog::cancel);
    QObject::connect(validButton, &QPushButton::pressed, this, &BiomeReplacementDialog::confirm);

    setLayout(vBoxLayout);
    setSizeGripEnabled(true);
}

void BiomeReplacementDialog::open()
{
    QDialog::open();
}

void BiomeReplacementDialog::cancel()
{
    selectedBiomeIndex = -1;
    setResult(-1);
    caller->tempIndex = -1;
    this->close();
}

void BiomeReplacementDialog::confirm()
{
    if (allAvailableBiomes->currentRow()) {
        selectedBiomeIndex = allAvailableBiomes->currentRow();
        setResult(selectedBiomeIndex);
        caller->tempIndex = selectedBiomeIndex;
        this->close();
    } else {
        cancel();
    }
}
