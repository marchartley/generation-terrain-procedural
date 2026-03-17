#include "GANPainterInterface.h"
#include "GUIElements/GanUIPainter.h"
#include "Interface/TerrainGenerationInterface.h"
#include "DataStructure/Image.h"

GANPainterInterface::GANPainterInterface(QWidget *parent)
    : ActionInterface("ganpainter", "cGAN Painter UI", "digging", "Draw islands using GANs", "open_button.png", parent)
{
    // ganModelNames["volcanic"] = "labels_to_terrain_volcanic_islands";
    // ganModelNames["classic"] = "labels_to_terrain_pacific_graphics_2025";
    // ganModelNames["reefs"] = "labels_to_terrain_pacific_graphics_2025_larger_reefs";
    // ganModelNames["undeformed"] = "labels_to_terrain_features_undeformed";
    // ganModelNames["no dendritic"] = "labels_to_terrain_no_dendritic";

    for (auto& [name, modelName] : ganModelNames) {
        ganProcesses[name] = new QProcess();
        ganModelWeights[name] = 1.f;
        ganComputedHeightfields[name] = GridF();
        auto& process = ganProcesses[name];
        process->setWorkingDirectory("..");

        process->start("python", QStringList() <<
                                "-u" <<
                                "generation-terrain-procedural/Python_tests/pytorch-CycleGAN-and-pix2pix/testSingleLabelsToHeightmaps.py" <<
                                "--fromCpp" <<
                                "--input" << "generation-terrain-procedural/cGANdata/input_label.png" <<
                                "--output" << "generation-terrain-procedural/cGANdata/result_height-" + QString::fromStdString(name) + ".png" <<
                                "--name" << QString::fromStdString(modelName) <<
                                "--model" << "pix2pix" <<
                                "--direction" << "AtoB");

        if (!process->waitForStarted(-1)) {
            std::cout << "Could not start cGAN '" + name + "'  process" << std::endl;
            qDebug() << "Error: " << process->errorString();
        }
    }
}

GANPainterInterface::~GANPainterInterface()
{
    for (auto& [name, process] : ganProcesses) {
        process->write("close\n");
        if (!process->waitForFinished()) {
            std::cout << "Could not finish cGAN '" + name + "' process" << std::endl;
        } else {
            std::cout << "Closed Python interpreter for cGAN '" + name + "'  process" << std::endl;
        }
    }
}

void GANPainterInterface::display(const Vector3 &camPos)
{

}

void GANPainterInterface::replay(nlohmann::json action)
{

}

void GANPainterInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);
    QObject::connect(&GanUIPainter::get("cGanPainter"), &GanUIPainter::clickedOnImage, this, [&](const Vector3& pos, const Vector3& value) {
        if (!pos.isValid()) { // On release
            this->runGANs();
        }
    });
}

void GANPainterInterface::runGANs()
{
    GridV3 labelImage = GanUIPainter::get("cGanPainter").dataModel->getImage() * subsidence;
    Image(labelImage).writeToFile("cGANdata/input_label.png");

    for (auto& [name, process] : ganProcesses) {
        process->write("infer\n");
        process->waitForBytesWritten();
        process->waitForReadyRead(2000);
        std::string resultErr = process->readAllStandardError().toStdString();
        if (!resultErr.empty()) {
            this->log("Python [error]: " + resultErr);
        }
        ganComputedHeightfields[name] = Image::readFromFile("cGANdata/result_height-" + name + ".png").getBwImage();
    }
}

GridV3 GANPainterInterface::computeNewHeightmap()
{
    GridF finalTerrain;
    float totalWeight = 0.f;

    for (auto& [name, heights] : ganComputedHeightfields) {
        auto& thisGanResult = heights;
        if (finalTerrain.empty()) {
            finalTerrain = thisGanResult * ganModelWeights[name];
        } else {
            finalTerrain += thisGanResult * ganModelWeights[name];
        }
        totalWeight += ganModelWeights[name];
        // break;
    }
    finalTerrain /= (255.f * totalWeight);

    displayProcessTime("Updating terrain...", [&]() {
            // this->heightmap->loadFromHeightmap("cGANdata/result_height.png", this->heightmap->getSizeX(), this->heightmap->getSizeY(), 30.f);
            this->heightmap->heights = finalTerrain.resize(this->heightmap->getSizeX(), this->heightmap->getSizeY(), 1) * 30.f;
            dynamic_cast<TerrainGenerationInterface*>(this->findOtherInterface("terraingeneration").get())->heightmapToAll();
        }, false);
    displayProcessTime("Calling update...", [&]() {
            Q_EMIT this->terrainUpdated();
            Q_EMIT this->updated();
        }, false);
    return finalTerrain;
}

GridV3 GANPainterInterface::randomIslandGenerator()
{
    GridV3 labels(256, 256, 1);
    float offset = random_gen::generate(10000.f);
    labels.iterateParallel([&] (const Vector3i& pos) {
        auto& p = pos;
        // float val = (random_gen::generate_perlin(p.x() / 150.f, p.y() / 150.f) * .8f + random_gen::generate_perlin(p.x() / 100.f, p.y() / 100.f) * .3f + random_gen::generate_perlin(p.x() / 50.f, p.y() / 50.f) * .2f) * 2.f - 1.f;
        float val = random_gen::generate_fbm(p.x() / 2.f, p.y() / 2.f, offset) * .5f + .7f;
        val = clamp(sign(val) * std::pow(abs(val), 1.5f), 0.f, .99f);
        val *= clamp(Vector3::distanceToBoundaries(p, Vector3(0, 0, 0), Vector3(256, 256), true), 0.f, 20.f) / 20.f;
        // val = val * std::sin(clamp(std::min({p.x(), p.y(), 256 - p.x(), 256 - p.y()}) / 256.f, 0.f, 1.f) * PI * .5f);
        if (val < .4f) val = 4;
        else if (val < .5f) val = 3;
        else if (val < .7f) val = 2;
        else if (val < .8f) val = 1;
        else val = 0;
        labels[pos] = GanUIPainter::getColorFromIndex(val);
    });
    return labels;
}

InterfaceUI* GANPainterInterface::createGUI()
{
    auto UI = new InterfaceUI();
    auto drawButton = new ButtonElement("Draw", [&]() {
        GanUIPainter::get("cGanPainter").show();
    });
    auto randomButton = new ButtonElement("Random", [&]() {
        GanUIPainter::get("cGANPainter").addImage(this->randomIslandGenerator()).show();
        this->runGANs();
        this->computeNewHeightmap();
    });
    auto subsidenceSlider = new SliderElement("Subsidence", 0.01f, .97f, 0.01f, this->subsidence);
    subsidenceSlider->setOnValueChanged([=](float newVal) { this->runGANs(); this->computeNewHeightmap(); });

    auto modelsUI = new InterfaceUI(InterfaceUI::VERTICAL);

    for (auto& [name, process] : ganProcesses) {
        auto slider = new SliderElement(name, 0.f, 1.f, 0.01f, ganModelWeights[name]);
        slider->setOnValueChanged([=](float newVal) { this->computeNewHeightmap(); });
        modelsUI->add(slider);
    }

    UI->add(std::vector<UIElement*>{drawButton, subsidenceSlider, randomButton, modelsUI});

    return UI;
}
