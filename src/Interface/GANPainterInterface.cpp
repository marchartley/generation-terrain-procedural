#include "GANPainterInterface.h"
#include "Graphics/GanUIPainter.h"
#include "Interface/TerrainGenerationInterface.h"

GANPainterInterface::GANPainterInterface(QWidget *parent)
    : ActionInterface("ganpainter", "cGAN Painter UI", "digging", "Draw islands using GANs", "open_button.png", parent)
{
    ganProc = new QProcess();
    ganProc->setWorkingDirectory("..");
    ganProc->start("python", QStringList() << "-u" << "generation-terrain-procedural/Python_tests/pytorch-CycleGAN-and-pix2pix/testSingleLabelsToHeightmaps.py" << "--fromCpp" << "--input" << "generation-terrain-procedural/cGANdata/input_label.png" << "--output" << "generation-terrain-procedural/cGANdata/result_height.png" << "--name" << "labels_to_terrain_pacific_graphics_2025_larger_reefs" << "--model" << "pix2pix" << "--direction" << "AtoB");

    if (!ganProc->waitForStarted(-1)) {
        std::cout << "Could not start cGAN process" << std::endl;
        qDebug() << "Error: " << ganProc->errorString();
    }
}

GANPainterInterface::~GANPainterInterface()
{
    ganProc->write("close\n");
    if (!ganProc->waitForFinished()) {
        std::cout << "Could not finish cGAN process" << std::endl;
    } else {
        std::cout << "Closed Python interpreter for cGAN process" << std::endl;
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
    QObject::connect(GanUIPainter::get("cGanPainter"), &GanUIPainter::clickedOnImage, this, [&](const Vector3& pos, const Vector3& value) {
        if (!pos.isValid()) { // On release
            GanUIPainter::get("cGanPainter")->saveFig("cGANdata/input_label.png");
            ganProc->write("infer\n");
            ganProc->waitForBytesWritten();
            ganProc->waitForReadyRead(2000);
            std::string resultErr = ganProc->readAllStandardError().toStdString();
            if (!resultErr.empty()) {
                this->log("Python [error]: " + resultErr);
            }
            displayProcessTime("Opening file...", [&]() {
                this->heightmap->loadFromHeightmap("cGANdata/result_height.png", this->heightmap->getSizeX(), this->heightmap->getSizeY(), 30.f);
            });
            displayProcessTime("Calling update...", [&]() {
                Q_EMIT this->terrainUpdated();
                Q_EMIT this->updated();
            });
        }
    });
}

QLayout *GANPainterInterface::createGUI()
{
    InterfaceUI* gui = new InterfaceUI(new QVBoxLayout());
    ButtonElement* btn = new ButtonElement("TEST", [&]() {
        GanUIPainter::get("cGanPainter")->show();
    });
    gui->add(btn);
    return gui->get()->layout();
}
