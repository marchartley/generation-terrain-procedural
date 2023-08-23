#include "TerrainSavingInterface.h"
#include "Interface/InterfaceUtils.h"

TerrainSavingInterface::TerrainSavingInterface(QWidget *parent):
    ActionInterface("TerrainSavingInterface", parent)
{

}

void TerrainSavingInterface::display(const Vector3& camPos)
{
    // Nothing to display
}

void TerrainSavingInterface::replay(nlohmann::json action)
{
    // Nothing to replay
}

QLayout *TerrainSavingInterface::createGUI()
{
    QVBoxLayout* layout = new QVBoxLayout();
    QLabel* selectedFilenameLabel = new QLabel((this->mainFilename == "" ? "No file selected" : QString::fromStdString(this->mainFilename).split("/").back().split("\\").back()));
    QPushButton* fileSelectionButton = new QPushButton("...");

    QCheckBox* saveHeightmapCheck = new QCheckBox("Heightmap");
    QCheckBox* saveVoxelsCheck = new QCheckBox("Voxels");
    QCheckBox* saveLayersCheck = new QCheckBox("Layers");

    QPushButton* saveButton = new QPushButton("Save");

    layout->addWidget(createHorizontalGroup({selectedFilenameLabel, fileSelectionButton}));
    layout->addWidget(createVerticalGroup({
                                              saveHeightmapCheck,
                                              saveVoxelsCheck,
                                              saveLayersCheck
                                          }));
    layout->addWidget(saveButton);

    QObject::connect(saveHeightmapCheck, &QCheckBox::toggled, this, [this](bool checked) { this->saveHeightmap = checked; });
    QObject::connect(saveVoxelsCheck, &QCheckBox::toggled, this, [this](bool checked) { this->saveVoxels = checked; });
    QObject::connect(saveLayersCheck, &QCheckBox::toggled, this, [this](bool checked) { this->saveLayers = checked; });

    QObject::connect(saveButton, &QPushButton::pressed, this, [this]() {
        time_t now = std::time(0);
        tm *gmtm = std::gmtime(&now);
        char s_time[80];
        std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);
        this->saveTerrainGeometry(this->mainFilename + "_" + std::string(s_time));
    });

    QObject::connect(fileSelectionButton, &QPushButton::pressed, this, [this, selectedFilenameLabel]() {
        #ifdef linux
            std::string path = "saved_maps/Geometry/";
        #else
            std::string path = "saved_maps/Geometry/";
        #endif
        QString fileSelection = QFileDialog::getSaveFileName(this, "Saving file", QString::fromStdString(path), "*.json", nullptr, QFileDialog::DontConfirmOverwrite);
        if (!fileSelection.isEmpty()) {
            this->mainFilename = fileSelection.toStdString();
            selectedFilenameLabel->setText(QString::fromStdString(this->mainFilename).split("/").back().split("\\").back());
        }
    });

    saveHeightmapCheck->setChecked(this->saveHeightmap);
    saveVoxelsCheck->setChecked(this->saveVoxels);
    saveLayersCheck->setChecked(this->saveLayers);

    return layout;
}

std::vector<std::string> TerrainSavingInterface::saveTerrainGeometry(std::string filename)
{
    bool verbose = true;
    if (filename == "")
        filename = this->mainFilename;

    std::ofstream file;
    if (verbose)
        std::cout << "Saving geometry..." << std::endl;
    Mesh m;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    std::vector<std::string> filenames;
    if (this->saveHeightmap) {
        start = std::chrono::system_clock::now();
        m = this->heightmap->getGeometry();
        file.open(filename + "-heightmap" + ".stl");
        filenames.push_back(filename + "-heightmap" + ".stl");
        file << m.toSTL();
        file.close();
        end = std::chrono::system_clock::now();
        if (verbose)
            std::cout << "Heightmap in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    }

    if (this->saveVoxels) {
        start = std::chrono::system_clock::now();
        m = this->voxelGrid->getGeometry();
        file.open(filename + "-voxels" + ".stl");
        filenames.push_back(filename + "-voxels" + ".stl");
        file << m.toSTL();
        file.close();
        end = std::chrono::system_clock::now();
        if (verbose)
            std::cout << "Voxels in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    }

    if (this->saveLayers) {
        start = std::chrono::system_clock::now();
        m = this->layerGrid->getGeometry();
        file.open(filename + "-layers" + ".stl");
        filenames.push_back(filename + "-layers" + ".stl");
        file << m.toSTL();
        file.close();
        end = std::chrono::system_clock::now();
        if (verbose)
            std::cout << "Layers in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    }

    if (verbose)
        std::cout << "Done." << std::endl;
    return filenames;
}

void TerrainSavingInterface::quickSaveAt(std::string folderName, std::string filePrefix, bool heightmap, bool voxels, bool layers)
{
    std::ofstream file;
    Mesh m;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    if (heightmap) {
        start = std::chrono::system_clock::now();
        m = this->heightmap->getGeometry();
        file.open(folderName + "/" + filePrefix + "-heightmap" + ".stl");
        file << m.toSTL();
        file.close();
        end = std::chrono::system_clock::now();
    }

    if (voxels) {
        start = std::chrono::system_clock::now();
        m = this->voxelGrid->getGeometry();
        file.open(folderName + "/" + filePrefix + "-voxels" + ".stl");
        file << m.toSTL();
        file.close();
        end = std::chrono::system_clock::now();
    }

    if (layers) {
        start = std::chrono::system_clock::now();
        m = this->layerGrid->getGeometry();
        file.open(folderName + "/" + filePrefix + "-layers" + ".stl");
        file << m.toSTL();
        file.close();
        end = std::chrono::system_clock::now();
    }
}
