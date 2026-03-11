#include "TerrainSavingInterface.h"
#include "GUIElements/InterfaceUtils.h"

TerrainSavingInterface::TerrainSavingInterface(QWidget *parent):
    ActionInterface("terrainsaving", "Terrain geometry saving", "view", "Save the geometry", "save_geometry.png", parent)
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

InterfaceUI* TerrainSavingInterface::createGUI()
{
    auto UI = new InterfaceUI();

    auto selectedFilenameLabel = new LabelElement((this->mainFilename == "" ? "No file selected" : QString::fromStdString(this->mainFilename).split("/").back().split("\\").back().toStdString()));
    auto fileSelectionButton = new ButtonElement("...", [&]() {
        std::string path = "saved_maps/Geometry/";
        QString fileSelection = QFileDialog::getSaveFileName(this, "Saving file", QString::fromStdString(path), "*.json", nullptr, QFileDialog::DontConfirmOverwrite);
        if (!fileSelection.isEmpty()) {
            this->mainFilename = fileSelection.toStdString();
            selectedFilenameLabel->setText(split(split(this->mainFilename, "/").back(), "\\").back());
        }
    });

    auto saveHeightmapCheck = new CheckboxElement("Heightmap", saveHeightmap);
    auto saveVoxelsCheck = new CheckboxElement("Voxels", saveVoxels);
    auto saveLayersCheck = new CheckboxElement("Layers", saveLayers);

    auto saveButton = new ButtonElement("Save", [this]() {
        time_t now = std::time(0);
        tm *gmtm = std::gmtime(&now);
        char s_time[80];
        std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);
        this->saveTerrainGeometry(this->mainFilename + "_" + std::string(s_time));
    });


    UI->add(std::vector<UIElement*>{
        createHorizontalGroupUI({selectedFilenameLabel, fileSelectionButton}),
        saveHeightmapCheck,
        saveVoxelsCheck,
        saveLayersCheck,
        saveButton
    });

    return UI;
}

std::vector<std::string> TerrainSavingInterface::saveTerrainGeometry(const std::string& _filename)
{
    std::string filename = _filename;
    bool verbose = true;
    if (filename == "")
        filename = this->mainFilename;

    std::ofstream file;
    if (verbose)
        std::cout << "Saving geometry..." << std::endl;
    Mesh m;
    float timing = 0;
    std::vector<std::string> filenames;
    if (this->saveHeightmap) {
        displayProcessTime("Saving heightmap geometry... ", [&]() {
            m = this->heightmap->getGeometry();
            file.open(filename + "-heightmap" + ".stl");
            filenames.push_back(filename + "-heightmap" + ".stl");
            file << m.toSTL();
            file.close();
        }, verbose);
    }

    if (this->saveVoxels) {
        displayProcessTime("Saving voxel grid geometry... ", [&]() {
            m = this->voxelGrid->getGeometry();
            file.open(filename + "-voxels" + ".stl");
            filenames.push_back(filename + "-voxels" + ".stl");
            file << m.toSTL();
            file.close();
        }, verbose);
    }

    if (this->saveLayers) {
        displayProcessTime("Saving layers geometry... ", [&]() {
            m = this->layerGrid->getGeometry();
            file.open(filename + "-layers" + ".stl");
            filenames.push_back(filename + "-layers" + ".stl");
            file << m.toSTL();
            file.close();
        }, verbose);
    }
    return filenames;
}

void TerrainSavingInterface::quickSaveAt(const std::string& folderName, std::string filePrefix, bool heightmap, bool voxels, bool layers)
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
