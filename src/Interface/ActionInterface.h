#ifndef ACTIONINTERFACE_H
#define ACTIONINTERFACE_H

#include "Interface/CustomInteractiveObject.h"
#include "TerrainModification/TerrainAction.h"
#include "TerrainGen/VoxelGrid.h"

class ActionInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    ActionInterface(std::string actionTypeName, QWidget* parent = nullptr);
//    ~ActionInterface();

    virtual void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid) {

    }
    virtual QLayout* createGUI() {
        return nullptr;
    }
    virtual void display() {

    }

    void addTerrainAction(nlohmann::json actionParameters) {
        jsonActionsHistory->push_back(nlohmann::json({
                                                         {"type", this->actionType},
                                                         {"parameters", actionParameters}
                                                     }));
//        saveAllActions();
    }

    void affectSavingFile(std::shared_ptr<std::vector<nlohmann::json>> history, std::shared_ptr<std::fstream> file, std::string filename) {
        this->jsonActionsHistory = history;
        this->savingFile = file;
        this->savingFilename = filename;
    }

    void saveAllActions(std::string filename = "") {
        savingFile->close();
        savingFile->open((filename.empty() ? savingFilename : filename), std::fstream::in | std::fstream::out | std::fstream::trunc);
        *savingFile << nlohmann::json({{"actions", *jsonActionsHistory}}).dump(1, '\t');
        savingFile->flush();
    }

    bool isConcerned(nlohmann::json& action) { return action.contains("type") && action.at("type").get<std::string>() == this->actionType; }

    virtual void replay(nlohmann::json action) = 0;
    std::string actionType;
    std::string savingFilename;
    std::shared_ptr<std::fstream> savingFile;
    std::shared_ptr<std::vector<nlohmann::json>> jsonActionsHistory;

Q_SIGNALS:
    void updated();
};

#include "Utils/BSpline.h"
#include "Karst/KarstHoleProfile.h"
nlohmann::json vec3_to_json(const Vector3& vec);
Vector3 json_to_vec3(nlohmann::json json);
nlohmann::json bspline_to_json(const BSpline& spline);
BSpline json_to_bspline(nlohmann::json json);
nlohmann::json karst_profile_to_json(KarstHoleProfile profile);
KarstHoleProfile json_to_karst_profile(nlohmann::json json);
#endif // ACTIONINTERFACE_H
