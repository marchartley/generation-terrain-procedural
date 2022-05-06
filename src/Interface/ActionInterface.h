#ifndef ACTIONINTERFACE_H
#define ACTIONINTERFACE_H

#include "Interface/CustomInteractiveObject.h"
#include "TerrainModification/TerrainAction.h"
// #include "Utils/RewritableFile.h"
#include "Utils/BSpline.h"
#include "Karst/KarstHoleProfile.h"

class ActionInterface : public CustomInteractiveObject
{
public:
    ActionInterface(std::string actionTypeName, QWidget* parent = nullptr);
    ~ActionInterface();

    /*
    void addTerrainAction(TerrainAction action) {
        actionsHistory->push_back(action);
    }

    void affectActionsHistory(std::shared_ptr<std::vector<TerrainAction>> history) {
        this->actionsHistory = history;
    }

    void saveAllActions(std::string filename) {
        auto allActionsArray = nlohmann::json::array();
        for (auto& action : *actionsHistory) {
            allActionsArray.push_back(action.serialize());
        }

        std::ofstream file(filename);
        file << nlohmann::json({{"actions", allActionsArray}}).dump(1, '\t');
    }

    void importActionsFromFile(std::string filename) {
        std::ifstream file(filename);
        nlohmann::json json = nlohmann::json::parse(file);
        actionsHistory->clear();
        for (auto& action : json.at("actions")) {
            actionsHistory->push_back(
                        Terrain)
        }
    }


    std::shared_ptr<std::vector<TerrainAction>> actionsHistory;
    */

    void addTerrainAction(nlohmann::json actionParameters) {
        jsonActionsHistory->push_back(nlohmann::json({
                                                         {"type", this->actionType},
                                                         {"parameters", actionParameters}
                                                     }));
        saveAllActions();
    }

    void affectSavingFile(std::shared_ptr<std::vector<nlohmann::json>> history, std::shared_ptr<std::fstream> file, std::string filename) {
        this->jsonActionsHistory = history;
        this->savingFile = file;
        this->savingFilename = filename;
    }

    void saveAllActions() {
        savingFile->close();
        savingFile->open(savingFilename, std::fstream::in | std::fstream::out | std::fstream::trunc);
        *savingFile << nlohmann::json({{"actions", *jsonActionsHistory}}).dump(1, '\t');
        savingFile->flush();
    }

    bool isConcerned(nlohmann::json& action) { return action.contains("type") && action.at("type").get<std::string>() == this->actionType; }

    virtual void replay(nlohmann::json action) = 0;
    std::string actionType;
    std::string savingFilename;
    std::shared_ptr<std::fstream> savingFile;
    std::shared_ptr<std::vector<nlohmann::json>> jsonActionsHistory;
};

nlohmann::json vec3_to_json(const Vector3& vec);
Vector3 json_to_vec3(nlohmann::json json);
nlohmann::json bspline_to_json(const BSpline& spline);
BSpline json_to_bspline(nlohmann::json json);
nlohmann::json karst_profile_to_json(KarstHoleProfile profile);
KarstHoleProfile json_to_karst_profile(nlohmann::json json);
#endif // ACTIONINTERFACE_H
