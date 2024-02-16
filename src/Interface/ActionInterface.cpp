#include "ActionInterface.h"

ActionInterface::ActionInterface(std::string actionTypeName,
                                 std::string interfaceName,
                                 std::string interfaceType,
                                 std::string mainActionDescription,
                                 std::string mainActionButtonLogo,  QWidget *parent)
    : CustomInteractiveObject(parent),
      actionType(actionTypeName),
      interfaceName(interfaceName),
      interfaceType(interfaceType),
      mainActionDescription(mainActionDescription),
      mainActionButtonLogo(mainActionButtonLogo)
{
    this->actionType = simplify(actionType);
    this->interfaceType = simplify(interfaceType);
}

std::shared_ptr<ActionInterface> ActionInterface::findOtherInterface(std::string name) const
{
    if (!viewer) {
        std::cerr << "The viewer from this interface (" << this->interfaceName << ") is not defined, cannot fetch interface '" << name << "'..." << std::endl;
        return nullptr;
    }
    if (viewer->interfaces.count(name) == 0) {
        std::cerr << "The viewer has no interface '" << name << "'... Here are all available interfaces :\n";
        for (auto& nameAndInterface : viewer->interfaces)
            std::cerr << "- " << nameAndInterface.first << "\n";
        std::cerr << std::endl;
        return nullptr;
    }
    return viewer->interfaces[name];
}

//ActionInterface::~ActionInterface()
//{
//}

