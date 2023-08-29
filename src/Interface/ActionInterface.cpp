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

//ActionInterface::~ActionInterface()
//{
//}

