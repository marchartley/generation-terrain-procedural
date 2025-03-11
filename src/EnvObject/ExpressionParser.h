#ifndef EXPRESSIONPARSER_H
#define EXPRESSIONPARSER_H

#include "DataStructure/Vector3.h"
#include <vector>
#include <string>
#include <map>
#include <functional>
#include <variant>
#include <cmath>
#include <memory>
#include <set>


struct TokenList;

using Token = std::variant<std::string, std::shared_ptr<TokenList>, Vector3>;
struct TokenList : public std::vector<Token> {
    using std::vector<Token>::vector;
};

using Variable = std::variant<float, Vector3>;
using VariableMap = std::map<std::string, Variable>;

using UnaryFloat = std::function<Variable(float)>;
using BinaryFloat = std::function<Variable(float, float)>;

using UnaryVec3 = std::function<Variable(Vector3)>;
using BinaryVec3 = std::function<Variable(Vector3, Vector3)>;

using BinaryVec3Float = std::function<Variable(Vector3, float)>;


class ExpressionParser {
public:
    ExpressionParser();

    std::function<float(const VariableMap&)> parse(const std::string& expression, const VariableMap &variables = {});

    bool validate(const std::string& expression, const VariableMap &variables = {}, bool raiseErrors = true);

    std::set<std::string> extractAllVariables(const std::string& expression);

    // Maps for user-defined operators
    std::map<std::string, std::function<float(float)>> userUnaryOperators;
    std::map<std::string, std::function<float(float, float)>> userBinaryOperators;

protected:
    std::vector<std::string> tokenizeExpression(const std::string& expression);

    Token groupTokensHierarchically(const std::vector<std::string>& tokens, const VariableMap &variables = {});

    std::function<Variable(const VariableMap&)> generateLambda(const Token &token, const VariableMap& variables = {});

    std::vector<std::vector<std::string>> extractObjectPropertyPatterns(const std::vector<std::string>& tokens);

    static VariableMap extendVariables(const VariableMap& _variables);

    // Maps for default operators
    std::map<std::string, BinaryFloat> binaryFloatOperators;
    std::map<std::string, UnaryFloat> unaryFloatOperators;
    std::map<std::string, BinaryVec3> binaryVectorOperators;
    std::map<std::string, UnaryVec3> unaryVectorOperators;
    std::map<std::string, BinaryVec3Float> binaryVectorFloatOperators;

    std::map<std::string, int> precedence;
};

template <typename T, typename U>
std::function<Variable(const VariableMap&)> getBinaryOperation(const std::string& tokenStr,
                                                   const std::function<Variable(const VariableMap&)>& leftLambda,
                                                   const std::function<Variable(const VariableMap&)>& rightLambda,
                                                   const std::map<std::string, std::function<Variable(T, U)>>& operationsMap) {
    auto operationFunc = operationsMap.at(tokenStr);
    return [=](const VariableMap& variables) -> Variable {
        return operationFunc(std::get<T>(leftLambda(variables)), std::get<U>(rightLambda(variables)));
    };
}
template <typename T, typename U>
std::function<Variable(const VariableMap&)> getUnaryOperation(const std::string& tokenStr,
                                                   const std::function<Variable(const VariableMap&)>& lambda,
                                                   const std::map<std::string, std::function<Variable(U)>>& operationsMap) {
    auto operationFunc = operationsMap.at(tokenStr);
    return [=](const VariableMap& variables) -> Variable {
        return operationFunc(std::get<U>(lambda(variables)));
    };
}

#endif // EXPRESSIONPARSER_H
