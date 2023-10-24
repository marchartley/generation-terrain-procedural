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


struct TokenList;

using Token = std::variant<std::string, std::shared_ptr<TokenList>, Vector3>;
struct TokenList : public std::vector<Token> {
    using std::vector<Token>::vector;
};

using Variable = std::variant<float, Vector3>;

using UnaryFloat = std::function<Variable(float)>;
using BinaryFloat = std::function<Variable(float, float)>;

using UnaryVec3 = std::function<Variable(Vector3)>;
using BinaryVec3 = std::function<Variable(Vector3, Vector3)>;

using BinaryVec3Float = std::function<Variable(Vector3, float)>;


class ExpressionParser {
public:
    ExpressionParser();

    std::function<float(const Vector3&)> parse(const std::string& expression);

    void validate(const std::string& expression);

    // Maps for user-defined operators
    std::map<std::string, std::function<float(float)>> userUnaryOperators;
    std::map<std::string, std::function<float(float, float)>> userBinaryOperators;

protected:
    std::vector<std::string> tokenizeExpression(const std::string& expression);

    Token groupTokensHierarchically(const std::vector<std::string>& tokens);

    std::function<Variable(const Vector3&)> generateLambda(const Token &token);

    std::vector<std::vector<std::string>> extractObjectPropertyPatterns(const std::vector<std::string>& tokens);

    // Maps for default operators
    std::map<std::string, BinaryFloat> binaryFloatOperators;
    std::map<std::string, UnaryFloat> unaryFloatOperators;
    std::map<std::string, BinaryVec3> binaryVectorOperators;
    std::map<std::string, UnaryVec3> unaryVectorOperators;
//    std::map<std::string, BinaryVec3Float> binaryVectorFloatOperators;

    std::map<std::string, int> precedence;
};

template <typename T, typename U>
std::function<Variable(const Vector3&)> getBinaryOperation(const std::string& tokenStr,
                                                   const std::function<Variable(const Vector3&)>& leftLambda,
                                                   const std::function<Variable(const Vector3&)>& rightLambda,
                                                   const std::map<std::string, std::function<Variable(T, U)>>& operationsMap) {
    auto operationFunc = operationsMap.at(tokenStr);
    return [=](const Vector3& p) -> Variable {
        return operationFunc(std::get<T>(leftLambda(p)), std::get<U>(rightLambda(p)));
    };
}
template <typename T, typename U>
std::function<Variable(const Vector3&)> getUnaryOperation(const std::string& tokenStr,
                                                   const std::function<Variable(const Vector3&)>& lambda,
                                                   const std::map<std::string, std::function<Variable(U)>>& operationsMap) {
    auto operationFunc = operationsMap.at(tokenStr);
    return [=](const Vector3& p) -> Variable {
        return operationFunc(std::get<U>(lambda(p)));
    };
}

#endif // EXPRESSIONPARSER_H
