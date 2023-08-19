#include "ExpressionParser.h"


ExpressionParser::ExpressionParser() {
    // Default binary operators
    binaryFloatOperators["+"] = [](float a, float b) { return a + b; };
    binaryFloatOperators["-"] = [](float a, float b) { return a - b; };
    binaryFloatOperators["*"] = [](float a, float b) { return a * b; };
    binaryFloatOperators["/"] = [](float a, float b) { return a / b; };
    binaryFloatOperators["^"] = [](float a, float b) { return std::pow(a, b); };

    // Default unary operators
    unaryFloatOperators["abs"] = [](float a) { return std::abs(a); };
    unaryFloatOperators["sqrt"] = [](float a) { return std::sqrt(a); };
    unaryFloatOperators["pow2"] = [](float a) { return a * a; };

    // Default binary operators
    binaryVectorOperators["+"] = [](const Vector3& a, const Vector3& b) { return a + b; };
    binaryVectorOperators["-"] = [](const Vector3& a, const Vector3& b) { return a - b; };
    binaryVectorOperators["*"] = [](const Vector3& a, const Vector3& b) { return a * b; };
    binaryVectorOperators["/"] = [](const Vector3& a, const Vector3& b) { return a / b; };

    // Default unary operators
    unaryVectorOperators["abs"] = [](const Vector3& a) { return std::abs(a); };
    unaryVectorOperators["norm2"] = [](const Vector3& a) { return a.norm2(); };
    unaryVectorOperators["norm"] = [](const Vector3& a) { return a.norm(); };

    precedence = {
            {"+", 1},
            {"-", 1},
            {"*", 2},
            {"/", 2},
            {"^", 3}
    };
}

std::function<float (const Vector3 &)> ExpressionParser::parse(const std::string &expression) {
    auto tokens = tokenizeExpression(expression);
    for (auto& t : tokens) {
        std::cout << t << " ";
    }
    std::cout << std::endl;
    auto groupedTokens = groupTokensHierarchically(tokens);
    auto returnedLambda = generateLambda(groupedTokens);

    // Create a dummy Vector3
    Vector3 dummyVec;

    // Invoke the lambda with the dummy Vector3
    std::variant<float, Vector3> result = returnedLambda(dummyVec);

    // Check the type of the result
    if (!std::holds_alternative<float>(result)) {
        throw std::runtime_error("The expression does not return a float value.");
    }

    // Convert the lambda to the desired type
    return [returnedLambda](const Vector3& vec) -> float {
        return std::get<float>(returnedLambda(vec));
    };
}

void ExpressionParser::validate(const std::string &expression) {
    // Empty for now
}

std::vector<std::string> ExpressionParser::tokenizeExpression(const std::string &expression) {
    std::vector<std::string> tokens;
    std::string currentToken;

    for (char ch : expression) {
        if (std::isspace(ch) || ch == ',') {
            if (!currentToken.empty()) {
                tokens.push_back(currentToken);
                currentToken.clear();
            }
        } else if (std::isalnum(ch) || ch == '.') {
            currentToken += ch;
        } else {
            if (!currentToken.empty()) {
                tokens.push_back(currentToken);
                currentToken.clear();
            }
            tokens.push_back(std::string(1, ch));
        }
    }

    if (!currentToken.empty()) {
        tokens.push_back(currentToken);
    }

    return tokens;
}


Token ExpressionParser::groupTokensHierarchically(const std::vector<std::string> &tokens) {
    TokenList groupedTokens;
    std::string currentToken;

    for (size_t i = 0; i < tokens.size(); ++i) {
        const auto& token = tokens[i];
        if (token == "(" || token == "<" || token == "{") {
            size_t depth = 1;
            size_t j = i + 1;
            for (; j < tokens.size() && depth > 0; ++j) {
                if (tokens[j] == "(" || tokens[j] == "<" || tokens[j] == "{") {
                    depth++;
                } else if (tokens[j] == ")" || tokens[j] == ">" || tokens[j] == "}") {
                    depth--;
                }
            }
            if (token == "{") {
                // Convert tokens inside {} to Vector3
                auto components = *std::get<std::shared_ptr<TokenList>>(groupTokensHierarchically(std::vector<std::string>(tokens.begin() + i + 1, tokens.begin() + j - 1)));
                if (components.size() == 3) {
                    auto _x = components[0];
                    auto _y = components[1];
                    auto _z = components[2];
                    std::string sx = std::get<std::string>(_x);
                    std::string sy = std::get<std::string>(_y);
                    std::string sz = std::get<std::string>(_z);
                    float x = std::stof(sx);
                    float y = std::stof(sy);
                    float z = std::stof(sz);
                    Vector3 vec(x, y, z);
                    groupedTokens.push_back(vec);
                } else {
                    // Handle error: Invalid Vector3 format
                    std::cerr << "Error..." << std::endl;
                }
            } else {
                auto innerTokens = groupTokensHierarchically(std::vector<std::string>(tokens.begin() + i + 1, tokens.begin() + j - 1));
                groupedTokens.push_back(innerTokens);
            }
            i = j - 1;  // Skip processed tokens
        } else {
            groupedTokens.push_back(token);
        }
    }

    return std::make_shared<TokenList>(groupedTokens);
}

std::function<Variable (const Vector3 &)> ExpressionParser::generateLambda(const Token &token) {
    if (std::holds_alternative<std::string>(token)) {
        float value = std::stof(std::get<std::string>(token));
        return [=](const Vector3&) { return value; };
    } else if (std::holds_alternative<Vector3>(token)) {
        Vector3 value = std::get<Vector3>(token);
        return [=](const Vector3&) { return value; };
    }

    const TokenList& tokens = *std::get<std::shared_ptr<TokenList>>(token);

    if (tokens.size() == 1)
        return generateLambda(tokens[0]);

    // Find the operator with the lowest precedence
    int min_precedence = INT_MAX;
    size_t op_index = 0;
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (std::holds_alternative<std::string>(tokens[i])) {
            const auto& tok = std::get<std::string>(tokens[i]);
            if (precedence.find(tok) != precedence.end() && precedence[tok] <= min_precedence) {
                min_precedence = precedence[tok];
                op_index = i;
            }
        }
    }

    // If an operator was found, process it
    if (min_precedence != INT_MAX) {
        const auto& token = std::get<std::string>(tokens[op_index]);
        auto leftLambda = generateLambda(std::make_shared<TokenList>(TokenList(tokens.begin(), tokens.begin() + op_index)));
        auto rightLambda = generateLambda(std::make_shared<TokenList>(TokenList(tokens.begin() + op_index + 1, tokens.end())));

        // Use dummy Vector3 to determine the return type of the lambdas
        Vector3 dummyVec;
        auto leftResult = leftLambda(dummyVec);
        auto rightResult = rightLambda(dummyVec);

        if (std::holds_alternative<float>(leftResult) && std::holds_alternative<float>(rightResult) && binaryFloatOperators.find(token) != binaryFloatOperators.end()) {
            return getBinaryOperation<float>(token, leftLambda, rightLambda, binaryFloatOperators);
        } else if (std::holds_alternative<Vector3>(leftResult) && std::holds_alternative<Vector3>(rightResult) && binaryVectorOperators.find(token) != binaryVectorOperators.end()) {
            return getBinaryOperation<Vector3>(token, leftLambda, rightLambda, binaryVectorOperators);
        }
//        auto binaryFunc = binaryFloatOperators[token];
//        return [=](const Vector3& p) -> Return { return binaryFunc(leftLambda(p), rightLambda(p)); };
    }

    for (size_t i = 0; i < tokens.size(); ++i) {
        if (std::holds_alternative<std::string>(tokens[i])) {
            const auto& token = std::get<std::string>(tokens[i]);

            // Unary functions
//            if (unaryFloatOperators.find(token) != unaryFloatOperators.end() && i + 1 < tokens.size()) {
                auto operandLambda = generateLambda(tokens[i + 1]);

                Vector3 dummyVec;
                auto lambdaResult = operandLambda(dummyVec);

                if (std::holds_alternative<float>(lambdaResult) && unaryFloatOperators.find(token) != unaryFloatOperators.end()) {
                    return getUnaryOperation<float>(token, operandLambda, unaryFloatOperators);
                } else if (std::holds_alternative<Vector3>(lambdaResult) && unaryVectorOperators.find(token) != unaryVectorOperators.end()) {
                    return getUnaryOperation<Vector3>(token, operandLambda, unaryVectorOperators);
                }
//                auto unaryFunc = unaryFloatOperators[token];
//                return [=](const Vector3& p) -> Variable { return unaryFunc(operandLambda(p)); };
//            }
        }
    }
    // If no recognized patterns, return a lambda that returns 0 (or handle appropriately)
    return [](const Vector3&) { return 0.0f; };
}

std::vector<std::vector<std::string> > ExpressionParser::extractObjectPropertyPatterns(const std::vector<std::string> &tokens) {
    std::vector<std::vector<std::string>> results;

    for (size_t i = 0; i < tokens.size() - 2; ++i) {
        if (tokens[i + 1] == "." && !std::isdigit(tokens[i][0]) && !std::isdigit(tokens[i + 2][0])) {
            results.push_back({tokens[i], tokens[i + 1], tokens[i + 2]});
            i += 2;  // Skip the next two tokens since we've processed them
        }
    }

    return results;
}
