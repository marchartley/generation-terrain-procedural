#include "ExpressionParser.h"
#include <stack>
#include "Utils/Utils.h"

ExpressionParser::ExpressionParser() {
    // Default binary operators
    binaryFloatOperators["+"] = [](float a, float b) { return a + b; };
    binaryFloatOperators["-"] = [](float a, float b) { return a - b; };
    binaryFloatOperators["*"] = [](float a, float b) { return a * b; };
    binaryFloatOperators["/"] = [](float a, float b) { return a / b; };
    binaryFloatOperators["^"] = [](float a, float b) { return std::pow(a, b); };
    binaryFloatOperators["max"] = [](float a, float b) { return std::max(a, b); };
    binaryFloatOperators["min"] = [](float a, float b) { return std::min(a, b); };
    binaryFloatOperators["%"] = [](float a, float b) { return std::fmod(a, b); };

    // Default unary operators
    unaryFloatOperators["abs"] = [](float a) { return std::abs(a); };
    unaryFloatOperators["sqrt"] = [](float a) { return std::sqrt(a); };
    unaryFloatOperators["pow2"] = [](float a) { return a * a; };
    unaryFloatOperators["normal"] = [](float a) -> float { return 2.506628275f * std::pow(M_E, -(a*a)*.5f); }; // standard normal:  1/sqrt(2 pi) * e^[-x^2/2]
    // For any other normal(x, mu, sigma), consider calling "normal((x - mu)/sigma)"

    // Comparator operators
    binaryFloatOperators[">"] = [](float a, float b) { return (a > b ? 1.f : 0.f); };
    binaryFloatOperators["<"] = [](float a, float b) { return (a < b ? 1.f : 0.f); };
    binaryFloatOperators["="] = [](float a, float b) { return (std::abs(a - b) < 1e-6 ? 1.f : 0.f); }; // Meh... epsilon?

    // Default binary operators
    binaryVectorOperators["+"] = [](const Vector3& a, const Vector3& b) { return a + b; };
    binaryVectorOperators["-"] = [](const Vector3& a, const Vector3& b) { return a - b; };
    binaryVectorOperators["*"] = [](const Vector3& a, const Vector3& b) { return a * b; };
    binaryVectorOperators["/"] = [](const Vector3& a, const Vector3& b) { return a / b; };
    binaryVectorOperators["dot"] = [](const Vector3& a, const Vector3& b) { return a.dot(b); };
    binaryVectorOperators["cross"] = [](const Vector3& a, const Vector3& b) { return a.cross(b); };
    binaryVectorOperators["?"] = [](const Vector3& a, const Vector3& b) { return (a.isValid() ? a : b); };

    // Default unary operators
    unaryVectorOperators["abs"] = [](const Vector3& a) { return std::abs(a); };
    unaryVectorOperators["norm2"] = [](const Vector3& a) { return a.norm2(); };
    unaryVectorOperators["d2"] = [](const Vector3& a) { return a.norm2(); };
    unaryVectorOperators["normalize"] = [](const Vector3& a) { return a.normalized(); };
    unaryVectorOperators["d"] = [](const Vector3& a) { return a.norm(); };
    unaryVectorOperators["random"] = [](const Vector3& a) { return random_gen::generate_perlin(a.x, a.y, a.z); };
    unaryVectorOperators["x"] = [](const Vector3& a) { return a.x; };
    unaryVectorOperators["y"] = [](const Vector3& a) { return a.y; };
    unaryVectorOperators["z"] = [](const Vector3& a) { return a.z; };

    binaryVectorFloatOperators["*"] = [](const Vector3& a, float b) { return a * b; };
    binaryVectorFloatOperators["/"] = [](const Vector3& a, float b) { return a / b; };

    precedence = {
            {"+", 1},
            {"-", 1},
            {"*", 2},
            {"/", 2},
            {"^", 3}
    };

    for (const auto& [symbol, func] : binaryFloatOperators) {
        if (precedence.count(symbol) == 0) {
            precedence[symbol] = 4;
        }
    }
    for (const auto& [symbol, func] : binaryVectorOperators) {
        if (precedence.count(symbol) == 0) {
            precedence[symbol] = 4;
        }
    }
}

std::function<float (const VariableMap &)> ExpressionParser::parse(const std::string &expression, const VariableMap& _variables) {
    VariableMap variables = extendVariables(_variables);
    auto tokens = tokenizeExpression(expression);
    /*for (auto& t : tokens) {
        std::cout << t << " ";
    }
    std::cout << std::endl;*/
    auto groupedTokens = groupTokensHierarchically(tokens, variables);
    auto returnedLambda = generateLambda(groupedTokens, variables);

    // Create a dummy Vector3
//    Vector3 dummyVec;

    // Invoke the lambda with the dummy Vector3
    std::variant<float, Vector3> result = returnedLambda(variables);

    // Check the type of the result
    if (!std::holds_alternative<float>(result)) {
        throw std::runtime_error("The expression does not return a float value.");
    }

    // Convert the lambda to the desired type
    return [returnedLambda](const VariableMap& vars) -> float {
        return std::get<float>(returnedLambda(ExpressionParser::extendVariables(vars)));
    };
}

bool ExpressionParser::validate(const std::string &expression, const VariableMap& _variables, bool raiseErrors) {

    VariableMap variables = extendVariables(_variables);
    std::vector<std::string> allErrors;
    // Tokenize the expression
    auto tokens = tokenizeExpression(expression);

    // Check for matching brackets
    std::stack<char> bracketStack;
    bracketStack.push('-'); // Limit! (used to avoid hitting out of memory when checking all brackets)
    for (const auto& token : tokens) {
        if (bracketStack.empty()) break; // There has been an error, at least '-' should have been kept
        if (token.size() == 1) {
            char ch = token[0];
            switch (ch) {
                case '(':
                case '{':
                case '[':
//                case '<':
                    bracketStack.push(ch);
                    break;
                case ')':
                    if (bracketStack.empty() || bracketStack.top() != '(') {
                        allErrors.push_back("Mismatched closing parenthesis.");
                    }
                    bracketStack.pop();
                    break;
                case '}':
                    if (bracketStack.empty() || bracketStack.top() != '{') {
                        allErrors.push_back("Mismatched closing curly brace.");
                    }
                    bracketStack.pop();
                    break;
                case ']':
                    if (bracketStack.empty() || bracketStack.top() != '[') {
                        allErrors.push_back("Mismatched closing square bracket.");
                    }
                    bracketStack.pop();
                    break;
                /*case '>':
                    if (bracketStack.empty() || bracketStack.top() != '<') {
                        allErrors.push_back("Mismatched closing angle bracket.");
                    }
                    bracketStack.pop();
                    break;*/
            }
        }
    }
//    if (!bracketStack.empty()) {
    if (bracketStack.size() != 1) { // Need to find the '-' here
        allErrors.push_back("Mismatched opening bracket.");
//        return false;
    }

    // Check that all variables have a type
    for (const auto& token : tokens) {
        if (precedence.count(token) > 0 || unaryFloatOperators.count(token) > 0 || unaryVectorOperators.count(token) > 0) continue;
        if (std::isalpha(token[0])) { // Assuming variables start with an alphabet character
            if (variables.find(token) == variables.end()) {
                allErrors.push_back("Variable '" + token + "' does not have a defined type.");
            }
        }
    }

    if (raiseErrors) {
        auto groupedTokens = groupTokensHierarchically(tokens, variables);
        VariableMap dummyVars;
        for (const auto& varTypePair : variables) {
            if (std::holds_alternative<float>(varTypePair.second)) {
                dummyVars[varTypePair.first] = 0.0f;
            } else if (std::holds_alternative<Vector3>(varTypePair.second)) {
                dummyVars[varTypePair.first] = Vector3{0.0f, 0.0f, 0.0f};
            }
        }
        auto lambda = generateLambda(groupedTokens, dummyVars);
        lambda(dummyVars); // Optionally invoke the lambda
    } else {

        // Attempt to generate the lambda and catch any exceptions
        try {
            auto groupedTokens = groupTokensHierarchically(tokens, variables);
            VariableMap dummyVars;
            for (const auto& varTypePair : variables) {
                if (std::holds_alternative<float>(varTypePair.second)) {
                    dummyVars[varTypePair.first] = 0.0f;
                } else if (std::holds_alternative<Vector3>(varTypePair.second)) {
                    dummyVars[varTypePair.first] = Vector3{0.0f, 0.0f, 0.0f};
                }
            }
            auto lambda = generateLambda(groupedTokens, dummyVars);
            lambda(dummyVars); // Optionally invoke the lambda
        } catch (const std::exception& e) {
            allErrors.push_back("Error in the evaluation of the function: " + std::string(e.what()));
        }
        if (!allErrors.empty()) {
            allErrors.insert(allErrors.begin(), "Checking " + expression);
//            if (raiseErrors)
//                throw std::runtime_error(join(allErrors, "\n"));
//            else
            std::cerr << join(allErrors, "\n") << std::endl;
            return false;
        }
    }
    return true;
}

std::set<std::string> ExpressionParser::extractAllVariables(const std::string &expression)
{
    std::set<std::string> variablesUsedInExpression;
    auto tokens = tokenizeExpression(expression);
    for (const auto& token : tokens) {
        if (precedence.count(token) > 0 || unaryFloatOperators.count(token) > 0 || unaryVectorOperators.count(token) > 0) continue;
        if (std::isalpha(token[0])) { // Assuming variables start with an alphabet character
            variablesUsedInExpression.insert(token);
        }
    }
    return variablesUsedInExpression;
}

std::vector<std::string> ExpressionParser::tokenizeExpression(const std::string &expression) {
    std::vector<std::string> tokens;
    std::string currentToken;

    for (size_t i = 0; i < expression.size(); i++) {
        char ch = expression[i];
        if (std::isspace(ch) || ch == ',') {
            if (!currentToken.empty()) {
                tokens.push_back(currentToken);
                currentToken.clear();
            }
        } else if (std::isalnum(ch) || ch == '.' || (ch == '-' && std::isalnum(expression[i + 1]))) {
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


Token ExpressionParser::groupTokensHierarchically(const std::vector<std::string> &tokens, const VariableMap& variables) {
    TokenList groupedTokens;
    std::string currentToken;

    for (size_t i = 0; i < tokens.size(); ++i) {
        const auto& token = tokens[i];
        if (token == "(" /*|| token == "<"*/ || token == "{") {
            size_t depth = 1;
            size_t j = i + 1;
            for (; j < tokens.size() && depth > 0; ++j) {
                if (tokens[j] == "(" /*|| tokens[j] == "<"*/ || tokens[j] == "{") {
                    depth++;
                } else if (tokens[j] == ")" /*|| tokens[j] == ">"*/ || tokens[j] == "}") {
                    depth--;
                }
            }
            if (token == "{") {
                // Convert tokens inside {} to Vector3
                auto components = *std::get<std::shared_ptr<TokenList>>(groupTokensHierarchically(std::vector<std::string>(tokens.begin() + i + 1, tokens.begin() + j - 1), variables));
                if (components.size() == 3) {
                    float x = std::get<float>(generateLambda(components[0], variables)(variables));
                    float y = std::get<float>(generateLambda(components[1], variables)(variables));
                    float z = std::get<float>(generateLambda(components[2], variables)(variables));

                    Vector3 vec(x, y, z);
                    groupedTokens.push_back(vec);
                } else {
                    // Handle error: Invalid Vector3 format
                    //std::cerr << "Error..." << std::endl;
                    throw std::runtime_error("Could not parse the vector ...");
                }
            } else {
                auto innerTokens = groupTokensHierarchically(std::vector<std::string>(tokens.begin() + i + 1, tokens.begin() + j - 1), variables);
                groupedTokens.push_back(innerTokens);
            }
            i = j - 1;  // Skip processed tokens
        } else {
            if (variables.count(token) > 0 && false) {
                if (std::holds_alternative<float>(variables.at(token)))
                    groupedTokens.push_back(std::get<Vector3>(variables.at(token))); // Transform the variable to float if possible
                if (std::holds_alternative<Vector3>(variables.at(token)))
                    groupedTokens.push_back(std::get<Vector3>(variables.at(token))); // Transform the variable to vector3 if possible
            } else {
                groupedTokens.push_back(token); // Otherwise, add the original token
            }
        }
    }

    return std::make_shared<TokenList>(groupedTokens);
}

bool parse_float(std::string in, double& res) {
    try {
        size_t read= 0;
        res = std::stod(in, &read);
        if (in.size() != read) {
            std::replace(in.begin(), in.end(), '.', ',');
            res = std::stod(in, &read);
            if (in.size() != read)
                return false;
        }
    } catch (std::invalid_argument) {
        return false;
    }
    return true;
}

std::function<Variable (const VariableMap &)> ExpressionParser::generateLambda(const Token &token, const VariableMap &variables) {
    if (std::holds_alternative<std::string>(token)) {
        std::string tokenStr = std::get<std::string>(token);
        double dVal;
        if (parse_float(tokenStr, dVal)) {
            float fVal = dVal;
            return [=](const VariableMap&) { return fVal; };
        } else {
            return [=](const VariableMap& vars) {
                if (vars.count(tokenStr)) return vars.at(tokenStr);
                else return Variable(); //0.f;
            };
        }
        /*try {
            float value = (tokenStr == "" ? 0.f : std::atof(tokenStr.c_str()));
            return [=](const VariableMap&) { return value; };
        } catch(const std::invalid_argument& e) {
            return [=](const VariableMap& vars) {
                if (vars.count(tokenStr)) return vars.at(tokenStr);
                else return Variable(); //0.f;
            };
        }*/
    } else if (std::holds_alternative<Vector3>(token)) {
        Vector3 value = std::get<Vector3>(token);
        return [=](const VariableMap&) { return value; };
    }

    const TokenList& tokens = *std::get<std::shared_ptr<TokenList>>(token);

    if (tokens.size() == 1)
        return generateLambda(tokens[0], variables);

    // Find the operator with the lowest precedence
    int min_precedence = INT_MAX;
    size_t op_index = 0;
    for (size_t i = 1; i < tokens.size(); ++i) {
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
        auto leftLambda = generateLambda(std::make_shared<TokenList>(TokenList(tokens.begin(), tokens.begin() + op_index)), variables);
        auto rightLambda = generateLambda(std::make_shared<TokenList>(TokenList(tokens.begin() + op_index + 1, tokens.end())), variables);

        auto leftResult = leftLambda(variables);
        auto rightResult = rightLambda(variables);

        if (std::holds_alternative<float>(leftResult) && std::holds_alternative<float>(rightResult) && binaryFloatOperators.find(token) != binaryFloatOperators.end()) {
//                auto opFunction = binaryFloatOperators[token];
            return getBinaryOperation<float>(token, leftLambda, rightLambda, binaryFloatOperators);
        } else if (std::holds_alternative<Vector3>(leftResult) && std::holds_alternative<Vector3>(rightResult) && binaryVectorOperators.find(token) != binaryVectorOperators.end()) {
            return getBinaryOperation<Vector3>(token, leftLambda, rightLambda, binaryVectorOperators);
        } else if (std::holds_alternative<Vector3>(leftResult) && std::holds_alternative<float>(rightResult) && binaryVectorFloatOperators.find(token) != binaryVectorFloatOperators.end()) {
            return getBinaryOperation<Vector3>(token, leftLambda, rightLambda, binaryVectorFloatOperators);
        } else {
            // std::cerr << "No operation found for this equation..." << std::endl;
            throw std::runtime_error("No operation found for this equation...");
            /*return [=](const VariableMap&) -> std::function<Variable(const VariableMap&)> {
                return generateLambda("0.0"); // [=](const VariableMap&) -> Variable { return 0.f; };
            };*/
        }
    }

    int nbTokens = int(tokens.size());
    int maxI = nbTokens - 1;
    for (int i = 0; i < maxI; ++i) {
        if (std::holds_alternative<std::string>(tokens[i])) {
            const auto& token = std::get<std::string>(tokens[i]);
            auto operandLambda = generateLambda(tokens[i + 1], variables);

            Vector3 dummyVec;
            auto lambdaResult = operandLambda(variables);

            if (std::holds_alternative<float>(lambdaResult) && unaryFloatOperators.find(token) != unaryFloatOperators.end()) {
                return getUnaryOperation<float>(token, operandLambda, unaryFloatOperators);
            } else if (std::holds_alternative<Vector3>(lambdaResult) && unaryVectorOperators.find(token) != unaryVectorOperators.end()) {
                return getUnaryOperation<Vector3>(token, operandLambda, unaryVectorOperators);
            }
        }
    }
    // If no recognized patterns, return a lambda that returns 0 (or handle it appropriately one day...)
//    std::cerr << "Careful! An operation hasn't be handled correctly and is transformed as 0.0" << std::endl;
//    return [](const VariableMap&) { return 0.0f; };
    throw std::runtime_error("Careful! An operation hasn't be handled correctly and is transformed as 0.0");
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

VariableMap ExpressionParser::extendVariables(const VariableMap &_variables)
{
    VariableMap variables = _variables;
    variables["e"] = float(M_E);
    variables["pi"] = float(M_PI);
    return variables;
    /*
    for (auto& [name, _var] : _variables) {
        if (std::holds_alternative<Vector3>(_var)) {
            Vector3 var = std::get<Vector3>(_var);
            variables[name + ".x"] = var.x;
            variables[name + ".y"] = var.y;
            variables[name + ".z"] = var.z;
        }
    }
    return variables;
    */
}
