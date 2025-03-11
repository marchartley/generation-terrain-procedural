#include "Table.h"

Table::Table(const std::vector<std::vector<float> > &data, const std::vector<std::string> &colNames, const std::vector<std::string> &rowNames)
    : colNames(colNames), rowNames(rowNames) {
    for (const auto& row : data) {
        std::vector<DataVariant> convertedRow;
        for (const auto& val : row) {
            convertedRow.push_back(val);
        }
        rows.push_back(convertedRow);
    }
}

Table::Table(const std::vector<std::vector<std::string> > &data, const std::vector<std::string> &colNames, const std::vector<std::string> &rowNames)
    : colNames(colNames), rowNames(rowNames) {
    for (const auto& row : data) {
        std::vector<DataVariant> convertedRow;
        for (const auto& val : row) {
            convertedRow.push_back(val);
        }
        rows.push_back(convertedRow);
    }
}

void Table::addRow(const std::vector<DataVariant> &rowData, const std::string &rowName) {
    if(rowData.size() != colNames.size()) {
        throw std::invalid_argument("Mismatched number of columns in the provided row data");
    }

    rows.push_back(rowData);
    rowNames.push_back(rowName);
}

Table& Table::sortBy(const std::string &columnName, bool reversed) {
    auto colIt = std::find(colNames.begin(), colNames.end(), columnName);
    if (colIt == colNames.end()) {
        throw std::invalid_argument("Column name not found");
    }

    size_t colIndex = std::distance(colNames.begin(), colIt);

    auto comparator = [&](const std::vector<DataVariant>& a, const std::vector<DataVariant>& b) -> bool {
        if (std::holds_alternative<float>(a[colIndex]) && std::holds_alternative<float>(b[colIndex])) {
            if (reversed) {
                return std::get<float>(a[colIndex]) > std::get<float>(b[colIndex]);
            } else {
                return std::get<float>(a[colIndex]) < std::get<float>(b[colIndex]);
            }
        } else if (std::holds_alternative<std::string>(a[colIndex]) && std::holds_alternative<std::string>(b[colIndex])) {
            if (reversed) {
                return std::get<std::string>(a[colIndex]) > std::get<std::string>(b[colIndex]);
            } else {
                return std::get<std::string>(a[colIndex]) < std::get<std::string>(b[colIndex]);
            }
        }
        throw std::runtime_error("Unsupported data type or mismatched column types");
    };

    std::vector<size_t> permutation(rowNames.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(), [&](size_t i, size_t j) { return comparator(rows[i], rows[j]); });

    std::vector<std::vector<DataVariant>> sortedRows(rows.size());
    std::vector<std::string> sortedRowNames(rowNames.size());
    for (size_t i = 0; i < permutation.size(); ++i) {
        sortedRows[i] = rows[permutation[i]];
        sortedRowNames[i] = rowNames[permutation[i]];
    }

    rows = std::move(sortedRows);
    rowNames = std::move(sortedRowNames);
    return *this;
}

std::string Table::displayTable() const {
    std::stringstream ss;
    if (rows.size() != rowNames.size() || (rows.size() > 0 && rows[0].size() != colNames.size())) {
        return "Invalid data or column/row names!";
    }

    std::vector<size_t> colWidths(colNames.size(), 0);
    for (size_t j = 0; j < colNames.size(); ++j) {
        colWidths[j] = colNames[j].size();
        for (size_t i = 0; i < rows.size(); ++i) {
            colWidths[j] = std::max(colWidths[j], variantToStr(rows[i][j]).size());
        }
    }

    ss << std::setw(8) << " ";
    for (size_t j = 0; j < colNames.size(); ++j) {
        ss << std::setw(colWidths[j] + 2) << colNames[j];
    }
    ss << '\n';

    for (size_t i = 0; i < rows.size(); ++i) {
        ss << std::setw(8) << rowNames[i];
        for (size_t j = 0; j < rows[i].size(); ++j) {
            ss << std::setw(colWidths[j] + 2) << variantToStr(rows[i][j]);
        }
        ss << '\n';
    }

    return ss.str();
}

std::string Table::toCSV() const {
    std::stringstream ss;
    ss << ",";  // Empty cell for the top-left corner
    for (const auto& col : colNames) {
        ss << col << ",";
    }
    ss << '\n';

    for (size_t i = 0; i < rows.size(); ++i) {
        ss << rowNames[i] << ",";
        for (const auto& cell : rows[i]) {
            ss << variantToStr(cell) << ",";
        }
        ss << '\n';
    }

    return ss.str();
}

std::string Table::variantToStr(const DataVariant &var) const {
    if (std::holds_alternative<float>(var)) {
        return std::to_string(std::get<float>(var));
    } else {
        return std::get<std::string>(var);
    }
}
