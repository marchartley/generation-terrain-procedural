#ifndef TABLE_H
#define TABLE_H

#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <variant>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <numeric>

class Table {
public:
    using DataVariant = std::variant<std::string, float>;

    template<typename... Ts>
    struct TableRow {
        std::tuple<Ts...> data;
        TableRow(Ts... values) : data(values...) {}
    };

    Table(const std::vector<std::vector<float>>& data,
          const std::vector<std::string>& colNames,
          const std::vector<std::string>& rowNames);

    template<typename... Ts>
    Table(const std::vector<TableRow<Ts...>>& tableRows,
          const std::vector<std::string>& colNames,
          const std::vector<std::string>& rowNames)
        : colNames(colNames), rowNames(rowNames) {
        for (const auto& row : tableRows) {
            std::vector<DataVariant> convertedRow = { std::get<Ts>(row.data)... };
            rows.push_back(convertedRow);
        }
    }

    void addRow(const std::vector<DataVariant>& rowData, const std::string& rowName);

    template<typename T>
    void addColumn(const std::vector<T>& newColumn, const std::string& colName) {
        if (newColumn.size() != rows.size()) {
            throw std::runtime_error("The new column size does not match the number of rows in the table.");
        }

        for (size_t i = 0; i < rows.size(); ++i) {
            rows[i].push_back(newColumn[i]);
        }

        colNames.push_back(colName);
    }

    Table& sortBy(const std::string& columnName, bool reversed = false);

    std::string displayTable() const;

    std::string toCSV() const;

private:
    std::string variantToStr(const DataVariant& var) const;

    std::vector<std::vector<DataVariant>> rows;
    std::vector<std::string> colNames;
    std::vector<std::string> rowNames;
};

#endif // TABLE_H
