#pragma once
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include "string"

using std::map;
using std::prev;
using std::vector;

struct Indexes {
    size_t number_row;
    size_t number_column;
};

bool operator<(const Indexes &first_indexes, const Indexes &second_indexes);

class CSRMatrix {
public:
    CSRMatrix(const map<Indexes, double> &matrix, const int &number_rows, const int &number_column);
    double GiveElement(const size_t &number_row, const size_t &number_column) const;
    friend bool operator==(const CSRMatrix &first_matrix, const CSRMatrix &second_matrix);
    vector<double> GiveValues() const;
    vector<size_t> GiveColumnsIndexes() const;
    vector<int> GiveNumberRowsNon0Elements() const;
    int GiveNumberColumns() const;
    int GiveNumberRows() const;

private:
    vector<double> values_;
    vector<size_t> columns_indexes_;
    vector<int> number_rows_non_0_elements_;
    int number_columns_;
    int number_rows_;

    static vector<int> CountNonZeroRowElements(const map<Indexes, double> &matrix);
    double FindElement(const size_t &number_row, const size_t &number_column) const;
};

vector<double> operator*(const CSRMatrix &matrix, const vector<double> &column);


vector<double> operator-(const vector<double> &first_vector, const vector<double> &second_vector);

double operator*(const vector<double> &first_vector, const vector<double> &second_vector);

vector<double> operator*(const double &number, const vector<double> &second_vector);

double GiveVectorLength(const vector<double> &vector);

using std::pair;
using std::string;

pair<vector<double>, int> MPI(const CSRMatrix &a_matrix, const vector<double> &x_0,
                              const vector<double> &b, const double &tolerance, const double &tau,
                              const string &file_path);