#ifndef MM_IO_H
#define MM_IO_H

#include <algorithm>
#include <iterator>
#include <string>
#include <fstream>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>

// This implementation supposes that the matrices are sparse, both in the constructors of the class and in the
// implemented functions
// The class was written so that functions and further extensions can be easily implemented to handle the dense
// formatting as well

using MM_typecode = std::string;

// Type of mathematical contruct
enum class FormatType {
    matrix
};

// Sparse or dense matrix
enum class SparsityType {
    dense, sparse
};

// Field of the matrix: "Determines the type and number of values listed for each matrix entry" -
//  The Matrix Market	Exchange Format: Initial Design
enum class FieldType {
    real, complex, integer, pattern
};

// Symmetry of the matrix
enum class SymmetryType {
    general, symmetric, hermitian, skewsymmetric
};

class MatrixMarketWriter {
private:

    //File stream specific members
    int outputFile;
    std::string outputFileName;

    //Details of the matrix to be written to file
    double* matrix;
    long rowNumber;
    long columnNumber;
    double cutoff;

    //MarketMatrix format specific information
    FormatType  formatType;
    SparsityType sparsityType;
    FieldType fieldType;
    SymmetryType symmetryType;

    //Writer functions responsible for recording matrix type
    void writeFormat();
    void writeSparsity();
    void writeField();
    void writeSymmetry();
    void writeSparseDetails();
    void writeCoordinates();

    //Writer function that simply dumps the X and Y coordinates and then the value
    void dumpCoordinates();

    //Using POSIX calls to writeToFileBinary binary formatted files
    void writeToFileBinary(std::string outputString);
    void writeToFileBinary(long number);
    void writeToFileBinary(double number);

public:

    // Standard matrices with Row × Column dimensions
    MatrixMarketWriter(std::string fileName, double* mat, long row, long column, FormatType ft = FormatType::matrix,
                       SparsityType st = SparsityType::sparse ) :
            outputFileName(fileName), matrix(mat), rowNumber(row), columnNumber(column),
            formatType (ft), sparsityType(st){}

    // Standard matrices with Row × Column dimensions, considering elements that are greater than the cutoff as nonzero
    MatrixMarketWriter(std::string fileName, double cutoff, double* mat, long row, long column, FormatType ft = FormatType::matrix,
                       SparsityType st = SparsityType::sparse ) :
            outputFileName(fileName), matrix(mat), rowNumber(row), columnNumber(column), cutoff(cutoff),
            formatType (ft), sparsityType(st) {}

    // Square matrices with Dim × Dim dimensions
    MatrixMarketWriter(std::string fileName, double* mat, long dim, FormatType ft = FormatType::matrix,
                       SparsityType st = SparsityType::sparse ) :
            outputFileName(fileName), matrix(mat), rowNumber(dim), columnNumber(dim),
            formatType (ft), sparsityType(st) {}

    // Square matrices with Dim × Dim dimensions, considering elements that are greater than the cutoff as nonzero
    MatrixMarketWriter(std::string fileName, double cutoff, double* mat, long dim, FormatType ft = FormatType::matrix,
                       SparsityType st = SparsityType::sparse ) :
            outputFileName(fileName), matrix(mat), rowNumber(dim), columnNumber(dim), cutoff(cutoff),
            formatType (ft), sparsityType(st) {}

    //Writes the matrix into the files in MatrixMarket format
    bool createBinaryFile();

    //Writes only the matrix coordinate and value sequences into the file
    bool createDumpFile();

/********************* MM_typecode modify functions ***************************/

    inline void setMatrix()   { formatType = FormatType::matrix; }

    inline void setDense()    { sparsityType = SparsityType::dense; }
    inline void setSparse()   { sparsityType = SparsityType::sparse; }

    inline void setReal()    { fieldType = FieldType ::real; }
    inline void setComplex() { fieldType = FieldType ::complex; }
    inline void setInteger() { fieldType = FieldType ::integer; }
    inline void setPattern() { fieldType = FieldType ::pattern; }

    inline void setGeneral()   { symmetryType = SymmetryType::general ; }
    inline void setSymmetric() { symmetryType = SymmetryType::symmetric ; }
    inline void setSkew()      { symmetryType = SymmetryType ::skewsymmetric; }
    inline void setHermitian() { symmetryType = SymmetryType ::hermitian ; }
};

#endif
