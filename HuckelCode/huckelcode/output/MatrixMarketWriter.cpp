#include <cmath>
#include "MatrixMarketWriter.h"

void MatrixMarketWriter::writeToFileBinary(std::string outputString)
{
    write(outputFile, outputString.c_str(),outputString.size());
}

void MatrixMarketWriter::writeToFileBinary(double number)
{
    write(outputFile, (const char *) &number, sizeof(double));
}

void MatrixMarketWriter::writeToFileBinary(long number)
{
    write(outputFile, (const char *) &number, sizeof(long));
}

//Writing the matrix to file in binary mode
bool MatrixMarketWriter::createBinaryFile()
{
    //Open file using POSIX system call for binary formatted writing
    outputFile = open(outputFileName.c_str(), O_RDWR | O_TRUNC | O_APPEND | O_CREAT, S_IRUSR | S_IWUSR);

    //If we managed to open the file, write the necessary data
    if (outputFile != -1)
    {
        //Writing the first line of the file containing the type of matrix
        writeToFileBinary("%%MatrixMarket");
        writeFormat();
        writeSparsity();
        writeField();
        writeSymmetry();
        writeToFileBinary("\n");

        //COMMENTS for the Matrix Market file COME HERE

        //Writing the dimensions and the number of nonzero elements
        writeSparseDetails();

        //Writing the coordinates of the matrix
        writeCoordinates();
    }
    else
    {
        std::cerr << "Cannot open output file for binary writing, please check if you have rights to write the file." << std::endl;
        return EXIT_FAILURE;
    }
    close(outputFile);
    return EXIT_SUCCESS;
}

bool MatrixMarketWriter::createDumpFile()
{
    //Open file using POSIX system call for binary formatted writing
    outputFile = open(outputFileName.c_str(), O_RDWR | O_TRUNC | O_APPEND | O_CREAT, S_IRUSR | S_IWUSR);

    //If we managed to open the file, write the necessary data
    if (outputFile != -1)
    {
        dumpCoordinates();
    }
    else
    {
        std::cerr << "Cannot open output file for binary writing, please check if you have rights to write the file." << std::endl;
        return EXIT_FAILURE;
    }
    close(outputFile);
    return EXIT_SUCCESS;
}

void MatrixMarketWriter::writeFormat()
{
    switch(formatType) {
        case FormatType::matrix:
            writeToFileBinary(" matrix");
            break;
    }
}

void MatrixMarketWriter::writeSparsity()
{
    switch(sparsityType) {
        case SparsityType::sparse:
            writeToFileBinary(" coordinate");
            break;
        case SparsityType::dense:
            writeToFileBinary(" array");
            break;
    }
}

void MatrixMarketWriter::writeField()
{
    switch(fieldType) {
        case FieldType::real:
            writeToFileBinary(" real");
            break;
        case FieldType::complex:
            writeToFileBinary(" complex");
            break;
        case FieldType::integer:
            writeToFileBinary(" integer");
            break;
        case FieldType::pattern:
            writeToFileBinary(" pattern");
            break;
    }
}

void MatrixMarketWriter::writeSymmetry()
{
    switch(symmetryType) {
        case SymmetryType::general:
            writeToFileBinary(" general");
            break;
        case SymmetryType::hermitian:
            writeToFileBinary(" hermitian");
            break;
        case SymmetryType::symmetric:
            writeToFileBinary(" symmetric");
            break;
        case SymmetryType::skewsymmetric:
            writeToFileBinary(" skew-symmetric");
            break;
    }
}

//Counting number of non-zero elements of the matrix using multiple threads if possible
void MatrixMarketWriter::writeSparseDetails()
{
    long numberOfNonZero = 0;
    #pragma omp parallel for shared(numberOfNonZero)
    for (long i = 0; i < rowNumber; i++) {
        long itab = i * columnNumber;
        for (long j = 0; j < columnNumber; j++){
            if (fabs(matrix[itab + j]) > cutoff)
            {
                #pragma omp atomic
                ++numberOfNonZero;
            }
        }
    }
    writeToFileBinary(rowNumber);
    writeToFileBinary(" ");
    writeToFileBinary(columnNumber);
    writeToFileBinary(" ");
    writeToFileBinary(numberOfNonZero);
    writeToFileBinary("\n");
}

void MatrixMarketWriter::writeCoordinates()
{
    // write the nonzero elements
    for (long i = 0; i < rowNumber; i++) {
        long itab = i * columnNumber;
        for (long j = 0; j < columnNumber; j++){
            if (fabs(matrix[itab + j]) > cutoff)
            {
                writeToFileBinary(i);
                writeToFileBinary(" ");
                writeToFileBinary(j);
                writeToFileBinary(" ");
                writeToFileBinary(matrix[itab + j]);
                writeToFileBinary("\n");
            }
        }
    }
}

void MatrixMarketWriter::dumpCoordinates()
{
    // write the nonzero elements
    for (long i = 0; i < rowNumber; i++) {
        long itab = i * columnNumber;
        for (long j = 0; j < columnNumber; j++){
            if (fabs(matrix[itab + j]) > cutoff)
            {
                writeToFileBinary(i);
                writeToFileBinary(j);
                writeToFileBinary(matrix[itab + j]);
            }
        }
    }
}