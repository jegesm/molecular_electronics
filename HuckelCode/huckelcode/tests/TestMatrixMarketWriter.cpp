#include "TestMatrixMarketWriter.h"


void TestMatrixMarketWriter::buildUp(std::string fileName)
{
    testOutputFileName = fileName;
    testMatrix = {1.1,0,0,0,0,2.2,3.3,0,0,0,0,0,4.4,0,0,5.21,0,0,0,0,0,0,0};
    objectToBeTested = new MatrixMarketWriter(fileName,&testMatrix[0],testMatrix.size());
}

void TestMatrixMarketWriter::buildUp(std::string fileName, double cutoff)
{
    testMatrix = {1.1,0,0,0,0,2.2,3.3,0,0,0,0,0,4.4,0,0,5.21,0,0,0,0,0,0,0};
    objectToBeTested = new MatrixMarketWriter(fileName,cutoff,&testMatrix[0],testMatrix.size());
}

void TestMatrixMarketWriter::writeSuccessFileNotempty()
{
    fileChecker.open(testOutputFileName);
    assert(fileChecker.peek() != std::ifstream::traits_type::eof());
    fileChecker.close();
}

void TestMatrixMarketWriter::tearDown(std::string fileName)
{
    fileChecker.close();
    remove(fileName.c_str());
    delete objectToBeTested;
}