#ifndef OUTPUT_TESTMARKETMATRIXWRITER_H
#define OUTPUT_TESTMARKETMATRIXWRITER_H

#include <vector>
#include "../output/MatrixMarketWriter.h"
#include "assert.h"

class TestMatrixMarketWriter {
private:
    MatrixMarketWriter* objectToBeTested;
    std::vector<double> testMatrix;
    std::string testOutputFileName;
    std::ifstream fileChecker;
public:
    void buildUp(std::string fileName);
    void buildUp(std::string fileName, double cutoff);

    void tearDown(std::string fileName);
    void writeSuccessFileNotempty();

    void createBinaryFile() {objectToBeTested->createBinaryFile();}
    void createDumpFile() {objectToBeTested->createDumpFile();}
};


#endif //OUTPUT_TESTMARKETMATRIXWRITER_H
