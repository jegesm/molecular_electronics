#include <unordered_map>
#include <cstdio>
#include <omp.h>

#include <iostream>
#include <assert.h>
#include "TestMatrixMarketWriter.h"

int main()
{
    std::string sampleOutputFileName = "someoutput.bin";

    //Testing MatrixMarketWriter class without using cutoff value
    TestMatrixMarketWriter testObject;
    testObject.buildUp(sampleOutputFileName);

    //Creating output binary file and checking if the file was not left empty
    testObject.createBinaryFile();
    testObject.writeSuccessFileNotempty();

    testObject.tearDown(sampleOutputFileName);

    //New test for dump function
    testObject.buildUp(sampleOutputFileName);

    //Creating dump file and checking if the file was not left empty
    testObject.createDumpFile();
    testObject.writeSuccessFileNotempty();

    testObject.tearDown(sampleOutputFileName);

    //Repeating previous tests for MatrixMarketWriter class with cutoff value
    testObject.buildUp(sampleOutputFileName,0.5);

    testObject.createBinaryFile();
    testObject.writeSuccessFileNotempty();
    testObject.tearDown(sampleOutputFileName);

    testObject.buildUp(sampleOutputFileName,0.5);

    testObject.createDumpFile();
    testObject.writeSuccessFileNotempty();

    testObject.tearDown(sampleOutputFileName);
    std::cout << "Tests were successful!" << std::endl;
    return 0;
}

