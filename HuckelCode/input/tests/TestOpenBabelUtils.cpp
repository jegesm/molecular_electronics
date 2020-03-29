#include "TestOpenBabelUtils.h"


void TestOpenBabelUtils::buildUp(std::string filename)
{
    objectToBeTested = new OpenBabelUtils(filename);

    //Open PDB mock input file that is going to be used for the test
    std::ofstream tempOutFile(filename);
    tempOutFile.close();
}

void TestOpenBabelUtils::tearDown(std::string fileName)
{
    //Delete previously created file
    remove(fileName.c_str());
    objectToBeTested->closeInputFile();
    objectToBeTested->closeOutputFile();
    delete objectToBeTested;
}

void TestOpenBabelUtils::obtainInputFileFormatSuccessFileName()
{
    //Getting first part of InputFileName, cutting ".<format>"
    std::string unformattedFileName1 =
                    objectToBeTested->getInputFileName().substr(0, objectToBeTested->getInputFileName().length()-
                    objectToBeTested->getInputFileName().find_last_of("."));

    //Getting first part of PdbFileName, cutting ".pdb"
    std::string unformattedFileName2 =
                    objectToBeTested->getPdbFileName().substr(0, objectToBeTested->getPdbFileName().length()-
                    objectToBeTested->getPdbFileName().find_last_of("."));

    assert(unformattedFileName1==unformattedFileName2);
}

void TestOpenBabelUtils::obtainInputFileFormatSuccessFormat()
{
    assert(objectToBeTested->getPdbFileName().substr(objectToBeTested->getPdbFileName().find_last_of("."),4)==".pdb");
}

void TestOpenBabelUtils::openBabelConversionSuccessInputClosed()
{
    assert(!objectToBeTested->getIsInputOpen());
}

void TestOpenBabelUtils::openBabelConversionSuccessOutputClosed()
{
    assert(!objectToBeTested->getIsOutputOpen());
}

void TestOpenBabelUtils::initInPutFileSuccessFileOpen(std::string filename)
{
    objectToBeTested->initInputFile(filename);
    assert(objectToBeTested->getIsInputOpen());
}

void TestOpenBabelUtils::initOutPutFileSuccessFileOpen(std::string filename)
{
    objectToBeTested->initOutPutFile(filename);
    assert(objectToBeTested->getIsOutputOpen());
}