#include "TestSettingsParser.h"

void TestSettingsParser::buildUp(std::string fileName)
{
    objectToBeTested = new SettingsParser(fileName, testProps);
    std::ofstream tempOutFile(fileName);

    tempOutFile << "something1=true " << std::endl;
    tempOutFile << "something2=true " << std::endl;
    tempOutFile << "something3=false " << std::endl;
    tempOutFile << "something4=true " << std::endl;
    tempOutFile << "something5=false " << std::endl;

    tempOutFile.close();
}

void TestSettingsParser::readSettingsSuccessPropsNotEmpty()
{
    assert(!objectToBeTested->getIsPropsEmpty());
}

void TestSettingsParser::writeSettingsSuccessFileNotEmpty()
{
    fileChecker.open(objectToBeTested->getConfigFileName());
    assert(fileChecker.peek() != std::ifstream::traits_type::eof());
    fileChecker.close();
}

void TestSettingsParser::tearDown(std::string fileName1,std::string fileName2)
{
    remove(fileName1.c_str());
    remove(fileName2.c_str());
    delete objectToBeTested;
}

