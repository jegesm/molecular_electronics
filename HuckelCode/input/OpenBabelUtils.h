#ifndef INPUT_PDB_OPENBABELUTILS_H
#define INPUT_PDB_OPENBABELUTILS_H

#include <iostream>
#include <fstream>
#include <openbabel/obconversion.h>
#include <string>
#include <locale>


class OpenBabelUtils {
private:
    std::string inputFileName;
private:
    std::string outputFileName;

    std::string inputFileFormat;
    std::string pdbFileName;

    std::ifstream inputFile;
    std::ofstream outputFile;
public:
    OpenBabelUtils(std::string fileName) : inputFileName(fileName) {}
    int convert();
    int OpenBabelConversion();
    int initInputFile(std::string filename);
    int initOutPutFile(std::string filename);
    void obtainInputFileFormat();

    void closeInputFile() {inputFile.close();}
    void closeOutputFile() {outputFile.close();}

    bool getIsInputOpen() {return inputFile.is_open();}
    bool getIsOutputOpen() {return inputFile.is_open();}

    const std::string &getPdbFileName() const;
    const std::string &getInputFileName() const;
    const std::string &getOutputFileName() const;
};


#endif //INPUT_PDB_OPENBABELUTILS_H
