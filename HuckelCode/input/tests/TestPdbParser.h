#ifndef INPUT_PDB_TESTPDBPARSER_H
#define INPUT_PDB_TESTPDBPARSER_H

#include "../PDBParser.h"

class TestPdbParser {
private:
    PDBParser* objectToBeTested;
    std::ifstream fileChecker;
public:

    //Creating test PDB file
    void buildUp(std::string filename, std::string molecule);

    //Destroying test PDB file
    void tearDown(std::string filename);
    void convert() {objectToBeTested->convert();}

    //Check if the output file was written after read() was invoked
    void readSuccessFileExists();
    void readSuccessFileNotempty();
};


#endif //INPUT_PDB_TESTPDBPARSER_H
