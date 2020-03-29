#ifndef INPUT_PDB_TESTSETTINGSPARSER_H
#define INPUT_PDB_TESTSETTINGSPARSER_H

#include "../SettingsParser.h"
#include "assert.h"

class TestSettingsParser {
private:
    std::map<std::string, std::string>& testProps;
    SettingsParser* objectToBeTested;
    std::ifstream fileChecker;
public:
    TestSettingsParser(std::map<std::string, std::string>& props) : testProps(props) {}
    void buildUp(std::string fileName);
    void tearDown(std::string fileName1,std::string fileName2);

    void readSettings() {objectToBeTested->readSettings();}
    void writeSettings(std::string outputFileName) {objectToBeTested->writeSettings(outputFileName);}

    void readSettingsSuccessPropsNotEmpty();
    void writeSettingsSuccessFileNotEmpty();
};


#endif //INPUT_PDB_TESTSETTINGSPARSER_H
