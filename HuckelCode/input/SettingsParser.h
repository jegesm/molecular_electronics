#ifndef INPUT_PROPERTIES_SETTINGSPARSER_H
#define INPUT_PROPERTIES_SETTINGSPARSER_H

#include <iostream>
#include <map>
#include <fstream>
#include "ErrorHandler.h"

class SettingsParser {
    std::map<std::string, std::string>& properties;
    std::string configFileName;
private:
    std::string inputMolecule;
public:
    SettingsParser(std::string configFileName, std::map<std::string, std::string>& props) :
            configFileName(configFileName), properties(props) {}
    ~SettingsParser() {}

    //Gets the settings specified in the configuration file
    int readSettings();

    //Writes the settings based on the data previously obtained from the configuration file
    int writeSettings(std::string outputFileName);

    const std::string &getInputMolecule() const;
    bool getIsPropsEmpty() {return properties.empty();}
    const std::string &getConfigFileName() const;
};


#endif //INPUT_PROPERTIES_SETTINGSPARSER_H
