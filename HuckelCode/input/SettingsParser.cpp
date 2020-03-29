#include <vector>
#include "SettingsParser.h"

int SettingsParser::readSettings()
{
    std::ifstream input_file(configFileName.c_str());
    std::string delimiter = "=";
    std::string line;
    if (input_file.good())
    {
        while (std::getline(input_file, line))
        {
            int pos = line.find(delimiter);
            std::string propKey = line.substr(0, pos);
            std::string propValue = line.substr(pos + 1, line.length());
            if(propKey=="molecule name")
                inputMolecule=propValue;
            properties[propKey] = propValue;
        }
        input_file.close();
        return EXIT_SUCCESS;
    }
    else {
        ErrorHandler::noSettingsInputFile(configFileName);
        return EXIT_FAILURE;
    }
}

int SettingsParser::writeSettings(std::string outputFileName)
{
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::fstream::app);
    outputFile << std::endl;
    outputFile << ";Printing options" << std::endl;
    outputFile << "Print" << std::endl;
    std::vector<std::string> dumps;
    if (outputFile.is_open())
    {
        for (auto const &item : properties) {
            if (item.second == "true" && item.first.find("dump") == std::string::npos)
            {
                outputFile << item.first << std::endl;
            }
            else if(item.second == "true")
            {
                dumps.push_back(item.first);
            }
        }
        outputFile << "End_Print \n" << std::endl;
        for(auto const &item : dumps)
        {
            outputFile << item << std::endl;
        }
        outputFile.close();
        std::cout << "Successfully created the yah file." << std::endl;
        return EXIT_SUCCESS;
    }
    else
    {
        ErrorHandler::settingsOutputFileNotOpen(outputFileName);
        return EXIT_FAILURE;
    }
}

const std::string &SettingsParser::getInputMolecule() const {
    return inputMolecule;
}

const std::string &SettingsParser::getConfigFileName() const {
    return configFileName;
}

