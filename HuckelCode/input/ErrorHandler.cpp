#include "ErrorHandler.h"
#include <iostream>

void ErrorHandler::invalidNumberOfArguments()
{
    std::cerr << "Please provide two filenames: one for your molecule file and one for your configuration file."  << std::endl;
}

void ErrorHandler::notAvailableFormats(const char * inputFormat, const char * outputFormat)
{
    std::cerr << "Formats not currently available in OpenBabel: " << inputFormat << " or " << outputFormat << std::endl;
}

void ErrorHandler::openBabelInputFileNotOpen(std::string fileName)
{
    std::cerr << "Cannot open input file(" <<  fileName << ") for conversion with OpenBabel, please check if it exists in the current folder." << std::endl;
}

void ErrorHandler::openBabelOutputFileNotOpen()
{
    std::cerr << "Cannot open output file for conversion with OpenBabel, please check if you have rights to write the file." << std::endl;
}

void ErrorHandler::yahInputFileNotOpen()
{
    std::cerr << "Cannot open PDB file to be converted, please check if it exists in the current folder." << std::endl;
};

void ErrorHandler::yahOutputFileNotOpen()
{
    std::cerr << "Cannot open output file for .yah conversion, please check if you have rights to write the file." << std::endl;
}

void ErrorHandler::noAtomsFound()
{
    std::cerr << "No atoms found in the PDB file. Please check if it has correct lines starting with the ATOM or HETATM keyword." << std::endl;
}

void ErrorHandler::noSettingsInputFile(std::string filename)
{
    std::cerr << "Configuration file " << filename << " does not exist." << std::endl;
}

void ErrorHandler::settingsOutputFileNotOpen(std::string filename)
{
    std::cerr << "Cannot open " << filename << " for writing the properties, please check if it exists in the current folder." << std::endl;
}
