#include "OpenBabelUtils.h"
#include "ErrorHandler.h"

int OpenBabelUtils::convert()
{
    //Getting inputFileFormat from input file name
    obtainInputFileFormat();

    //if already PDB, just leave it
    if (inputFileFormat != "PDB")
    {
        //Opening input file
        if(initInputFile(inputFileName))
        {
            return EXIT_FAILURE;
        }

        //Opening output file
        if(initOutPutFile(pdbFileName))
        {
            return EXIT_FAILURE;
        }
        return OpenBabelConversion();
    }
    return EXIT_SUCCESS;
}

int OpenBabelUtils::initInputFile(std::string filename)
{
    inputFile.open(filename);
    if (!inputFile.is_open()) {
        ErrorHandler::openBabelInputFileNotOpen(filename);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int OpenBabelUtils::initOutPutFile(std::string filename)
{
    outputFile.open(filename);
    if (!outputFile.is_open()) {
        ErrorHandler::openBabelOutputFileNotOpen();
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int OpenBabelUtils::OpenBabelConversion()
{
    //Creating object responsible for conversion
    OpenBabel::OBConversion conv(&inputFile, &outputFile);
    const char *iFileFormat = inputFileFormat.c_str();
    const char *oFileFormat = "PDB";

    //Setting conversion formats, checking if they are available in OpenBabel
    if (!conv.SetInAndOutFormats(iFileFormat, oFileFormat))
    {
        ErrorHandler::notAvailableFormats(iFileFormat,oFileFormat);
        return EXIT_FAILURE;
    }
    else
    {
        int n = conv.Convert();
        std::cout << n << " molecules converted" << std::endl;
        inputFile.close();
        outputFile.close();
        return EXIT_SUCCESS;
    }
};

void OpenBabelUtils::obtainInputFileFormat()
{
    //Obtaining file format from inputfile name
    auto const pos = inputFileName.find_last_of('.');
    inputFileFormat = inputFileName.substr(pos + 1);

    //Creating output filename
    pdbFileName = inputFileName.substr(0, pos);
    pdbFileName.append(".pdb");
    //Converting format to uppercase so that it conforms to OpenBabel syntax
    std::locale loc;
    for (std::string::size_type i=0; i<inputFileFormat.length(); ++i)
        inputFileFormat.at(i) = std::toupper(inputFileFormat.at(i),loc);
}

const std::string &OpenBabelUtils::getPdbFileName() const {
    return pdbFileName;
}

const std::string &OpenBabelUtils::getInputFileName() const {
    return inputFileName;
}

const std::string &OpenBabelUtils::getOutputFileName() const {
    return outputFileName;
}
