#include "OpenBabelUtils.h"
#include "SettingsParser.h"
#include "PDBParser.h"

int main(int argc, char** argv)
{
        //Check number of arguments given by the user
        if(argc < 3)
        {
            ErrorHandler::invalidNumberOfArguments();
            return EXIT_FAILURE;
        }

        //Convert file to PDB format if not given so
        OpenBabelUtils toPdbConversion(argv[1]);
        if(toPdbConversion.convert())
        {
            return EXIT_FAILURE;
        }

        //Getting information from properties file
        std::map<std::string, std::string> props;
        SettingsParser settingsParser{argv[2], props};
        if(settingsParser.readSettings())
        {
            return EXIT_FAILURE;
        }

        //Parsing molecule from PDB file
        PDBParser moleculeParser(toPdbConversion.getPdbFileName(),settingsParser.getInputMolecule());
        if(moleculeParser.convert())
        {
            return EXIT_FAILURE;
        }

        //Creating outputfile with the parsed molecule and the given properties
        if(settingsParser.writeSettings(moleculeParser.getOutputFileName()))
        {
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
}