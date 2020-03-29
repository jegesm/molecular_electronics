#ifndef INPUT_PDB_ERRORHANDLER_H
#define INPUT_PDB_ERRORHANDLER_H

#include <iostream>

namespace ErrorHandler {
    //Handling error of unknown filename, too few arguments were given
    extern void invalidNumberOfArguments();

    //Handling error of unknown formats given for OpenBabel
    extern void notAvailableFormats(const char * inputFormat, const char * outputFormat);

    //Informing user about inability to open input file for OpenBabel conversion

    extern void openBabelInputFileNotOpen(std::string fileName);

    //Informing user about inability to open output file for OpenBabel conversion
    extern void openBabelOutputFileNotOpen();

    //Informing user about inability to open input file for .yah conversion
    extern void yahInputFileNotOpen();

    //Informing user about inability to open output file for .yah conversion
    extern void yahOutputFileNotOpen();

    //No atoms found in PDB file
    extern void noAtomsFound();

    //Informing user about non existent settings file
    extern void noSettingsInputFile(std::string filename);

    //Informing user about inability to open settings output file
    extern void settingsOutputFileNotOpen(std::string);

};


#endif //INPUT_PDB_ERRORHANDLER_H
