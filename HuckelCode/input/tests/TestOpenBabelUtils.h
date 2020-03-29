#ifndef INPUT_PDB_TESTOPENBABELUTILS_H
#define INPUT_PDB_TESTOPENBABELUTILS_H

#include "../OpenBabelUtils.h"

class TestOpenBabelUtils
{
private:
    OpenBabelUtils* objectToBeTested;
public:

    void buildUp(std::string filename);
    void tearDown(std::string fileName);

    //Redefining private member's functions due to inaccessibility
    void obtainInputFileFormat()  { objectToBeTested->obtainInputFileFormat(); }
    void openBabelConversion()  { objectToBeTested->OpenBabelConversion(); }

    //Check if convert function failed in execution
    void convert() { assert(!objectToBeTested->convert()); }

    //Redefining I/O functions for accessibility and testing purposes
    void initInPutFileSuccessFileOpen(std::string filename);
    void initOutPutFileSuccessFileOpen(std::string filename);

    //Testing OpenBabelUtils obtainInputFileFormat function
    void obtainInputFileFormatSuccessFileName();
    void obtainInputFileFormatSuccessFormat();

    //Testing OpenBabelUtils openBabelConversion function (most of it consists of OpenBabel API calls)
    void openBabelConversionSuccessInputClosed();
    void openBabelConversionSuccessOutputClosed();
};


#endif //INPUT_PDB_TESTOPENBABELUTILS_H
