#include "TestOpenBabelUtils.h"
#include "TestPdbParser.h"
#include "TestSettingsParser.h"

int main(int argc, char** argv)
{
    //Values needed for testing
    std::string sampleOtherFormatterFileName = "some.sdf";
    std::string samplePDBFileName = "some.pdb";

    //Testing OpenBabelUtils class
    TestOpenBabelUtils testObabelObject;
    testObabelObject.buildUp(sampleOtherFormatterFileName);

    //Checking if the init input and init output functions open the files
    testObabelObject.initInPutFileSuccessFileOpen(sampleOtherFormatterFileName);
    testObabelObject.initOutPutFileSuccessFileOpen(sampleOtherFormatterFileName);

    //Checking if the file format is well acquired
    testObabelObject.obtainInputFileFormat();
    testObabelObject.obtainInputFileFormatSuccessFileName();
    testObabelObject.obtainInputFileFormatSuccessFormat();

    //Checking if files are closed after OpenBabel conversion
    testObabelObject.openBabelConversion();
    testObabelObject.openBabelConversionSuccessInputClosed();
    testObabelObject.openBabelConversionSuccessOutputClosed();

    //Checking if convert function uses the previously unit tested functions well
    testObabelObject.convert();

    testObabelObject.tearDown(sampleOtherFormatterFileName);

    //Checking if it is already PDB
    testObabelObject.buildUp(samplePDBFileName);
    testObabelObject.convert();
    testObabelObject.tearDown(samplePDBFileName);

    //Testing PdbParser class
    TestPdbParser testPdbParserObject;
    testPdbParserObject.buildUp(samplePDBFileName,"somemolecule");

    //After parsing the PDB file, checking if the .yah file was created and is non-empty
    testPdbParserObject.convert();
    testPdbParserObject.readSuccessFileExists();
    testPdbParserObject.readSuccessFileNotempty();
    testPdbParserObject.tearDown(samplePDBFileName);

    //Testing the SettingsParser class
    std::string sampleConfigFileName = "someconfig.txt";
    std::string sampleOutputFileName = "someoutputfile.yah";

    std::map<std::string, std::string> testProps;
    TestSettingsParser testSettingsParserObject(testProps);
    testSettingsParserObject.buildUp(sampleConfigFileName);

    //Testing if properties were parsed from configuration file
    testSettingsParserObject.readSettings();
    testSettingsParserObject.readSettingsSuccessPropsNotEmpty();

    //Testing if file was created in the process of appending
    testSettingsParserObject.writeSettings(sampleOutputFileName);
    testSettingsParserObject.writeSettingsSuccessFileNotEmpty();
    testSettingsParserObject.tearDown(sampleConfigFileName,sampleOutputFileName);

    std::cout << "Tests were successful" << std::endl;
    return EXIT_SUCCESS;
}
