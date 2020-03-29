#ifndef INPUT_PDB_PDB_PARSER_H
#define INPUT_PDB_PDB_PARSER_H

#include <iostream>
#include <ESBTL/default.h>
#include <ESBTL/weighted_atom_iterator.h>
#include "ErrorHandler.h"

class PDBParser {
private:
    struct CustomMandatoryFields;
    typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<CustomMandatoryFields> > Accept_none_occupancy_policy;

    //Members used for file manipulation and providing information for the output file
    std::string inputFileName;
    std::string moleculeName;
    std::string outputFileName;

    //Output file containing the .yah file format
    std::ofstream outputFile;

    ESBTL::PDB_line_selector line_selector;
    std::vector<ESBTL::Default_system> systems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system> builder;

    template<ESBTL::Reading_mode mode,class mandatoryFields, class Line_selector,class Builder,class Occupancy_handler>
    bool readPdb(const std::string& filename,Line_selector& sel,Builder& builder,const Occupancy_handler& occupancy,char altloc=' ');

    template<class mandatoryFields, class Line_selector,class Builder,class Occupancy_handler>
    bool readPdb(const std::string& filename,Line_selector& sel,Builder& builder,const Occupancy_handler& occupancy,char altloc=' ');

    //In this struct you can specify what attributes your PDB has
    struct CustomMandatoryFields
    {
        static const bool record_name=true;
        static const bool atom_serial_number=true;
        static const bool atom_name=true;
        static const bool alternate_location=false;
        static const bool residue_name=false;
        static const bool chain_identifier=false;
        static const bool residue_sequence_number=false;
        static const bool insertion_code=false;
        static const bool x=true;
        static const bool y=true;
        static const bool z=true;
        static const bool occupancy=false;
        static const bool temperature_factor=true;
        static const bool element=false;
        static const bool charge_str=false;
        static const bool model_number=true;
    };
public:
    PDBParser(std::string fileName, std::string moleculeName) :
            line_selector(), systems(), builder(systems,line_selector.max_nb_systems()), inputFileName(fileName), moleculeName(moleculeName)
    {
        auto const pos = inputFileName.find_last_of('.');
        outputFileName = inputFileName.substr(0,pos);
        outputFileName.append(".yah");
    }

    ~PDBParser() {}

    //Reading PDB file
    bool readPdb();

    //Convert PDB file to .yah file
    bool convert();

    const std::string &getOutputFileName() const { return outputFileName; }
};



#endif //INPUT_PDB_PDB_PARSER_H
