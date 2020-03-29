#include "PDBParser.h"

template<ESBTL::Reading_mode mode,class mandatoryFields, class Line_selector,class Builder,class Occupancy_handler>
bool PDBParser::readPdb(const std::string& filename,Line_selector& sel,Builder& builder,const Occupancy_handler& occupancy,char altloc){
    ESBTL::Line_reader<ESBTL::PDB::Line_format<mandatoryFields>,Line_selector,Builder> lineReader(sel,builder);
    return lineReader.template read<mode>(filename,occupancy,altloc);
}

template<class mandatoryFields, class Line_selector,class Builder,class Occupancy_handler>
bool PDBParser::readPdb(const std::string& filename,Line_selector& sel,Builder& builder,const Occupancy_handler& occupancy,char altloc){
    return readPdb<ESBTL::ASCII, mandatoryFields>(filename,sel,builder,occupancy,altloc);
}

bool PDBParser::readPdb()
{
    return readPdb<ESBTL::ASCII, CustomMandatoryFields>(inputFileName, line_selector, builder, Accept_none_occupancy_policy());
}

bool PDBParser::convert()
{
    outputFile.open(outputFileName);
    if (outputFile.is_open())
    {
        if (readPdb())
        {
            if (systems.empty() || systems[0].has_no_model())
            {
                ErrorHandler::noAtomsFound();
                remove(outputFileName.c_str());
                return EXIT_FAILURE;
            }
            outputFile << moleculeName << std::endl;
            outputFile << "Molecular" << std::endl;
            outputFile << "Geometry" << std::endl;
            outputFile << line_selector.getNblines() << std::endl;

            //Consider only the first model of the first system
            const ESBTL::Default_system::Model &model = *systems[0].models_begin();
            for (ESBTL::Default_system::Model::Atoms_const_iterator it_atm = model.atoms_begin();
                 it_atm != model.atoms_end(); ++it_atm) {
                outputFile << it_atm->atom_name().at(0) << " ";
                outputFile << it_atm->x() << " " << it_atm->y() << " " << it_atm->z() << std::endl;
            }
            outputFile << "Charge" << std::endl;
            outputFile << "0" << std::endl;
            outputFile.close();
        }
        else
        {
            //Error messaging made through ESBTL
            return EXIT_FAILURE;
        }
    }
    else
    {
        ErrorHandler::yahOutputFileNotOpen();
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
