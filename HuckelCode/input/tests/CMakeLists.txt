cmake_minimum_required(VERSION 2.6)
project(input_pdb)

set(CMAKE_CXX_STANDARD 14)

include_directories(../include)

set(SOURCE_FILES main.cpp TestOpenBabelUtils.cpp TestOpenBabelUtils.h TestPdbParser.cpp TestPdbParser.h TestSettingsParser.cpp TestSettingsParser.h)
add_executable(tests ${SOURCE_FILES} ../PDBParser.cpp ../PDBParser.h ../SettingsParser.cpp ../SettingsParser.h
        ../OpenBabelUtils.cpp ../OpenBabelUtils.h ../ErrorHandler.h ../ErrorHandler.cpp)

target_link_libraries(tests /home/toncsi/openbabel/lib/libopenbabel.so)
