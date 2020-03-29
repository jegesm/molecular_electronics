cmake_minimum_required(VERSION 2.6)
project(input_pdb)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

set(SOURCE_FILES main.cpp PDBParser.cpp PDBParser.h SettingsParser.cpp SettingsParser.h OpenBabelUtils.cpp OpenBabelUtils.h ErrorHandler.h ErrorHandler.cpp)
add_executable(input_pdb ${SOURCE_FILES})

target_link_libraries(input_pdb /home/toncsi/openbabel/lib/libopenbabel.so)
