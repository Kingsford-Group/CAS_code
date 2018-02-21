#!/bin/bash

# Get SDSL
if [ ! -d "sdsl-lite" ]; then
  git clone https://github.com/simongog/sdsl-lite
  cd sdsl-lite
  ./install.sh .
cd ..
fi
# Get seqan
if [ ! -d "seqan" ]; then
  git clone https://github.com/seqan/seqan
  cd seqan
  git checkout develop
cd ..
fi
# Compile release folder
if [ ! -d "release" ]; then
  mkdir release
  cd release
  cmake .. \
    -DCMAKE_PREFIX_PATH="seqan/util/cmake" \
    -DSEQAN_INCLUDE_PATH="seqan/include"
  make
  cd ..
fi
# Compile debug folder
if [ ! -d "debug" ]; then
  mkdir debug
  cd debug
  cmake .. \
    -DCMAKE_PREFIX_PATH="seqan/util/cmake" \
    -DSEQAN_INCLUDE_PATH="seqan/include" \
    -DCMAKE_BUILD_TYPE=Debug
  make
fi
# Print welcome message
echo -e "\e[1;31mAll compiled and ready for \e[1;96mMingfu \e[1;31mto test!!!!!\e[0m"

