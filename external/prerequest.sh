#!/bin/bash

# fmt https://fmt.dev/latest/index.html


if [ ! -f fmt-8.1.1.zip ]; then
  echo "Download fmt ... "
  curl -O https://github.com/fmtlib/fmt/releases/download/8.1.1/fmt-8.1.1.zip
  echo "done."
fi

if [ ! -f FMT_FLAG ]; then
  echo "Build fmt ... "
  unzip -qo fmt-8.1.1.zip
  cd ./fmt-8.1.1
  mkdir -p build
  cd build
  pwd
  cmake -DCMAKE_INSTALL_PREFIX=../.. ..
  cmake --build .
  cmake --install .
  cd ../..
  touch FMT_FLAG
  echo "done."
fi


# JSON for modern C++ https://json.nlohmann.me/
if [ ! -f json.hpp ]; then
  echo "Download json ... "
  curl -O https://github.com/nlohmann/json/releases/download/v3.7.3/json.hpp
  mkdir -p include
  cp json.hpp include
  echo "done."
fi


