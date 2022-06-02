#!/bin/bash

DLD=wget

if [ ! -e bin ]; then
  mkdir bin include lib64 share
  ln -s lib64 lib
fi

# fmt https://fmt.dev/latest/index.html

if [ ! -f fmt-8.1.1.zip ]; then
  echo "Download fmt ... "
  ${DLD} https://github.com/fmtlib/fmt/releases/download/8.1.1/fmt-8.1.1.zip
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
  ${DLD} https://github.com/nlohmann/json/releases/download/v3.7.3/json.hpp
  mkdir -p include
  cp json.hpp include
  echo "done."
fi

# cppitertools https://twitter.com/cppitertools

if [ ! -f cppitertools-2.1.zip ]; then
  echo "Download cppitertools ... "
  ${DLD} -O cppitertools-2.1.zip https://github.com/ryanhaining/cppitertools/archive/refs/tags/v2.1.zip 
  echo "done."
fi

if [ ! -f ITERTOOLS_FLAG ]; then
  echo "Build cppitertools ... "
  unzip -qo cppitertools-2.1.zip
  cd ./cppitertools-2.1
  mkdir -p build
  cd build
  pwd
  cmake -DCMAKE_INSTALL_PREFIX=../.. ..
  cmake --build .
  cmake --install .
  cd ../..
  touch ITERTOOLS_FLAG
  echo "done."
fi

# NLOpt https://nlopt.readthedocs.io/en/latest/
if [ ! -f nlopt-2.7.1.tar.gz ]; then
  echo "Download NLOpt ... "
  ${DLD} -O nlopt-2.7.1.tar.gz https://github.com/stevengj/nlopt/archive/v2.7.1.tar.gz
  echo "done."
fi

if [ ! -f NLOPT_FLAG ]; then
  echo "Build NLOpt ... "
  tar xzf nlopt-2.7.1.tar.gz
  cd ./nlopt-2.7.1
  mkdir -p build
  cd build
  pwd
  cmake -DCMAKE_INSTALL_PREFIX=../.. ..
  cmake --build .
  cmake --install .
  cd ../..
  touch NLOPT_FLAG
  echo "done."
fi

# DataFrame https://github.com/hosseinmoein/DataFrame
if [ ! -f DataFrame-1.20.0.tar.gz ]; then
  echo "Download DataFrame ... "
  ${DLD} -O DataFrame-1.20.0.tar.gz https://github.com/hosseinmoein/DataFrame/archive/refs/tags/1.20.0.tar.gz
  echo "done."
fi

if [ ! -f DATAFRAME_FLAG ]; then
  echo "Build DataFrame ... "
  tar xzf DataFrame-1.20.0.tar.gz
  cd ./DataFrame-1.20.0
  mkdir -p build
  cd build
  pwd
  cmake -DCMAKE_INSTALL_PREFIX=../.. -DCMAKE_BUILD_TYPE=Release ..
  cmake --build .
  cmake --install .
  cd ../..
  touch DATAFRAME_FLAG
  echo "done."
fi

