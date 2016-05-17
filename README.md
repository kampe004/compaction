# README #

This directory contains a dry firn densification model, written in C++, which is used to evaluate different parametrizations of snow metamorphism, compaction and fresh snow density. 

## How do I get set up? ##

First, download the code from github:

    git clone git@github.com:kampe004/compaction.git
 
Compile the firn model: 

    mkdir build && cd build
    cmake -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ -D CMAKE_CXX_FLAGS=-fPIC ..
    make -j

## Running the model ##

Prepare an INI file with settings. Find example INI under examples/

    ./dfdm.exe settings.ini

## Contribution guidelines ##

Please use <a href="http://gitimmersion.com/lab_24.html" target="_blank">git branches</a> for feature development. 

## Who do I talk to? ##

Leo van Kampenhout (L.vankampenhout@uu.nl)
