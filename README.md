# Introduction

This is yet another example of Cutting Stock Problem (CSP) solving with C++ using Google OR-Tools.
The code for this project was taken from repository https://github.com/emadehsan/csp and
"Practical Python AI Projects: Mathematical Models of Optimization Problems with Google OR-Tools"
book written by Serge Kruk. 
Original Python code:
https://github.com/sgkruk/Apress-AI/blob/master/cutting_stock.py
https://github.com/emadehsan/csp

# Project Build

This build support only x64 Windows machines.

You'll need:

* CMake >= 3.8
* A C++17 compiler

To generate build files including dependencies in a new
subdirectory called 'build', run:

```sh
cmake -S. -Bbuild
```
Be patient, this will take some time as cmake will download OR-Tools and all required dependencies.

Then build with:

```sh
cmake --build build --config Release
```