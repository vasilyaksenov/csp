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

In repository directory run:

```sh
git submodule update --init --recursive
```

To generate build files including dependencies in a new
subdirectory called 'build', run:

```sh
cmake -S. -Bbuild
```
This will take some time as cmake will download OR-Tools and all required dependencies.
Then build with:

```sh
cmake --build build --config Release
```

# Run example

To run the example navigate to the target directory:

```sh
cd build/RELEASE
```

And run:

```sh
csp.exe -f example_data/example_0.csv -w 2300
csp.exe -f example_data/example_1.csv -w 2300
csp.exe -f example_data/example_2.csv -w 3000
```

Solution example:

```sh
[406 * 3, 1623 * 1] : [ 159 ]
[1084 * 1, 1575 * 1, 333 * 1] : [ 8 ]
[1084 * 1, 1575 * 1, 333 * 1] : [ 8 ]
[2084 * 1, 824 * 1, 23 * 4] : [ 0 ]
```

Here, each line is a cutting pattern for one workpiece, the length of each is 3000.
In total, 4 blanks are used in the solution.
The first blank here is cut to 3 pieces of length 406 and 1 piece of length 1623,
the last value, separated by colon, 159 - is waste.
The second blank is cut to 1 piece of length 1084, 1 piece of length 1575, 1 piece of length 333 and 8 waste.
The third blank is cut in the same way as the second.
The last blank is cut to 1 piece of length 2084, 1 piece of length 824, 4 pieces of length 23 and have no waste.




