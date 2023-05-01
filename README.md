# LeafSim

Elementary leaf simulation built by Finley Webb in C/C++ using GLFW

Based on FRAPc by Francois Nedelec

## How to Build

Make sure the terminal is in the leafsim directory and cmake/make is installed
Then run

```
mkdir build
cd build
cmake ..
make
```

To run, ensure you are in the /build directory

```
./leafsim
```

To run a Genetic Algorithm parameter search, copy the executable to the /GA
directory

```
cd build
cp leafsim ../GA
cd GA
./evolve.py
``

## Authors and acknowledgment

FJN, 13.11.2021
FW, 01.05.2023

## License

This is Open Source
