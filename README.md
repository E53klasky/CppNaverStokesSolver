# CppNaverStokesSolver

I will be taking a finite difference python code and converting it to cpp

Here is the code that I will be converting: https://github.com/Ceyron/machine-learning-and-simulation/blob/main/english/simulation_scripts/lid_driven_cavity_python_simple.py

How to build

First install protobuf for MGARD
Then install MGARD
Then install ADIOS2

How to run
the exe is in the same dir as it is built in so run it there as ./main

to change the parms you must change them inside the code

the it output the bp dir in the same place as the exe

then there is a python file which reads it in and does a simple plot of it

Here is the difference in time btw the debug mode and non debug mode on my computer

DEBUG mode:
Elapsed time: 1.21396 seconds

Realse mode:
Elapsed time: 0.13192 seconds
