# Bachelor_source_code
This is the source code of a reaction-diffusion differential equation, coded in c++ using amandus  ( https://bitbucket.org/guidokanschat/amandus/src/master/ ) and deal.II ( https://www.dealii.org/ )

The important code is hidden in the *.cc files, the physics of the system is defined in the explicit.h, implicit.h, matrix.h.
The outpuut of the compiled program depends dynamically on the parameter files and outputs *.vtk files for a 3D Visualization.

Here are some examples of one reactant starting randomly distributed in some area resulting in a steady-state pattern. The visualitation is done with the use of VisIT.
Initial state:
![First]("theta0_5gamma1000zeit0balken.png")
After 30 timesteps:
![First]("theta0_5gamma1000zeit30balken.png")
After 300 timesteps:
![First]("theta0_5gamma1000zeit300balken.png")
