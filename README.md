# Finite Element Solver of a Reaction-Diffusion differential equation
This is the source code of a reaction-diffusion differential equation, coded in c++ using amandus  ( https://bitbucket.org/guidokanschat/amandus/src/master/ ) and deal.II ( https://www.dealii.org/ ). It uses the discontinuous Galerkin Finite Element Method together with the Rothe-Method and Crank-Nicolson. Furthermore the equations are solved with the use of a solver using the GMRES-Method.

For the actual equation, specific coefficients and further questions please do not hesitate to contact me.

The important code is hidden in the *.cc files, the physics of the system is defined in the explicit.h, implicit.h, matrix.h.
The outpuut of the compiled program depends dynamically on the parameter files and outputs *.vtk files for a 3D Visualization.

Here are some examples of one reactant starting randomly distributed in some area resulting in a steady-state pattern. The visualitation is done with the use of VisIT.
Initial state:
![Initial State](theta0_5gamma1000zeit0balken.png)
After 30 timesteps:
![30 timesteps](theta0_5gamma1000zeit30balken.png)
After 300 timesteps:
![300 timesteps](theta0_5gamma1000zeit300balken.png)
