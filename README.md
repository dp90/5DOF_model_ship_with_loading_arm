# 5DOF_model_ship_with_loading_arm
Nonlinear 5 d.o.f. dynamic model of a ship with loading arm loaded by hydrostatic forces in Matlab.

The Powerpoint file shows a mechanical scheme of the situation studied. 
The Maple files contain derivations for equations that are the input in the Matlab functions.

The Matlab files perform the following functions:
1. Main.m solves the equations for the motions of the vessel.
2. compute_loads.m computes the loads generated by the waves acting on the vessel.
3. compute_stresses.m computes hydrostatic pressures and is used only as a verification.
4. solve_statespace_vector.m computes the state-space vector of the degrees of freedom at each time step.
5. determine_vessel_motion.m combines the forces to apply in the solve_spacestate_vector function to determine the actual motion of the vessel. 
6. Simulation.m can be run after Main.m to show a simulation of the vessel's motion.
