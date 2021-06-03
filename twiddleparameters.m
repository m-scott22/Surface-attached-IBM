%Parameters for twiddle only model

clear all
close all

t0=0; % Start time of simulation
tN=100; % End time
dt=0.1; % time step in numerical simulation (lambda dt should be << 1 for numerical convergence)
Ncells=100; %Number of cells

%work with dimensional values 
Dr=0; %rotational diffusion measured in units of ***
vs=1;%Speed measured in units of ***

ChiT=0.5; %chemotactic bias
kappa=0.5; %twiddle exit bias
lambdaTb=0.6; %basal twiddle rate


