%Parameters for reversal only model

clear all
close all

t0=0; % Start time of simulation
tN=100; % End time
dt=0.1; % time step in numerical simulation (lambda dt should be << 1 for numerical convergence)
Ncells=100; %Number of cells

%work with dimensional values 
Dr=0.01; %rotational diffusion measured in units of ***
vs=1;%Speed measured in units of ***

ChiR=0.5; %Chemotactic bias
lambdaRb=0.6; %Basal reversal rate
