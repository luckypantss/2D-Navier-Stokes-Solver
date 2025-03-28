clear 
clc
close all
isPlot = true;      % Turn to false to disable plots outside of simulation
runOnGPU = false;    % Turn to false to run on the CPU
%% Physical properties:
% For room temperature at 20deg celsius. Data taken from "Fysika".
rho = 997;          % Density                   [kg/m^3]
mu = 1057e-6;       % Dynamic Viscocity         [Pa*s]
nu = mu/rho;        % Kinematic Viscocity       [m^2/s]
F_y = - 9.823;      % Gravitatonal force        [m/s^2]

%% Discretizinfg the domain 
% Length in x and y directrions respectevley    [m]
L_y = 1;
L_x = 1;

% Number of gridpoints to be assigned           [-]
n_y = 10;
n_x = 10;

[X, Y, dx, dy] = SpatialDisc(L_x, L_y, n_x, n_y, isPlot, runOnGPU);

%% Trying Baysein to get values for Beta and dt

%[beta, dt] = BayesianParameterSweep(X, Y, rho, nu, F_y, runOnGPU, dx, dy);

%% Initializing the 2D Navier Stokes solver using FDM:
% Setting up the time and time steps used in the solver
close all 
t       = 10000;                               % Simulation time step      [s]
%D_h     = 4*L_y*L_x/(2*(L_x+L_y));          % Hydraluic diameter        [m]
beta    = 0.015;          
%% 
dt      = 0.01;
tic
%NSsolver(X, Y, dx, dy, t, dt, rho, nu, F_y, D_h, isPlot, runOnGPU);
toc
    
%NSsolverTest(X, Y, dx, dy, t, dt, rho, nu, F_y, D_h, Beta, isPlot, runOnGPU);

[X,Y] = meshgrid(linspace(0,99,10));
%a = NSsolver(X, Y, dx, dy, t, dt, rho, nu, F_y, D_h, 1,1);

%[U_max] = NSsolverTest2(X, Y, dx, dy, t, dt, rho, nu, F_y, beta, D_h, isPlot);
 NSsolver(X,Y,1,1,2000,1,1,0.05,0,0.0125,1,1);
%%

[X,Y] = meshgrid(linspace(0,99,25));
U = NSsolverTest2(X,Y,1,1,10000,1,1,0.05,0,0.0125,1,1);

%% Fetching data for comparison 
fluidData = readtable("ghia_data.txt");
newVariableNames = {'x', 'Re100', 'Re400', 'Re1000', 'Re3200', 'Re5000', 'Re7500', 'Re10000'};

fluidTable = array2table(fluidData);
fluidTable.Properties.VariableNames = newVariableNames;
