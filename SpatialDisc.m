function [X,Y,dx,dy] = SpatialDisc(L_x,L_y,n_x,n_y, isPlot,runOnGPU)
% SpatialDisc - Discretizes the spatial domain into grid points.
%
% Inputs:
%   L_x: Length of the domain in the x-direction.
%   L_y: Length of the domain in the y-direction.
%   n_x: Number of desired grid points in the x-direction.
%   n_y: Number of desired grid points in the y-direction.
%   isPlot: Boolean indicating whether to plot the grid points.
%   runOnGPU: Boolean indicating whether to use GPU computation.
%
% Outputs:
%   X: Grid points along the x-axis.
%   Y: Grid points along the y-axis.
%   dx: Spatial step size in the x-direction.
%   dy: Spatial step size in the y-direction.
%
% Usage:
%   [X, Y, dx, dy] = SpatialDisc(L_x, L_y, n_x, n_y, isPlot, runOnGPU);
%
% Description:
%   This function discretizes the spatial domain into grid points with
%   specified lengths and number of grid points along each axis. It returns
%   the grid points along with the spatial step sizes in both directions.
%   Additionally, it can plot the grid points if specified, and it supports
%   computation on either CPU or GPU.
%
% Author: Elias Larsson
% Date: 

    % Defining the domain of x and y
    x       = linspace(0, L_x, n_x);
    y       = linspace(0, L_y, n_y);
    % Spacing between the gridpoints
    dx      = L_x / (n_x-1);
    dy      = L_y / (n_y-1);
    
    % Seetting the domain into gridpoints
    [X_temp, Y_temp] = meshgrid(x,y);

    % Checking if the script is going to run on GPU or CPU
    if runOnGPU == true
        X = gpuArray(X_temp);
        Y = gpuArray(Y_temp);
    else
        X = X_temp;
        Y = Y_temp;
    end
    
    % Plotting the gridpoints
    if isPlot == true 
        figure
        scatter(X(:), Y(:), 'filled')
        xlabel('X')
        ylabel('Y')
        title(['Discretized points of x and y with ', num2str(length(X) * length(Y)), ' grid points']);
        grid minor
    end
end