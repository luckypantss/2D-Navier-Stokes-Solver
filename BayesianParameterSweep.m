function [beta, dt] = BayesianParameterSweep(X, Y, rho, nu, F_y,runOnGPU, dx, dy)
    rng(1); % Set the seed for reproducibility

    % Define the parameter ranges
    dtRange = [0.001, 0.01];
    betaRange = [0.1, 2];
    t_max = 100;
    U_size = [t_max, size(X, 1), size(X, 2), 2]; % Assuming 2D grid for U_temp
    u_diff = zeros(2, 1);
    u_conv = zeros(2, 1);
   
    % Initialize simulation parameters
    if runOnGPU
        U_temp  = gpuArray(zeros(U_size));
        U       = gpuArray(zeros(U_size));
        P_temp  = gpuArray(zeros([t_max, size(X, 1), size(X, 2)]));
        P       = gpuArray(zeros([t_max, size(X, 1), size(X, 2)]));
    else
        U_temp  = zeros(U_size);
        U       = zeros(U_size);
        P_temp  = zeros(t_max, size(X, 1), size(X, 2));
        P       = zeros(t_max, size(X, 1), size(X, 2));
    end

    % Bayesian Optimization setup
    results = []; 

    % Main simulation loop
    for t_n = 1: t_max
        % Extract the current best parameters
        if ~isempty(results)
            bestParams  = bestPoint(results);
            dt          = bestParams.dt;
            beta        = bestParams.Beta;
        else
            % Initial parameters
            dt      = 0.005;
            beta    = 1;
        end

        % Update the simulation using the current parameters
        for i = 2: size(X, 1) - 1
            for j = 2: size(Y, 2) - 1
                for A = 1: 2
                % Convective and diffusive terms of U
                    u_conv(A, 1) = -0.5 * ( ...
                        U_temp(t_n, i, j, 1) * (U_temp(t_n, i+1, j, A) - U_temp(t_n, i-1, j, A)) / dx + ...
                        U_temp(t_n, i, j, 2) * (U_temp(t_n, i, j+1, A) - U_temp(t_n, i, j-1, A)) / dy ...
                    );
                    u_diff(A, 1) = nu * ( ...
                        (U_temp(t_n, i+1, j, A) - 2 * U_temp(t_n, i, j, A) + U_temp(t_n, i-1, j, A)) / dx^2 + ...
                        (U_temp(t_n, i, j+1, A) - 2 * U_temp(t_n, i, j, A) + U_temp(t_n, i, j-1, A)) / dy^2 ...
                    );
                end

                % Pressure force and momentum terms of P 
                p_force = 0.5 * beta * ( ...
                    (U_temp(t_n, i+1, j, 1) - U_temp(t_n, i-1, j, 1)) / dx + ...
                    (U_temp(t_n, i, j+1, 2) - U_temp(t_n, i, j-1, 2)) / dy ...
                );
                p_moment = -0.5 * ( ...
                    U_temp(t_n, i, j, 1) * (P(t_n, i+1, j) - P(t_n, i-1, j)) / dx + ...
                    U_temp(t_n, i, j, 2) * (P(t_n, i, j+1) - P(t_n, i, j-1)) / dy ...
                );

                % Updateing pressure field
                P_temp(t_n + 1, i, j) = P_temp(t_n, i, j) + p_force + p_moment;

                % Update velocity components
                U_temp(t_n + 1, i, j, 1)    = U_temp(t_n, i, j, 1) + dt * ( ...
                    u_conv(1, 1) + u_diff(1, 1) - ...
                    1/(2*rho*dx) *(P(t_n, i+1, j) - P(t_n, i-1, j)) + ...
                    1 / rho * F_y ...
                    );
                
                U_temp(t_n + 1, i, j, 2)    = U_temp(t_n, i, j, 2) + dt * ( ...
                    u_conv(2, 1) + u_diff(2, 1) - ...
                    1/(2*rho*dy) *(P(t_n, i, j+1) - P(t_n, i, j-1)) + ...
                    1 / rho * F_y ...
                    );
            end
        end
        % Applying the boundary conditions
        U_temp      = applyNoSlipBoundary(U_temp, t_n);
        P_temp      = applyNeumannBoundary(P_temp, t_n);

        % Updating the matrices from the 
        U   = U_temp;
        P   = P_temp;


        % Calculate error metric for current timestep
        error = calculate_error(U, P);

        % Update Bayesian Optimization model
        if mod(t_n, 10) == 0 % Run optimization every 10 steps
            results = optimize_params(dtRange, betaRange, error);
        end
    end

    % Return the best Beta and dt values
    if ~isempty(results)
        bestParams = bestPoint(results);
        beta = bestParams.Beta;
        dt = bestParams.dt;
    else
        beta = beta;
        dt = dt;
    end
end

function error = calculate_error(U_temp, P)
    % Calculate an error metric based on U_temp and P
    % Implement error calculation
    error = norm(U_temp(:)) + norm(P(:)); % Example metric, adjust as needed
end

function results = optimize_params(dtRange, betaRange, currentError)
    % Bayesian Optimization setup
    optimVars = [
        optimizableVariable('dt', dtRange),
        optimizableVariable('Beta', betaRange)
    ];

    % Objective function for Bayesian Optimization
    objFunc = @(params) simulation_error(params.dt, params.Beta, currentError);

    % Run Bayesian Optimization
    results = bayesopt(objFunc, optimVars, 'Verbose', 0, 'UseParallel', false);
end

function error = simulation_error(dt, Beta, currentError)
    % Placeholder function to simulate error calculation for Bayesian Optimization
    % In real implementation, this would rerun the simulation with new parameters
    % For simplicity, let's assume the error improves over time in a mock manner
    error = currentError * (0.95 + 0.1 * rand); % Simulated improvement
end

% Apply no-slip boundary conditions
function [U_bound] = applyNoSlipBoundary(U_field, t_n)
    U_field(t_n + 1, 1, :, :)           = 0; % Top boundary
    U_field(t_n + 1, end, :, :)         = 0; % Bottom boundary
    U_field(t_n + 1, :, 1, :)           = 0; % Left boundary
    U_field(t_n + 1, :, end, :)         = 0; % Right boundary
    U_bound = U_field;
end

% Apply Neumann boundary conditions (zero gradient)
function [P_bound] = applyNeumannBoundary(P_temp, t_n)
    P_temp(t_n + 1, 1, :)           = P_temp(t_n + 1, 2, :); % Top boundary
    P_temp(t_n + 1, end, :)         = P_temp(t_n + 1, end-1, :); % Bottom boundary
    P_temp(t_n + 1, :, 1)           = P_temp(t_n + 1, :, 2); % Left boundary
    P_temp(t_n + 1, :, end)         = P_temp(t_n + 1, :, end-1); % Right boundary
    P_bound = P_temp;
end
