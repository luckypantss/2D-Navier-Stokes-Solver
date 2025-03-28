function [U] = NSsolverTest2(X, Y, dx, dy, t, dt, rho, nu, F_y, beta, D_h, isPlot)
    format long
    % Initialize matrices
    U = zeros([size(X), 2]);
    P = zeros(size(X));
    U_temp = zeros(size(U));
    P_temp = zeros(size(P));

    % Variables for convective and diffusive terms
    u_conv = zeros(2, 1);
    u_diff = zeros(2, 1);

    % Variables for pressure gradient and velocity divergence
    p_force = 0;
    p_moment = 0;
    alpha = 0.15;

    % Variables for Reynolds number calculation
    Re = 100; 
    L = 1; 
    U_lid = 0.1;
    nu = U_lid*L/Re;



    % Initialize plots if required
    if isPlot
        figure(1);
        uQuiver = quiver(X, Y, U(:, :, 1), U(:, :, 2));
        title('Velocity Field [$\vec{U}$]', 'Interpreter', 'latex');
        xlabel('X');
        ylabel('Y');
        axis tight manual;

        figure(2);
        hContour = contourf(X, Y, P, 50, 'LineColor', 'none');
        colorbar;
        title('Pressure Field');

        figure(3);
        hStreamSlice = streamslice(X, Y, U(:, :, 1), U(:, :, 2), 'noarrows');
        title('Streamlines of the Flow');
        xlabel('X');
        ylabel('Y');
    end

    % Main loop for updating values
    for t_n = 1:t
        for i = 2:size(X, 1) - 1
            for j = 2:size(Y, 2) - 1
                for A = 1:2
                    u_conv(A, 1) = -0.5 * ( ...
                        U(i, j, 1) * (U(i+1, j, A) - U(i-1, j, A)) / dx + ...
                        U(i, j, 2) * (U(i, j+1, A) - U(i, j-1, A)) / dy ...
                    );
                    u_diff(A, 1) =  ( ...
                        (U(i+1, j, A) - 2 * U(i, j, A) + U(i-1, j, A)) / dx^2 + ...
                        (U(i, j+1, A) - 2 * U(i, j, A) + U(i, j-1, A)) / dy^2 ...
                    ); 
                end

                p_force = 0.5 * beta * ( ...
                    (U(i+1, j, 1) - U(i-1, j, 1)) / dx + ...
                    (U(i, j+1, 2) - U(i, j-1, 2)) / dy ...
                );
                p_moment = -0.5 * ( ...
                    U(i, j, 1) * (P(i+1, j) - P(i-1, j)) / dx + ...
                    U(i, j, 2) * (P(i, j+1) - P(i, j-1)) / dy ...
                );

                P_temp(i, j) = P(i, j) + dt * (p_force + p_moment);

                %Smudge the pressure field
                P_temp(i, j) = alpha * P_temp(i, j) + ((1-alpha)/4) * (P_temp(i-1, j) + P_temp(i+1, j) + P_temp(i, j-1) + P_temp(i, j+1));

                U_temp(i, j, 1) = U(i, j, 1) + dt * ( ...
                    u_conv(1, 1) + nu * u_diff(1, 1) - ...
                    1 / (2 * rho * dx) * (P(i+1, j) - P(i-1, j)) + ...
                    1 / rho * F_y ...
                );
                U_temp(i, j, 2) = U(i, j, 2) + dt * ( ...
                    u_conv(2, 1) + nu * u_diff(2, 1) - ...
                    1 / (2 * rho * dy) * (P(i, j+1) - P(i, j-1)) + ...
                    1 / rho * F_y ...
                );

            end
        end

        U_temp = applyNoSlipBoundary(U_temp);

        U = U_temp;
        P = P_temp;
        

        U_max = max(U(:,:,1),[],"all");
        V_max = max(U(:,:,2),[],"all");
        P_max = max(P,[],"all");

        disp(['t = ', num2str(t_n), ' Re = ', num2str(Re), ' U_max: ', num2str(U_max)...
            , ' V_max: ', num2str(V_max), ' P_max', num2str(P_max)]);

        if isPlot
            clf; % Clear the current figure window
            uData = U(:, :, 1);
            vData = U(:, :, 2);
            PData = P; 

            figure(1);
            uQuiver = quiver(X, Y, uData, vData);
            title('Velocity Field [$\vec{U}$]', 'Interpreter', 'latex');
            xlabel('X');
            ylabel('Y');
            drawnow;

            figure(2);
            hContour = contourf(X, Y, PData, 50, 'LineColor', 'none');
            colorbar;
            title('Pressure Field');
            drawnow;

            figure(3);
            streamslice(X, Y, uData, vData, 'noarrows');
            title('Streamlines of the Flow');
            xlabel('X');
            ylabel('Y');
            drawnow;
        end

        u_tot = 0;
    end

end

function [U_bound] = applyNoSlipBoundary(U_field)
    % X-Direction
    U_field(1, :, 1) = 0; % Top boundary
    U_field(end, :, 1) = 0.1; % Bottom boundary
    U_field(:, 1, 1) = 0; % Left boundary
    U_field(:, end, 1) = 0; % Right boundary
    % Y-Direction
    U_field(1, :, 2) = 0; % Top boundary
    U_field(end, :, 2) = 0; % Bottom boundary
    U_field(:, 1, 2) = 0; % Left boundary
    U_field(:, end, 2) = 0; % Right boundary
    U_bound = U_field;
end
