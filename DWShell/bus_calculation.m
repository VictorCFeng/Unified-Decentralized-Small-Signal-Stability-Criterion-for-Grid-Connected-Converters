function gf = bus_calculation(A, plt_flag)
% This function calculates the numerical range and Davis-Wielandt shell 
% of a given matrix A and optionally visualizes the results.
% 
% Inputs:
%   A - Input matrix
%   plt_flag - Flag for plotting (1 to enable visualization, 0 otherwise)
%
% Output:
%   gf - Structure containing calculated properties

    n = size(A, 1);  % Get the size of matrix A (assumed to be square)
    num_samples = 300000;  % Number of random samples for visualization
    P = 500;  % Number of sampling points for the boundary calculation

    % Compute extreme points of the numerical range
    count = 1;
    [eigvec, ~] = eig(A);  % Compute eigenvectors of A
    for j = 1:n
        u = eigvec(:, j) / norm(eigvec(:, j));  % Normalize eigenvector
        x(count) = real(u' * A * u);  % Real part of u^* A u
        y(count) = imag(u' * A * u);  % Imaginary part of u^* A u
        z(count) = real(u' * (A' * A) * u);  % Real part of u^* A^* A u
        count = count + 1;
    end  

    % Compute additional boundary points
    count = 1;
    for theta = 0:P
        % Rotate A by a unit complex number on the unit circle
        T = (cos(theta * 2 * pi / P) + 1i * sin(theta * 2 * pi / P)) * A;
        ReT = 0.5 * (T + T');  % Symmetric part of T
        lambda = max(eig(ReT));  % Maximum eigenvalue of Re(T)
        [eigvec, eigval] = eig(ReT);  % Compute eigenvalues and eigenvectors of Re(T)

        for j = 1:n
            value = eigval(j, j);
            if value >= lambda
                u = eigvec(:, j) / norm(eigvec(:, j));  % Normalize eigenvector
                x1(count) = real(u' * A * u);
                y1(count) = imag(u' * A * u);
                z1(count) = real(u' * (A' * A) * u);
                count = count + 1;
            end
        end
    end

    % Combine computed points
    x = [x, x1]; 
    y = [y, y1]; 
    z = [z, z1];

    % Compute smoothed boundary
    shrink_factor = 0.01;  % Controls the smoothness of the boundary
    K1 = boundary([x', y'], shrink_factor);  % Compute 2D boundary

    % Compute singular values of A
    gf.min_sigma = min(svd(A));  % Minimum singular value
    gf.max_sigma = max(svd(A));  % Maximum singular value

    % Check if the origin (0,0) is inside the computed boundary
    flag = inpolygon(0, 0, x, y);
    if flag == 0
        % Compute angular range of the numerical range boundary
        angles = atan2d(y, x);
        gf.min_theta = min(angles);
        gf.max_theta = max(angles);
        if gf.max_theta - gf.min_theta >= 180
            angles(angles >= 0) = angles(angles >= 0) - 360;
            gf.min_theta = min(angles);
            gf.max_theta = max(angles);       
        end
    else
        % If the origin is inside, set the angle range to full circle
        gf.min_theta = -180;
        gf.max_theta = 180; 
    end

    % Visualization if plt_flag is set
    if plt_flag == 1
        % Generate random samples on the complex unit sphere
        parfor i = 1:num_samples
            pts = zeros(1, 3);
            rng(i);  % Set random seed for reproducibility
            u = randn(n, 1) + 1i * randn(n, 1);  % Generate random complex vector
            u = u / norm(u);  % Normalize to unit norm
            uAu = u' * A * u;
            pts(1) = real(uAu);  % Real part of u^* A u
            pts(2) = imag(uAu);  % Imaginary part of u^* A u
            pts(3) = real(u' * (A' * A) * u);  % Real part of u^* A^* A u
            points(i, :) = pts;
        end

        % Combine sampled points with computed boundary points
        points = [points; [x', y', z']];
        K = boundary(points, shrink_factor);  % Compute 3D boundary

        % Plot the Davis-Wielandt shell
        trisurf(K, points(:, 1), points(:, 2), points(:, 3), ...
                'FaceColor', 'interp', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on

        % Plot 2D numerical range boundary
        patch(x(K1), y(K1), zeros(size(K1, 1), 1), ...
              'FaceColor', 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on

        % Label axes using LaTeX notation
        xlabel('$\mathrm{Re}(u^*\mathbf{A}u)$', 'Interpreter', 'Latex');
        ylabel('$\mathrm{Im}(u^*\mathbf{A}u)$', 'Interpreter', 'Latex');
        zlabel('$u^*\mathbf{A}^*\mathbf{A}u$', 'Interpreter', 'Latex');

        % Improve lighting for better visualization
        camlight;
        lighting gouraud;
    end
end
