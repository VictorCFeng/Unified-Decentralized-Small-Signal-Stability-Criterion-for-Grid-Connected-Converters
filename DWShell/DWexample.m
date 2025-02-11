clear; clc;  % Clear workspace and command window
% IMPORTANT: The function `bus_calculation` uses parallel computing (`parfor`).
% If Parallel Computing Toolbox is NOT available, modify `bus_calculation.m`:
% - Replace `parfor` with `for` to ensure compatibility.
%% Define matrix A
A = [0.58-0.21i  -0.92+0.41i  0.35-0.90i;
     0.91+0.31i   0.69-0.93i  0.51-0.80i;
     0.31-0.65i   0.86-0.44i  0.48+0.64i];
%% Call the bus_calculation function
% Compute the numerical range and Davis-Wielandt shell for -A^(-1)
GainPhase = bus_calculation(A, 1); 
%% Adjust plot settings
h.XAxis.Label.Rotation = -13;  % Rotate X-axis label for better visibility
h.YAxis.Label.Rotation = 60;   % Rotate Y-axis label for better visibility
view(20, 30);                  % Set 3D view angle (Azimuth=20°, Elevation=30°)
grid minor;                     % Enable minor grid for better visualization

