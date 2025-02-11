clear;clc;
%% Network Base Values 
Sb = 50*10^3;                          % Base apparent power [VA]
Vb = 690;                              % Base high voltage, line-to-line RMS [V]
fb = 60;                               % Base frequency [Hz]
w_b = fb*2*pi;                         % Base angular frequency [rad/s]
%% Calculating Base Quantities
Vb_p_ph = Vb*sqrt(2/3);                % Base low voltage, phase-to-neutral peak [V]
Ib_p_ph = Sb/Vb/sqrt(3)*sqrt(2);       % Base peak current for converter, phase-to-neutral [A]
Zb = Vb^2/Sb;                          % Base impedance for converter [Ohm]
Yb = 1/Zb;
Lb = Zb/w_b;                           % Base inductance [H]
Cb = 1/(Zb*w_b);                       % Base capacitance [F]
%% setting points
P_set_pu = 1;                          % Active power setpoint [pu]
Q_set_pu = 1.34*10^4/Sb;               % Reactive power setpoint [pu]
theta_0 = mod(31.93,2*pi);             % Power angle [rad]
%% Converter Parameters (All values in per unit (pu))
J_pu = 10 / (Sb / w_b);                % Inertia constant, affecting frequency dynamics
D_p_pu = 1270 / (Sb / w_b);            % Damping for active power control
D_q_pu = 500 / (Sb / Vb_p_ph);         % Damping for reactive power control
kp_vc_pu = 2;                          % Proportional gain for voltage control
ki_vc_pu = 10;                         % Integral gain for voltage control
kp_i_pu = 0.3;                         % Proportional gain for current control
ki_i_pu = 10;                          % Integral gain for current control
r_f = 0.01;                            % Filter resistance [pu]
c_f = 0.06;                            % Output filter capacitor [pu]
l_f = 0.05;                            % Output filter inductance [pu]
%%
YGFL_num = fun00_GFM_Small_Signal(Sb,Vb,fb,P_set_pu,Q_set_pu,J_pu,...
    D_p_pu,D_q_pu,kp_vc_pu,ki_vc_pu,kp_i_pu,ki_i_pu,r_f,l_f,c_f,theta_0);
%%
% Define frequency range for scanning
ScanFt = logspace(0, 3, 10000);
% Preallocate arrays for impedance calculations
Yddtest = zeros(1, length(ScanFt));
Ydqtest = zeros(1, length(ScanFt));
Yqdtest = zeros(1, length(ScanFt));
Yqqtest = zeros(1, length(ScanFt));
% Compute impedances across frequency range
for ff = 1:length(ScanFt)
    omega = 1i * 2 * pi * ScanFt(ff);
    YdqGFLScanFt = YGFL_num(omega);
    Yddtest(ff) = YdqGFLScanFt(1,1);
    Ydqtest(ff) = YdqGFLScanFt(1,2);
    Yqdtest(ff) = YdqGFLScanFt(2,1);
    Yqqtest(ff) = YdqGFLScanFt(2,2);
end
% Load measured impedance data
load Ydq.mat;
freq = fix(logspace(0, 3, 10));
% Extract real impedance data
Yddreal = squeeze(Ydq(1,1,:));
Ydqreal = squeeze(Ydq(1,2,:));
Yqdreal = squeeze(Ydq(2,1,:));
Yqqreal = squeeze(Ydq(2,2,:));
% Define impedance names and corresponding data
ImpedanceNames = {'Ydd', 'Ydq', 'Yqd', 'Yqq'};
TestData = {Yddtest, Ydqtest, Yqdtest, Yqqtest};
RealData = {Yddreal, Ydqreal, Yqdreal, Yqqreal};
% Create figure with 4x2 subplot layout
figure;
for i = 1:length(ImpedanceNames)
    % Get impedance name and corresponding data
    currentName = ImpedanceNames{i};
    testZ = TestData{i};
    realZ = RealData{i};
    
    % Plot magnitude response
    subplot(4, 2, 2*i-1);
    semilogx(ScanFt, abs(testZ/Yb), 'LineWidth', 1.5);
    hold on;
    semilogx(freq, abs(realZ/Yb), 'LineWidth', 1.5, 'Marker', '*', 'LineStyle', 'none');
    hold off;
    grid on;
    set(gca, 'FontSize', 8, 'GridLineStyle', ':', 'GridColor', '#000000', 'GridAlpha', 0.8);
    xlabel('Frequency [Hz]', 'FontName', 'Times New Roman');
    ylabel('Magnitude [abs]', 'FontName', 'Times New Roman');
    legend('Test', 'Real', 'Location', 'best');
    title(sprintf('%s Magnitude Comparison', currentName), 'FontSize', 10, 'FontName', 'Times New Roman');
    
    % Plot phase response
    subplot(4, 2, 2*i);
    phase_noisy = mod(180/pi * angle(testZ), 360);
    phase_real = mod(180/pi * angle(realZ), 360);
    
    semilogx(ScanFt, phase_noisy, 'LineWidth', 1.5);
    hold on;
    semilogx(freq, phase_real, 'LineWidth', 1.5, 'Marker', '*', 'LineStyle', 'none');
    hold off;
    grid on;
    set(gca, 'FontSize', 8, 'GridLineStyle', ':', 'GridColor', '#000000', 'GridAlpha', 0.8);
    xlabel('Frequency [Hz]', 'FontName', 'Times New Roman');
    ylabel('Phase [deg]', 'FontName', 'Times New Roman');
    legend('Test', 'Real', 'Location', 'best');
    title(sprintf('%s Phase Comparison', currentName), 'FontSize', 10, 'FontName', 'Times New Roman');
end

% Adjust figure properties
set(gcf, 'WindowStyle', 'normal');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 30, 40]); % Adjust figure size
set(gcf, 'Color', 'w'); % Set background color to white
