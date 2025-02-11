clear; clc;
%% Network Base Values 
Sb = 50*10^3;                          % Base apparent power [VA]
Vb = 690;                              % Base high voltage, line-to-line RMS [V]
fb = 60;                               % Base frequency [Hz]
w_b = fb*2*pi;                         % Base angular frequency [rad/s]
%% Calculating Base Quantities
Vb_p_ph = Vb*sqrt(2/3);                % Base low voltage, phase-to-neutral peak [V]
Ib_p_ph = Sb/Vb/sqrt(3)*sqrt(2);       % Base peak current for converter, phase-to-neutral [A]
Zb = Vb^2/Sb;                          % Base impedance for converter [Ohm]
Lb = Zb/w_b;                           % Base inductance [H]
Cb = 1/(Zb*w_b);                       % Base capacitance [F]
%% Setting points for power and reactive power (in per-unit system)
p_set_pu = 1;                          % Active power setpoint [pu]
q_set_pu = 1.135*10^4/Sb;              % Reactive power setpoint [pu]
%% Converter Parameters (non-per-unit values)
% Phase-Locked Loop (PLL) Controller
kp_pll = 27.5;                         % Proportional gain for PLL
ki_pll = 377.7;                        % Integral gain for PLL
% Inner Current Control Loop
kp_i_pu = 0.3;                         % Proportional gain for current PI controller
ki_i_pu = 10;                          % Integral gain for current PI controller
% Outer Power and Voltage PI-Droop Controller
kp_q_pu = 0.1;                         % Proportional gain for reactive power PI controller
ki_q_pu = 5;                           % Integral gain for reactive power PI controller
kp_p_pu = 0.1;                         % Proportional gain for active power PI controller
ki_p_pu = 5;                           % Integral gain for active power PI controller
% Output Filter Parameters
r_f = 0.01;                            % Output filter resistance [Ohm]
c_f = 0.06;                            % Output filter capacitance [F]
l_f = 0.05;                            % Output filter inductance [H]
%% Convert per-unit values to nominal values
[Vdc,pset,qset,kp_pll,ki_pll,...
    kp_i,ki_i,kp_q,ki_q,kp_p,ki_p,r_f,c_f,l_f]...
    = fun00_GFL_pu_2_normial(Sb,Vb,fb,p_set_pu,q_set_pu,kp_pll,ki_pll,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f,c_f,l_f);
%% Initialization
ln = logspace(0,3,10);                 % Generate logarithmic spaced frequency points
nn = size(ln);                         % Get the size of ln
nn = max(nn);                          % Get the number of columns in ln
freq = zeros(1,nn);                    % Initialize frequency array
%% Frequency Response Analysis
for ii = 1:nn
    tptf = 5;                          % Time to reach steady-state before analysis [s]
    Upd = 10;                          % d-axis disturbance voltage magnitude (<= 5%)
    Upq = 0;                           % q-axis disturbance voltage magnitude (<= 5%)
    freq(ii) = fix(ln(ii));            % Round harmonic frequency to nearest integer
    fp1 = freq(ii);                    % Set fundamental perturbation frequency
    fpt = fp1;                         % Assign frequency for perturbation
    sim('GFL_mask',[0 15]);            % Run simulation from 0 to 15 seconds
    % Extract Vd component
    FFTDATA = power_fftscope(Vdq);
    FFTDATA.startTime = tptf;          % Start Fourier analysis at tptf seconds
    FFTDATA.fundamental = 1;           % Set fundamental frequency for analysis
    FFTDATA.maxFrequency = 1200;       % Maximum frequency for analysis [Hz]
    FFTDATA.signal = 1;                % Select d-axis voltage signal
    FFTDATA = power_fftscope(FFTDATA);
    Vmd = FFTDATA.mag;                 % Extract magnitude
    Vphd = FFTDATA.phase - 90;         % Extract phase angle (adjust for sine-to-cosine shift)
    % Extract Vq component
    FFTDATA=power_fftscope(Vdq);
    FFTDATA.startTime=tptf;
    FFTDATA.fundamental=1;
    FFTDATA.maxFrequency=1200; 
    FFTDATA.signal=2; 
    FFTDATA=power_fftscope(FFTDATA);
    Vmq=FFTDATA.mag;
    Vphq=FFTDATA.phase-90;
    % Extract Id component
    FFTDATA=power_fftscope(Idq);
    FFTDATA.startTime=tptf;
    FFTDATA.fundamental=1;
    FFTDATA.maxFrequency=1200; 
    FFTDATA.signal=1; 
    FFTDATA=power_fftscope(FFTDATA);
    Imd = FFTDATA.mag;
    Iphd = FFTDATA.phase - 90;
    % Extract Iq component
    FFTDATA=power_fftscope(Idq);
    FFTDATA.startTime=tptf;
    FFTDATA.fundamental=1;
    FFTDATA.maxFrequency=1200; 
    FFTDATA.signal=2; 
    FFTDATA=power_fftscope(FFTDATA);
    Imq=FFTDATA.mag;
    Iphq=FFTDATA.phase-90;
    
    % Compute phasor values for fundamental frequency
    Vd1(ii) = Vmd(fp1+1)*cos(Vphd(fp1+1)*pi/180) + 1i*Vmd(fp1+1)*sin(Vphd(fp1+1)*pi/180);
    Id1(ii) = Imd(fp1+1)*cos(Iphd(fp1+1)*pi/180) + 1i*Imd(fp1+1)*sin(Iphd(fp1+1)*pi/180);
    
    Vq1(ii) = Vmq(fp1+1)*cos(Vphq(fp1+1)*pi/180) + 1i*Vmq(fp1+1)*sin(Vphq(fp1+1)*pi/180);
    Iq1(ii) = Imq(fp1+1)*cos(Iphq(fp1+1)*pi/180) + 1i*Imq(fp1+1)*sin(Iphq(fp1+1)*pi/180);

    %% Second Disturbance Injection (q-axis)
    Upd = 0;                           % Ensure disturbances are linearly independent
    Upq = 10;                          % Apply disturbance in q-axis
    
    sim('GFL_mask',[0 15]);            % Run simulation again
    
    % Repeat FFT extraction for second perturbation
    FFTDATA=power_fftscope(Vdq);
    FFTDATA.startTime=tptf;
    FFTDATA.fundamental=1;
    FFTDATA.maxFrequency=1200;
    FFTDATA.signal=1; 
    FFTDATA=power_fftscope(FFTDATA);
    Vmd=FFTDATA.mag;
    Vphd=FFTDATA.phase-90;
    
    FFTDATA=power_fftscope(Vdq);
    FFTDATA.startTime=tptf;
    FFTDATA.fundamental=1;
    FFTDATA.maxFrequency=1200;
    FFTDATA.signal=2; 
    FFTDATA=power_fftscope(FFTDATA);
    Vmq=FFTDATA.mag;
    Vphq=FFTDATA.phase-90;

    FFTDATA=power_fftscope(Idq);
    FFTDATA.startTime=tptf;
    FFTDATA.fundamental=1;
    FFTDATA.maxFrequency=1200;
    FFTDATA.signal=1; 
    FFTDATA=power_fftscope(FFTDATA);
    Imd=FFTDATA.mag;
    Iphd=FFTDATA.phase-90;

    FFTDATA=power_fftscope(Idq);
    FFTDATA.startTime=tptf;
    FFTDATA.fundamental=1;
    FFTDATA.maxFrequency=1200;
    FFTDATA.signal=2; 
    FFTDATA=power_fftscope(FFTDATA);
    Imq=FFTDATA.mag;
    Iphq=FFTDATA.phase-90;
    
    Vd2(ii)=Vmd(fp1+1)*cos(Vphd(fp1+1)*pi/180)+1i*Vmd(fp1+1)*sin(Vphd(fp1+1)*pi/180);
    Id2(ii)=Imd(fp1+1)*cos(Iphd(fp1+1)*pi/180)+1i*Imd(fp1+1)*sin(Iphd(fp1+1)*pi/180);

    Vq2(ii)=Vmq(fp1+1)*cos(Vphq(fp1+1)*pi/180)+1i*Vmq(fp1+1)*sin(Vphq(fp1+1)*pi/180);
    Iq2(ii)=Imq(fp1+1)*cos(Iphq(fp1+1)*pi/180)+1i*Imq(fp1+1)*sin(Iphq(fp1+1)*pi/180);

    %% Compute Admittance Matrix
    Vdq2(:,:,ii) = [Vd1(ii),Vd2(ii);Vq1(ii),Vq2(ii)];
    Idq2(:,:,ii) = [Id1(ii),Id2(ii);Iq1(ii),Iq2(ii)];
    
    Ydq(:,:,ii) = Idq2(:,:,ii) * inv(Vdq2(:,:,ii));

    Yddm(ii)=abs(Ydq(1,1,ii));
    Yddph(ii)=angle(Ydq(1,1,ii))*180/pi;

    Ydqm(ii)=abs(Ydq(1,2,ii));
    Ydqph(ii)=angle(Ydq(1,2,ii))*180/pi;

    Yqdm(ii)=abs(Ydq(2,1,ii));
    Yqdph(ii)=angle(Ydq(2,1,ii))*180/pi;

    Yqqm(ii)=abs(Ydq(2,2,ii));
    Yqqph(ii)=angle(Ydq(2,2,ii))*180/pi;
end
%%
save('Ydq.mat', 'Ydq');  % Save admittance matrix data
%% Plot admittance characteristics for different frequency responses
figure
subplot(2,1,1); 
semilogx(freq,Yddm(1:nn));
subplot(2,1,2); 
semilogx(freq,Yddph(1:nn));
sgtitle('Ydd scan')

figure
subplot(2,1,1);
semilogx(freq,Ydqm(1:nn));
subplot(2,1,2);
semilogx(freq,Ydqph(1:nn));
sgtitle('Ydq scan')

figure
subplot(2,1,1);
semilogx(freq,Yqdm(1:nn));
subplot(2,1,2);
semilogx(freq,Yqdph(1:nn));
sgtitle('Yqd scan')

figure
subplot(2,1,1);
semilogx(freq,Yqqm(1:nn));
subplot(2,1,2);
semilogx(freq,Yqqph(1:nn));
sgtitle('Yqq scan')
