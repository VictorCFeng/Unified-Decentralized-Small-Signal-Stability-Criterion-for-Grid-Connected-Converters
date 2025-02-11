%%
clear; clc;
Ts = 1e-5;
%% Network Base Values 
Sb = 50*10^3;                          %base apparent power [VA]
Vb = 690;                              %base high voltage, l-l rms [V]
fb = 60;                               %base frequency [Hz]
w_b = fb*2*pi;                         %base angular frequency [rad/s]
%% 
Vb_p_ph = Vb*sqrt(2/3);                %Base low voltage, p-n peak [V]
Ib_p_ph = Sb/Vb/sqrt(3)*sqrt(2);       %L-n peak Converter base current [A]
Zb = Vb^2/Sb;                          %Converter base impedance [Ohm]
Lb = Zb/w_b;                           %base inductance [H]
Cb = 1/(Zb*w_b);                       %base capacitance [F]
%% setting points
P_set_pu = 1;                          % Active power setpoint [pu]
Q_set_pu = 1.34*10^4/Sb;               % Reactive power setpoint [pu]
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
[Vdc,P_set,Q_set,J,D_p,D_q,...
    kp_vc,ki_vc,kp_i,ki_i,r_f,c_f,l_f]...
    = fun00_GFM_pu_2_normial(Sb,Vb,fb,P_set_pu,Q_set_pu,J_pu,...
    D_p_pu,D_q_pu,kp_vc_pu,ki_vc_pu,kp_i_pu,ki_i_pu,r_f,l_f,c_f);
%% Initialization
ln = logspace(0,3,10);                 % Generate logarithmic spaced frequency points
nn = size(ln);                         % Get the size of ln
nn = max(nn);                          % Get the number of columns in ln
freq = zeros(1,nn);                    % Initialize frequency array
%%
for ii=1:nn
    tptf = 5;                          % Time to reach steady-state before analysis [s]
    Upd = 10;                          % d-axis disturbance voltage magnitude (<= 5%)
    Upq = 0;                           % q-axis disturbance voltage magnitude (<= 5%)
    freq(ii) = fix(ln(ii));            % Round harmonic frequency to nearest integer
    fp1 = freq(ii);                    % Set fundamental perturbation frequency
    fpt = fp1;                         % Assign frequency for perturbation
    sim('GFM_mask',[0 10]);            % Run simulation from 0 to 15 seconds
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
    
    sim('GFM_mask',[0 10]);            % Run simulation again
    
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
