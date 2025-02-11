%%
% case 3: weak GFM + slow GFL case
%%
clear;clc;
Ts = 1e-5;
%% Line data Format (line)
% All values are given on the same system base MVA
% 1: From bus   2: To bus   3: Resistance (pu) 4: Reactance  (pu)
% 5: Charge     (pu)  6: Transformer Tap Amplitute 7: base MVA
% 8: Nomonal Voltage (KV) 
%   1   2    3        4       5     6    7    8  
line=[...
    1	2	0.0035	0.0411	0.6987	0	100	345
    1	39	0.001	0.025	0.75	0	100	345
    2	3	0.0013	0.0151	0.2572	0	100	345
    2	25	0.007	0.0086	0.146	0	100	345
    2	30	0	    0.0181	0	    0	100	345
    3	4	0.0013	0.0213	0.2214	0	100	345
    3	18	0.0011	0.0133	0.2138	0	100	345
    4	5	0.0008	0.0128	0.1342	0	100	345
    4	14	0.0008	0.0129	0.1382	0	100	345
    5	8	0.0008	0.0112	0.1476	0	100	345
    6	5	0.0002	0.0026	0.0434	0	100	345
    6	7	0.0006	0.0092	0.113	0	100	345
    6	11	0.0007	0.0082	0.1389	0	100	345
    7	8	0.0004	0.0046	0.078	0	100	345
    8	9	0.0023	0.0363	0.3804	0	100	345
    9	39	0.001	0.025	1.2	    0	100	345
    10	11	0.0004	0.0043	0.0729	0	100	345
    10	13	0.0004	0.0043	0.0729	0	100	345
    10	32	0	    0.02	0	    0	100	345
    12	11	0.0016	0.0435	0	    0	100	345
    12	13	0.0016	0.0435	0	    0	100	345
    13	14	0.0009	0.0101	0.1723	0	100	345
    14	15	0.0018	0.0217	0.366	0	100	345
    15	16	0.0009	0.0094	0.171	0	100	345
    16	17	0.0007	0.0089	0.1342	0	100	345
    16	19	0.0016	0.0195	0.304	0	100	345
    16	21	0.0008	0.0135	0.2548	0	100	345
    16	24	0.0003	0.0059	0.068	0	100	345
    17	18	0.0007	0.0082	0.1319	0	100	345
    17	27	0.0013	0.0173	0.3216	0	100	345
    19	33	0.0007	0.0142	0	    0	100	345
    19	20	0.0007	0.0138	0	    0	100	345
    20	34	0.0009	0.018	0	    0	100	345
    21	22	0.0008	0.014	0.2565	0	100	345
    22	23	0.0006	0.0096	0.1846	0	100	345
    22	35	0	    0.0143	0	    0	100	345
    23	24	0.0022	0.035	0.361	0	100	345
    23	36	0.0005	0.0272	0	    0	100	345
    25	26	0.0032	0.0323	0.513	0	100	345
    25	37	0.0006	0.0232	0	    0	100	345
    26	27	0.0014	0.0147	0.2396	0	100	345
    26	28	0.0043	0.0474	0.7802	0	100	345
    26	29	0.0057	0.0625	1.029	0	100	345
    28	29	0.0014	0.0151	0.249	0	100	345
    29	38	0.0008	0.0156	0	    0	100	345
    31	6	0	    0.025	0	    0	100	345];  
line_pu = line;
%% Bus data (Bus)
% 1. Bus number 2. Nominal phase-phase voltage 
%    1    2
Bus = [ ...
    1   345
    2   345
    3   345
    4   345
    5   345
    6   345
    7   345
    8   345
    9   345
    10   345
    11   345
    12   230
    13   345
    14   345
    15   345
    16   345
    17   345
    18   345
    19   345
    20   345
    21   345
    22   345
    23   345
    24   345
    25   345
    26   345
    27   345
    28   345
    29   345
    30    22
    31    22
    32    22
    33    22
    34    22
    35    22
    36    22
    37    22
    38   345
    39   345];
%% Bases and Actual Values
s=1;
zbase=(line(:,8).^2)./line(:,7);        % Zbase for each line (O)
fn = 60;                                % fbase (Hz)
wn = 2*pi*fn;                           % wbase (rad)
line(:,3)=line(:,3).*zbase;             % Actual R (O)
line(:,4)=line(:,4).*zbase/wn;          % Actual L (H)
line(:,5)=line(:,5)./zbase/wn;          % Actual C (F)
I_base = 100./Bus(:,2)*1000/sqrt(3);    % I base for each line (A)
%% Synchronous machine parameters
% Nominal power and frequency
num_gen = 10;                                                % number of generators
Pn_sm = 1000*ones(1,num_gen);                                % Power Bases [MW]
cosPhi_sm = ones(1,num_gen);                                 % nominal rated reactive power [MVar]
Un_sm = [345,22,22,22,22,22,22,22,345,22];                   % nominal rated voltage (line-to-line) [kV]
fn_sm = fn*ones(1,10);                                       % nominal rated frequency [Hz]
Kd_sm = 0*ones(1,10);                                        % damping factor [p.u.]
H_sm = [5,3.03,3.58,2.86,2.6,3.48,2.64,2.43,3.45,4.2];       % inertia coefficient [s]
wn_sm = 1;                                                   % omega set points [p.u.]
wb_sm = wn_sm*wn;                                            % omega set points (rad/s)
Pref_sm = [500 250 600 978 510 650 560 540 830 550]./Pn_sm; % Active Power Generation [p.u.]
Vref_sm = 1*ones(1,10);                                      % excitation voltage refrence [p.u.]
% Computed parameters
Sn_sm = Pn_sm./cosPhi_sm;           % nominal rated apparent power [MVA]
Qn_sm = sqrt(Sn_sm.^2-Pn_sm.^2);    % reactive power [MVA]
In_sm = Sn_sm./Un_sm;               % nominal rated current [kA]
Zb_sm = Un_sm.^2./Sn_sm;            % R/X Base [O]
Lb_sm = Zb_sm./(2*pi*fn_sm);        % L Base [H]
% Per unit reactances
xd = 1.305;                         % synchronous reactance in d-axis [p.u.]
xd_t = 0.296;                       % transient reactance in d-axis [p.u.]
xd_st = 0.252;                      % subtransient reactance in d-axis [p.u.]
xq = 0.474;                         % synchronous reactance in q-axis [p.u.]
xq_t = 0.5;                         % transient reactance in q-axis [p.u.]
xq_st = 0.243;                      % subtransient reactance in q-axis [p.u.]
xl = 0.18;                          % armature leakage reactance [p.u.]
% Time constants
Td_t = 1.01;                        % transient short-circuit time constant [s]
Td_st = 0.053;                      % subtransient short-circuit time constant [s]
Tq_t = 0.6;                         % transient short-circuit time constant [s]
Tq_st = 0.1;                        % subtransient short-circuit time constant [s]
% Other parameters
rs = 0.28;                          % stator resistance [p.u.]
p = 1;                              % number of pole pairs
xd = xd_t;
%% Governor
R = 0.05;            % controller droop gain [p.u.]
T1 = 0.1;            % governor time constant [sec]
T2 = 0.5;            % turbine derivative time constant [sec]
T3 = 0.5;            % turbine delay time constant [sec]
Tg = [10,10,10];
Kg = [0.95,0.95,0.95];
Fg = [0.1,0.1,0.1];
Rg = [0.05,0.05,0.05];    % Turbines
Dt = 0;                   % frictional losses factor
Vmin = 0;                 % minimum gate limit
Vmax = 2;                 % maximum gate limit
wref = 1;       
%% Automatic Voltage Regulator (AVR)
Ke = 1000;       % AVR controller gain [150-300]
Te = 0.05;      % exciter time constant [sec]
Tfa = 3;        % filter derivative time constant [sec]
Tfb = 10;       % filter delay time [sec]
Emin = 0.5;     % controller minimum output
Emax = 1.3;     % controller maximum output
Kas = 2.16e-4;  % AVR scaling gain
%% Power System Stabilizer (PSS)
Ts1 = 0.76; %0.76
Ts2 = 0.1; %20.1
Ts3 = 0.76; %0.76
Ts4 = 0.1; %0.1
Ks = 1; %1
Ts5 = 10; %10
Epss = 0.09;
%% Initial conditions for SMs
dw0_sm = 0;                      % frequency deviation [%]
theta0_sm = 0;                   % rotor angle [deg]
ia0_sm = 0;                      % stator currents in phase a [p.u.]
ib0_sm = 0;                      % stator currents in phase b [p.u.]
ic0_sm = 0;                      % stator currents in phase c [p.u.]
pha0_sm = 0;                     % angle of phase a [deg]
phb0_sm = 0;                     % angle of phase b [deg]
phc0_sm = 0;                     % angle of phase c [deg]
vf0_sm = 1.29;                   % excitation voltage [p.u.]
%% ********************************************************
%% ********************************************************
%% Base Values used for converters
Sb = 1000*10^6;                         %base apparent power [VA]
fb = fn;                               %base frequency [Hz]
w_b = fb*2*pi;                         %base angular frequency [rad/s]
wb = w_b;
load('PQmat.mat');
%% GFL Converter Parameters (p.u.)
w_pll = 10;                                            % bandwidth of PLL
ki_pll_pu = w_pll^2/(2+sqrt(5));                       % integral gain PLL
kp_pll_pu = sqrt(2*ki_pll_pu);                         % proportional gain PLL 
kp_i_pu = 0.3;                                         % proportional gain of PI
ki_i_pu = 10;                                          % integral gain of PI 
kp_q_pu = 0.1;                                         % proportional gain of PI
ki_q_pu = 1;                                           % integral gain of PI
kp_p_pu = 0.1;                                         % proportional gain of PI
ki_p_pu = 1;                                           % integral gain of PI
r_f_pu = 0.01;                                         % output filter resistor 
c_f_pu = 0.005;                                        % output filter capacitor 
l_f_pu = 0.9;                                          % output filter inductance
%% *********************** Converter 1 (GFL1) - Replace G1 **********************
p_set_pu(1) = 0.5012; v_set_pu(1) =1; q_set_pu(1) = 0.5291;
Vb(1) = 345*10^3;                                      % base high voltage, l-l rms [V]
[Vdc(1),pset(1),qset(1),kp_pll(1),ki_pll(1),...
    kp_i(1),ki_i(1),kp_q(1),ki_q(1),kp_p(1),ki_p(1),r_f(1),c_f(1),l_f(1),Init_GFL(1)]...
    = fun00_GFL_pu_2_normial(Sb,Vb(1),fb,p_set_pu(1),q_set_pu(1),kp_pll_pu,ki_pll_pu,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%
w_pll_s = 5;                                            % bandwidth of PLL
ki_pll_pu_s = w_pll_s^2/(2+sqrt(5));                     % integral gain PLL
kp_pll_pu_s = sqrt(2*ki_pll_pu_s);                       % proportional gain PLL 
[~,~,~,kp_pll_s(1),ki_pll_s(1),~,~,~,~,~,~,~,~,~,~]...
    = fun00_GFL_pu_2_normial(Sb,Vb(1),fb,p_set_pu(1),q_set_pu(1),kp_pll_pu_s,ki_pll_pu_s,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%% *********************** Converter 2 (GFL2) - Replace G2 **********************
p_set_pu(2) = 0.2505; v_set_pu(2) =1; q_set_pu(2) = 0.2846;
Vb(2) = 22*10^3;                                       % base high voltage, l-l rms [V]
[Vdc(2),pset(2),qset(2),kp_pll(2),ki_pll(2),...
    kp_i(2),ki_i(2),kp_q(2),ki_q(2),kp_p(2),ki_p(2),r_f(2),c_f(2),l_f(2),Init_GFL(2)]...
    = fun00_GFL_pu_2_normial(Sb,Vb(2),fb,p_set_pu(2),q_set_pu(2),kp_pll_pu,ki_pll_pu,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%
w_pll_s = 5;                                            % bandwidth of PLL
ki_pll_pu_s = w_pll_s^2/(2+sqrt(5));                     % integral gain PLL
kp_pll_pu_s = sqrt(2*ki_pll_pu_s);                       % proportional gain PLL 
[~,~,~,kp_pll_s(2),ki_pll_s(2),~,~,~,~,~,~,~,~,~,~]...
    = fun00_GFL_pu_2_normial(Sb,Vb(2),fb,p_set_pu(2),q_set_pu(2),kp_pll_pu_s,ki_pll_pu_s,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%% *********************** Converter 3 (GFL3)  - Replace G6 **********************
p_set_pu(3) = 0.6506; v_set_pu(3) =1; q_set_pu(3) = 0.2908;
Vb(3) = 22*10^3;                                       %base high voltage, l-l rms [V]
[Vdc(3),pset(3),qset(3),kp_pll(3),ki_pll(3),...
    kp_i(3),ki_i(3),kp_q(3),ki_q(3),kp_p(3),ki_p(3),r_f(3),c_f(3),l_f(3),Init_GFL(3)]...
    = fun00_GFL_pu_2_normial(Sb,Vb(3),fb,p_set_pu(3),q_set_pu(3),kp_pll_pu,ki_pll_pu,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%
w_pll_s = 5;                                            % bandwidth of PLL
ki_pll_pu_s = w_pll_s^2/(2+sqrt(5));                     % integral gain PLL
kp_pll_pu_s = sqrt(2*ki_pll_pu_s);                       % proportional gain PLL 
[~,~,~,kp_pll_s(3),ki_pll_s(3),~,~,~,~,~,~,~,~,~,~]...
    = fun00_GFL_pu_2_normial(Sb,Vb(3),fb,p_set_pu(3),q_set_pu(3),kp_pll_pu_s,ki_pll_pu_s,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%% *********************** Converter 4 (GFL4)  - Replace G8  **********************
p_set_pu(4) = 0.5407; v_set_pu(4) =1; q_set_pu(4) = 0.1341;
Vb(4) = 22*10^3;                                       % base high voltage, l-l rms [V]
[Vdc(4),pset(4),qset(4),kp_pll(4),ki_pll(4),...
    kp_i(4),ki_i(4),kp_q(4),ki_q(4),kp_p(4),ki_p(4),r_f(4),c_f(4),l_f(4),Init_GFL(4)]...
    = fun00_GFL_pu_2_normial(Sb,Vb(4),fb,p_set_pu(4),q_set_pu(4),kp_pll_pu,ki_pll_pu,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%
w_pll_s = 5;                                            % bandwidth of PLL
ki_pll_pu_s = w_pll_s^2/(2+sqrt(5));                     % integral gain PLL
kp_pll_pu_s = sqrt(2*ki_pll_pu_s);                       % proportional gain PLL 
[~,~,~,kp_pll_s(4),ki_pll_s(4),~,~,~,~,~,~,~,~,~,~]...
    = fun00_GFL_pu_2_normial(Sb,Vb(4),fb,p_set_pu(4),q_set_pu(4),kp_pll_pu_s,ki_pll_pu_s,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%% *********************** Converter 5 (GFL5)  - Replace G9  **********************
p_set_pu(5) = 0.6860; v_set_pu(5) =1; q_set_pu(5) = 0.1800;
Vb(5) = 345*10^3;                                       % base high voltage, l-l rms [V]
[Vdc(5),pset(5),qset(5),kp_pll(5),ki_pll(5),...
    kp_i(5),ki_i(5),kp_q(5),ki_q(5),kp_p(5),ki_p(5),r_f(5),c_f(5),l_f(5),Init_GFL(5)]...
    = fun00_GFL_pu_2_normial(Sb,Vb(5),fb,p_set_pu(5),q_set_pu(5),kp_pll_pu,ki_pll_pu,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%
w_pll_s = 5;                                            % bandwidth of PLL
ki_pll_pu_s = w_pll_s^2/(2+sqrt(5));                     % integral gain PLL
kp_pll_pu_s = sqrt(2*ki_pll_pu_s);                       % proportional gain PLL 
[~,~,~,kp_pll_s(5),ki_pll_s(5),~,~,~,~,~,~,~,~,~,~]...
    = fun00_GFL_pu_2_normial(Sb,Vb(5),fb,p_set_pu(5),q_set_pu(5),kp_pll_pu_s,ki_pll_pu_s,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%% *********************** Converter 6 (GFL6)  - Replace G5  **********************
p_set_pu(6) = 0.4460; v_set_pu(6) =1; q_set_pu(6) = 0.1560;
Vb(6) = 22*10^3;                                       % base high voltage, l-l rms [V]
[Vdc(6),pset(6),qset(6),kp_pll(6),ki_pll(6),...
    kp_i(6),ki_i(6),kp_q(6),ki_q(6),kp_p(6),ki_p(6),r_f(6),c_f(6),l_f(6),Init_GFL(6)]...
    = fun00_GFL_pu_2_normial(Sb,Vb(6),fb,p_set_pu(6),q_set_pu(6),kp_pll_pu,ki_pll_pu,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%
w_pll_s = 5;                                            % bandwidth of PLL
ki_pll_pu_s = w_pll_s^2/(2+sqrt(5));                     % integral gain PLL
kp_pll_pu_s = sqrt(2*ki_pll_pu_s);                       % proportional gain PLL 
[~,~,~,kp_pll_s(6),ki_pll_s(6),~,~,~,~,~,~,~,~,~,~]...
    = fun00_GFL_pu_2_normial(Sb,Vb(6),fb,p_set_pu(6),q_set_pu(6),kp_pll_pu_s,ki_pll_pu_s,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%% *********************** Converter 7 (GFL7)  - Replace G3  **********************
p_set_pu(7) = 0.6110; v_set_pu(7) =1; q_set_pu(7) = 0.3221;
Vb(7) = 22*10^3;                                       % base high voltage, l-l rms [V]
[Vdc(7),pset(7),qset(7),kp_pll(7),ki_pll(7),...
    kp_i(7),ki_i(7),kp_q(7),ki_q(7),kp_p(7),ki_p(7),r_f(7),c_f(7),l_f(7),Init_GFL(7)]...
    = fun00_GFL_pu_2_normial(Sb,Vb(7),fb,p_set_pu(7),q_set_pu(7),kp_pll_pu,ki_pll_pu,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%
w_pll_s = 5;                                            % bandwidth of PLL
ki_pll_pu_s = w_pll_s^2/(2+sqrt(5));                     % integral gain PLL
kp_pll_pu_s = sqrt(2*ki_pll_pu_s);                       % proportional gain PLL 
[~,~,~,kp_pll_s(7),ki_pll_s(7),~,~,~,~,~,~,~,~,~,~]...
    = fun00_GFL_pu_2_normial(Sb,Vb(7),fb,p_set_pu(7),q_set_pu(7),kp_pll_pu_s,ki_pll_pu_s,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f_pu,c_f_pu,l_f_pu);
%%
Init_GFL(1).theta0 = theta_0(1);
Init_GFL(2).theta0 = theta_0(2);
Init_GFL(3).theta0 = theta_0(4);
Init_GFL(4).theta0 = theta_0(6);
Init_GFL(5).theta0 = 4;
Init_GFL(6).theta0 = 4;
Init_GFL(7).theta0 = 4;
%% GFM Converter Parameters (p.u.)
Vb_gfm = 22*10^3;                         %base high voltage, l-l rms [V]
Vb_p_ph = Vb_gfm*sqrt(2/3);               %Base low voltage, p-n peak [V]
% VSM
J_pu = 5;                                %inertia 
D_p_pu = 100;                             %active power damping
D_q_pu = 0.01;                            %reactive power damping
kp_vc_pu = 0.01;
ki_vc_pu = 0.1;
kp_i_pu = 0.3;                            %proportional gain of PI
ki_i_pu = 10;                             %integral gain of PI 
r_f_pu = 0.01;
c_f_pu = 0.005;                            
l_f_pu = 0.9;  
%% *********************** Converter x (GFM1) - Replace G10 **********************
P_set_pu_gfm(1) = 0.5610; Q_set_pu_gfm(1) = 0.2628;                      %Set Points [pu]
[Vdc_gfm(1),P_set_gfm(1),Q_set_gfm(1),J_gfm(1),D_p_gfm(1),D_q_gfm(1),...
    kp_vc_gfm(1),ki_vc_gfm(1),kp_i_gfm(1),ki_i_gfm(1),r_f_gfm(1),c_f_gfm(1),l_f_gfm(1),Init(1)]...
    = fun00_GFM_pu_2_normial(Sb,Vb_gfm,fb,P_set_pu_gfm(1),Q_set_pu_gfm(1),J_pu,...
    D_p_pu,D_q_pu,kp_vc_pu,ki_vc_pu,kp_i_pu,ki_i_pu,r_f_pu,l_f_pu,c_f_pu);
%
D_p_pu_s = 35;                             %active power damping
[~,~,~,~,D_p_gfm_s(1),~,~,~,~,~,~,~,~,~]...
    = fun00_GFM_pu_2_normial(Sb,Vb_gfm,fb,P_set_pu_gfm(1),Q_set_pu_gfm(1),J_pu,...
    D_p_pu_s,D_q_pu,kp_vc_pu,ki_vc_pu,kp_i_pu,ki_i_pu,r_f_pu,l_f_pu,c_f_pu);
%% *********************** Converter x (GFM2) - Replace G7 ********************** 
P_set_pu_gfm(2) = 0.5760; Q_set_pu_gfm(2) = 0.1883;                      %Set Points [pu]
[Vdc_gfm(2),P_set_gfm(2),Q_set_gfm(2),J_gfm(2),D_p_gfm(2),D_q_gfm(2),...
    kp_vc_gfm(2),ki_vc_gfm(2),kp_i_gfm(2),ki_i_gfm(2),r_f_gfm(2),c_f_gfm(2),l_f_gfm(2),Init(2)]...
    = fun00_GFM_pu_2_normial(Sb,Vb_gfm,fb,P_set_pu_gfm(2),Q_set_pu_gfm(2),J_pu,...
    D_p_pu,D_q_pu,kp_vc_pu,ki_vc_pu,kp_i_pu,ki_i_pu,r_f_pu,l_f_pu,c_f_pu);
%
D_p_pu_s = 35;                             %active power damping
[~,~,~,~,D_p_gfm_s(2),~,~,~,~,~,~,~,~,~]...
    = fun00_GFM_pu_2_normial(Sb,Vb_gfm,fb,P_set_pu_gfm(2),Q_set_pu_gfm(2),J_pu,...
    D_p_pu_s,D_q_pu,kp_vc_pu,ki_vc_pu,kp_i_pu,ki_i_pu,r_f_pu,l_f_pu,c_f_pu);
%%
Init(1).theta0 = theta_0(7);
Init(2).theta0 = theta_0(5);
Init(3).theta0 = theta_0(3);
%%
t_s = 70;
t_gfl = 50;
t_gfm = 25;
sim('IEEE_39_case3.slx',[0 t_s]);        