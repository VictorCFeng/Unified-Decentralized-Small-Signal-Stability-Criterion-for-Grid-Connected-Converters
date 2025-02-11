function YGFL_num...
    = fun01_GFL_Small_Signal_Q(Sb,Vb,fb,p_set_pu,q_set_pu,kp_pll,ki_pll,...
    kp_i_pu,ki_i_pu,kp_q_pu,ki_q_pu,kp_p_pu,ki_p_pu,r_f,c_f,l_f,theta_0)

%%
w_b = fb*2*pi;                         %base angular frequency [rad/s]
%% 
Vb_p_ph = Vb*sqrt(2/3);                %Base low voltage, p-n peak [V]
Ib_p_ph = Sb/Vb/sqrt(3)*sqrt(2);       %L-n peak Converter base current [A]
Zb = Vb^2/Sb;                          %Converter base impedance [Ohm]
Lb = Zb/w_b;                           %base inductance [H]
Cb = 1/(Zb*w_b);                       %base capacitance [F]
%% setting points
pset = p_set_pu*Sb; qset = q_set_pu*Sb;   %Set points [W] and [V]
%% Converter Parameters (already transformed into no-pu values)
%PLL
kp_pll = kp_pll/Vb_p_ph;                              %proportional gain PLL 
ki_pll = ki_pll/Vb_p_ph;                             %integral gain PLL
%inner current control loop 
kp_i = kp_i_pu*Vb_p_ph/Ib_p_ph;                           %proportional gain of PI
ki_i = ki_i_pu*Vb_p_ph/Ib_p_ph;                            %integral gain of PI 
%outer power and voltage PI-droop controller
kp_q = kp_q_pu/Sb*Ib_p_ph;                             %proportional gain of PI
ki_q = ki_q_pu/Sb*Ib_p_ph;                            %integral gain of PI
kp_p = kp_p_pu/Sb*Ib_p_ph;                                %proportional gain of PI
ki_p = ki_p_pu/Sb*Ib_p_ph;                                 %integral gain of PI
r_f = r_f*Zb;                                           %output filter resistor 
c_f = c_f*Cb;                                           %output filter capacitor 
l_f = l_f*Lb;                                          %output filter inductance
%% P = 3/2*(Vod*Iod+Voq*Ioq), Q = 3/2*(-Vod*Ioq+Voq*Iod), Voq =0
K1 = 0;
Vod = Vb_p_ph; Voq =0;
%
b = [pset;qset]; Func = 3/2*[Vod,Voq;Voq,-Vod];
x = inv(Func)*b;
Iod = x(1); Ioq = x(2);
%
LCL2 = [0 -w_b*K1*l_f;w_b*K1*l_f 0];
x = [Vod;Voq]+LCL2*[Iod;Ioq];
Vod1 = x(1); Voq1 = x(2); 
%
%
Func = [1           0           -r_f          w_b*l_f;
     0           1           -w_b*l_f     -r_f;
     0           0           1             0;
     0           0           0             1];
b = [1 0 0 0;
    0 1 0 0;
    0 -w_b*c_f 1 0;
    w_b*c_f 0 0 1]*[Vod1; Voq1; Iod; Ioq];
x = Func \ b;
Vcd = x(1); Vcq = x(2); Icd = x(3); Icq = x(4);
I = [1 0;0 1];
%%
syms s;
% 定义方程
PI_PLL = kp_pll + ki_pll/s; % PLL%
GPLL = PI_PLL/(s+Vod*PI_PLL);%
%%转角
Gvodq_csv = I+[0 Voq*GPLL;0 (-1)*Vod*GPLL];%
Gvcdq_csv = [0 Vcq*GPLL;0 (-1)*Vcd*GPLL];%
Giodq_csv = [0 Ioq*GPLL;0 (-1)*Iod*GPLL];%
Gicdq_csv = [0 Icq*GPLL;0 (-1)*Icd*GPLL];%
kdq = 0.5;
PI_CC = kp_i + ki_i/s; % Inner PI%
GPI_CC = [PI_CC 0;0 PI_CC];

PI_PC = kp_p + ki_p/s; % PC%
PI_VC = kp_q + ki_q/s; % VC%

GRL_f = [r_f -w_b*l_f;w_b*l_f r_f];
Z1 = [r_f+s*l_f -w_b*l_f;w_b*l_f r_f+s*l_f];
YC = [s*c_f -w_b*c_f;w_b*c_f s*c_f];

A = [-3/2*PI_PC*Vod -3/2*PI_PC*Voq;3/2*PI_VC*Voq -3/2*PI_VC*Vod];
B = [-3/2*PI_PC*Iod -3/2*PI_PC*Ioq;-3/2*PI_VC*Ioq 3/2*PI_VC*Iod];

Gvcdq_svo = (-GPI_CC*I+kdq*GRL_f)*Gicdq_csv + ...
    (GPI_CC*B+kdq*I)*Gvodq_csv + ...
    GPI_CC*A*Giodq_csv - Gvcdq_csv; 
Gvcdq_sic = GPI_CC*(A-I)+kdq*GRL_f;

% 定义旋转矩阵
e_Jdelta = [cos(theta_0), -sin(theta_0); sin(theta_0), cos(theta_0)];
Y = (Gvcdq_sic*I-Z1)^(-1)*((I+Z1*YC)-Gvcdq_svo-Gvcdq_sic*YC);
Y = -Y;

YGFL = e_Jdelta*Y*(e_Jdelta)^(-1);

YGFL_num = matlabFunction(YGFL);
end