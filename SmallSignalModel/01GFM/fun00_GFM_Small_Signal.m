function YGFL_num = fun00_GFM_Small_Signal(Sb,Vb,fb,P_set_pu,Q_set_pu,J_pu,...
    D_p_pu,D_q_pu,kp_vc_pu,ki_vc_pu,kp_i_pu,ki_i_pu,r_f,l_f,c_f,theta_0)

%%
w_b = fb*2*pi;                         %base angular frequency [rad/s]
%% 
Vb_p_ph = Vb*sqrt(2/3);                %Base low voltage, p-n peak [V]
Ib_p_ph = Sb/Vb/sqrt(3)*sqrt(2);       %L-n peak Converter base current [A]
Zb = Vb^2/Sb;                          %Converter base impedance [Ohm]
Lb = Zb/w_b;                           %base inductance [H]
Cb = 1/(Zb*w_b);                       %base capacitance [F]
%% setting points
P_set = P_set_pu*Sb; Q_set = Q_set_pu*Sb;   %Set points [W] and [V]
%% DC side parameters
Vdc = 3*Vb_p_ph;                               %DC-voltage reference [V]
%% Converter Parameters (already transformed into no-pu values)
% VSM
J = J_pu*(Sb/w_b);                     % inertia 
D_p = D_p_pu * (Sb/w_b);               % active power damping
D_q = D_q_pu * (Sb/Vb_p_ph);           % reactive power damping
%inner voltage control
kp_vc = kp_vc_pu/Vb_p_ph*Ib_p_ph;
ki_vc = ki_vc_pu/Vb_p_ph*Ib_p_ph;
%inner current control loop 
kp_i = kp_i_pu*Vb_p_ph/Ib_p_ph;                           %proportional gain of PI
ki_i = ki_i_pu*Vb_p_ph/Ib_p_ph;                            %integral gain of PI 
%
r_f = r_f*Zb;
c_f = c_f*Cb;                                           %output filter capacitor [pu]
l_f = l_f*Lb;  
%% P = 3/2*(Vod*Iod+Voq*Ioq), Q = 3/2*(-Vod*Ioq+Voq*Iod), Voq =0
Vod = Vb_p_ph; Voq =0;
%
b = [P_set;Q_set]; Func = 3/2*[Vod,Voq;Voq,-Vod];
x = inv(Func)*b;
Iod = x(1); Ioq = x(2);
%
Func = [1           0           -r_f          w_b*l_f;
     0           1           -w_b*l_f     -r_f;
     0           0           1             0;
     0           0           0             1];
b = [1 0 0 0;
    0 1 0 0;
    0 -w_b*c_f 1 0;
    w_b*c_f 0 0 1]*[Vod; Voq; Iod; Ioq];
x = Func \ b;
Vcd = x(1); Vcq = x(2); Icd = x(3); Icq = x(4);

I = [1 0;0 1];
%
syms s;
% synchronization
Gp = -1/s/(J*s+D_p);
Giop = 3/2*[Vod,Voq]; Guop = 3/2*[Iod,Ioq];
% axis transformation
Giovo = [Voq;-Vod]*Gp*Giop;
Gvovo = [Voq;-Vod]*Gp*Guop;
G1 = (I-Gvovo)^(-1);
G2 = (I-Gvovo)^(-1)*Giovo;
Gioio = [Ioq;-Iod]*Gp*Giop;
Gvoio = [Ioq;-Iod]*Gp*Guop;
G3 = (I-Gioio)^(-1);
G4 = (I-Gioio)^(-1)*Gvoio;
Gusuc = (I-G2*G4)^(-1)*G1;
Gisuc = (I-G2*G4)^(-1)*G2*G3;
Gisic = (I-G4*G2)^(-1)*G3;
Gusic = (I-G4*G2)^(-1)*G4*G1;
% more axis transformation
Gioic = [Icq;-Icd]*Gp*Giop;
Gvoic = [Icq;-Icd]*Gp*Guop;
Gisicc = Gioic*Gisic+Gvoic*Gisuc;
Gusicc = Gioic*Gusic+Gvoic*Gusuc;
Giovc = [Vcq;-Vcd]*Gp*Giop;
Gvcvc = [Vcq;-Vcd]*Gp*Guop;
Gisvcc = Giovc*Gisic+Gvcvc*Gisuc;
Gvsvcc = Giovc*Gusic+Gvcvc*Gusuc;
% reactive power control
Gq = -1/D_q;
Gioq = 3/2*[Voq,-Vod]; 
Guoq = 3/2*[-Ioq,Iod];
Gqioq = Gq*Gioq; Gquoq = Gq*Guoq;
G11 = Gqioq*Gusic+Gquoq*Gusuc;
G12 = Gqioq*Gisic+Gquoq*Gisuc;
% voltage control loop
PI_VC = kp_vc + ki_vc/s;
H1 = [PI_VC*G11;0,0]-[PI_VC,0;0,PI_VC]*Gusuc;
H2 = [PI_VC*G12;0,0]-[PI_VC,0;0,PI_VC]*Gisuc;
% current control loop
pi_CC = kp_i + ki_i/s;
PI_CC = [pi_CC,0;0,pi_CC];
H11 = PI_CC*H1-PI_CC*Gusicc+Gusuc;
H12 = PI_CC*H2-PI_CC*Gisicc+Gisuc;
% filter
GRLf = [r_f+s*l_f -w_b*l_f;w_b*l_f r_f+s*l_f];
Y_C = [s*c_f -w_b*c_f;w_b*c_f s*c_f];
% global axis
e_Jdelta = [cos(theta_0), -sin(theta_0); sin(theta_0), cos(theta_0)];
Y = (H12-Gisvcc-PI_CC-GRLf)^(-1)*(H11-Gvsvcc-PI_CC*Y_C-I-GRLf*Y_C);
YGFM = e_Jdelta*Y*(e_Jdelta)^(-1);
YGFL_num = matlabFunction(YGFM);
end