clear all
clc

%Parameter%
KRt=1; %ohm/K
TRt=1/250; %s
KR=1/100; %V/ohm
Ky=3/40; %bar/V
Ty=1/100; %s
Kz=2/100; %m/bar
Kh=1000; %(m^3/min)/m
Kpv=1; %(m^3/min)/bar
TD=1/2; %s
KD=20; %K/(m^3/min)
KW=1; %K/(1/min)
TW1=20; %s
TW2=50; %s
Tpv = 100; %s
Tmw = 300; %s
Tan = 30; %s
Taus = 50; %s

%%%%%%%%%%%%%%%%%% SYSTEMANALYSE %%%%%%%%%%%%%%%%%%%

%Zustandsraummodell%

%x=Py,m.D,ta,t.a,Rt
A=[-1/Ty 0 0 0 0;
   -(Kh*Kz)/TD -1/TD 0 0 0;
    0 0 0 1 0;
    0 KD/(TW1*TW2) -1/(TW1*TW2) -(TW1+TW2)/(TW1*TW2) 0;
    0 0 KRt/TRt 0 -1/TRt];

Api = [-1/Ty 0 0 0 0;
   -(Kh*Kz)/TD -1/TD 0 0 0;
    0 0 0 1 0;
    0 KD/(TW1*TW2) -1/(TW1*TW2) -(TW1+TW2)/(TW1*TW2) 0;
    0 0 KRt/TRt 0 -1/TRt]; %it is the same as A, because in our case,
%workaround method was not possible, since the system included no
%I0 strecke. 

%u=uy
B=[Ky/Ty;
   0;
   0;
   0;
   0];

Bpi = [Ky/Ty;
   0;
   0;
   0;
   0]; %it is the same as B, because in our case,
%workaround method was not possible, since the system included no
%I0 strecke. 


%x=Py,m.D,ta,t.a,Rt
C = [0 0 0 0 KR];

C_R = [0 0 1 0 0];

Ch_R = C_R;

%u=uy
D = [0];

%z=Pv,m.w
E=[0 0;
   Kpv/TD 0;
   0 0;
   0 -KW/(TW1*TW2);
   0 0];

%z=Pv,m.w
F=[0 0];

sysORORK=ss(A,B,C,D);

%pzmap%
figure(1)
pzmap(sysORORK)

%step%
figure(2)
step(sysORORK)

%%%%%%%%%%%%% REGLERSTRUKTURFESTLEGUNG %%%%%%%%%%%%%%%%

%%%%%% CONTROLLABILTY %%%%%%%%

Qs = ctrb(A,B); %ctrb computes a controllability matrix from state
% matrices or from a state-space model. One can use this matrix 
% to determine controllability.

r = rank(Qs);  %The rank function provides an estimate of the 
% number of linearly independent rows or columns of a matrix.


% Check if the system is fully controllable
if rank(Qs) == length(A)
    disp('The system is fully controllable.');
else
    disp('The system is not fully controllable.');
end

%%%%%%%%% OBSERVABILITY %%%%%%%%%%%

Qb = obsv(A,C);  %obsv computes an observability matrix from 
%the state matrices.

rb = rank(Qb);

l = length(A);

% Check if the system is observable
if rank(Qb) == length(A)
    disp('The system is observable.');
else
    disp('The system is not observable.');
end

%%%%%%%%%%%%% Zustandsraumregelung %%%%%%%%%%%%%

%%%%%% Determination of the state space feedback R %%%%%

eig_A = eig(A); %eigenvalues are calculated (they are the poles)

P_set = eig_A-1; %since some pole zeroes were too close to the
% y axis, P_set = eig_A-1 was deemed necessary.

R_place = place(A,B,P_set);  % Calculates the controller matrix R.

%%%%%% Determination of Vorfilter Matrix M %%%%%%

M = (C*(B*R_place-A)^(-1)*B)^(-1);

%%%%%%%%%%% Zustandraumreglung mit Beobachter %%%%%%%%%%%%

%%% Determination of the matrixes %%%

Ah = A;
Bh = B;
Ch = C;

%%% Observability Matrix L %%%

P_set_obs = eig_A;

L_luen_T = place(Ah',Ch',P_set_obs); % L (observability matrix) 
% is calculated but transposed.

L_luen = L_luen_T'; % L = L transposed transposed.

%%%%%%%% Zustandraumreglung mit Beobachter und Störgrößen %%%%%%%%%%%

%%% Determination of Störgrößenmatrix Z %%%%%%

Z = ((B')*B)^(-1)*(B')*E;

%%%%% PI ZUSTANDSGRÖßENREGLER %%%%%%%

A_PI_real = [A zeros(5,1);
      -C_R 0];

B_PI_real = [B;
            0];

P_set_PI = eig(A_PI_real)-10;
R_place_PI = place(A_PI_real,B_PI_real,P_set_PI);

R_I = -R_place_PI(:,6);
R_p = -(C_R*Api^-1*Bpi)^-1;
R_x = R_place_PI(:,1:5)-R_p*C_R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% REGLERPARAMETERBERECHNUNG %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Controller parameter calculation
R = diag(1); % R is the control cost matrix (an Einheitsmatrix), Matrix
% R should be square and have the same number of rows as matirx B

%%%%%%%%%% Minimum Damping Maximum Time Constant %%%%%%%%%

%%%%% Determination of D %%%%%%%
K = 50; 
T = 2; 
D_var=[0.7:0.01:0.75]; 
for n=1:1:length(D_var) 
PT2_varD = tf([K],[T^2 2*D_var(n)*T 1]); 
[yD_varD{n},tOutD_varD{n}] = step(PT2_varD, 500);
end

figure(3)
hold on
for n=1:1:length(yD_varD) 
plot(tOutD_varD{n},yD_varD{n})
end
legend('D=0,7','D=0,71','D=0,72','D=0,73','D=0,74','D=0,75','Location','northwest')

% Toleranzband
yline(52, '--r', 'Oberes Toleranzband','HandleVisibility', 'off');
yline(48, '--r', 'Unteres Toleranzband','HandleVisibility', 'off');
yline(50, '--p','HandleVisibility', 'off');
xlim([5, 12]);
ylim([47, 53]);

title("Dmin");

hold off

Dmin = 0.72;

%%%%% Determination of T %%%%%%%

T_var = [8:0.5:10];

figure(4)
hold on
for n=1:1:length(T_var)
 T_max = T_var(n);
 PT2_varT = tf([K],[T_max^2 2*Dmin*T_max 1]);
 step(PT2_varT);
end
legend('T=8s','T=8,5s','T=9s','T=9,5s','T=10s','Location','northwest')

% Tan und Taus

xline(30, '--r', 'Tan','HandleVisibility', 'off');
xline(50, '--r', 'Taus','HandleVisibility', 'off');
yline(48, '--r', 'Unteres Toleranzband','HandleVisibility', 'off');
ylim([45, 55]);
xlim([25, 55]);
title("Tmax");
hold off

Tmax = 9.5;

%%%%%%% GRK %%%%%%%%

figure(5)
hold on
GRK1 = tf([K],[Tmax^2 2*Dmin*Tmax 1]);
step(GRK1)
grid on
xline(30, '--r', 'Tan');
xline(50, '--r', 'Taus');
yline(52, '--r', 'Oberes Toleranzband');
yline(48, '--r', 'Unteres Toleranzband');
ylim([45, 55]);
xlim([25, 60]);
title("GRK mit Tmax und Dmin");

%%%%%%%%%%%%%%%%%%%% Polplatzierung %%%%%%%%%%%%%%%%%%%

%%%%%%%%% D = 0.72, T = 9.5 %%%%%%%%%%%

D_set1 = Dmin;
P11 = -1/Tmax*(cos(acos(D_set1)) + j*sin(acos(D_set1)));
P21 = -1/Tmax*(cos(acos(D_set1)) - j*sin(acos(D_set1)));

P_set_min1 = [-251; P11 ; P21 ;-3; -101];
R_place_opt1 = place(A,B,P_set_min1);
M1 = (C_R*(B*R_place_opt1-A)^(-1)*B)^(-1);

%%%%%%%%% D = 0.72, T = 5 %%%%%%%%%%%

D_set2 = Dmin;
P12 = -1/5*(cos(acos(D_set2)) + j*sin(acos(D_set2)));
P22 = -1/5*(cos(acos(D_set2)) - j*sin(acos(D_set2)));

P_set_min2 = [-251; P12 ; P22 ;-3; -101];
R_place_opt2 = place(A,B,P_set_min2);
M2 = (C_R*(B*R_place_opt2-A)^(-1)*B)^(-1);

%%%%%%%%% D = 0.9, T = 9.5 %%%%%%%%%%%

D_set3 = 0.9;
P13 = -1/Tmax*(cos(acos(D_set3)) + j*sin(acos(D_set3)));
P23 = -1/Tmax*(cos(acos(D_set3)) - j*sin(acos(D_set3)));

P_set_min3 = [-251; P13 ; P23 ;-3; -101];
R_place_opt3 = place(A,B,P_set_min3);
M3 = (C_R*(B*R_place_opt3-A)^(-1)*B)^(-1);

%%%%%%%%% D = 0.9, T = 5 %%%%%%%%%%%

D_set4 = 0.9;
P14 = -1/5*(cos(acos(D_set4)) + j*sin(acos(D_set4)));
P24 = -1/5*(cos(acos(D_set4)) - j*sin(acos(D_set4)));

P_set_min4 = [-251; P14 ; P24 ;-3; -101];
R_place_opt4 = place(A,B,P_set_min4);
M4 = (C_R*(B*R_place_opt4-A)^(-1)*B)^(-1);

%%%%%% Bestmögliche Kombination von D und T ist D=0.9 und T=5!  %%%%%%%

%%%%%%%%%%%%%%% Matrix Ricatti Regler (LQR) %%%%%%%%%%%%%

Q1 = diag([10^0 10^0 10^0 10^0 10^0]);
Q2 = diag([10^0 10^0 10^1 10^4 10^5]);
Q3 = diag([10^0 10^0 10^1 10^2 10^2]);
Q4 = diag([10^0 10^1 10^-1 10^0 10^0]);

R_place_lqr1 = lqr(A,B,Q1,R);
M_lqr1 = (C_R*(B*R_place_lqr1-A)^(-1)*B)^(-1);

R_place_lqr2 = lqr(A,B,Q2,R);
M_lqr2 = (C_R*(B*R_place_lqr2-A)^(-1)*B)^(-1);

R_place_lqr3 = lqr(A,B,Q3,R);
M_lqr3 = (C_R*(B*R_place_lqr3-A)^(-1)*B)^(-1);

R_place_lqr4 = lqr(A,B,Q4,R);
M_lqr4 = (C_R*(B*R_place_lqr4-A)^(-1)*B)^(-1);

%%%%%%%%%%%%%%%%%%%Q1 is the best%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Luenburger %%%%%%%%%%%%%%%%%


% P1 = -1/Tmax*(cos(acos(Dmin)) + j*sin(acos(Dmin)))-0.15;
% P2 = -1/Tmax*(cos(acos(Dmin)) - j*sin(acos(Dmin)))-0.15;

P1 = -1/5*(cos(acos(0.9)) + j*sin(acos(0.9)));
P2 = -1/5*(cos(acos(0.9)) - j*sin(acos(0.9)));

% Luenberger-Beobachter:

% Setting up the specific pole points for Luenberger-Beobachter

P_set_L1 = [-250; P1; P2; -2; -100];

% Controller matrix Calculation
RL = place(A,B,P_set_L1);

% Prefilter Matrix M Calculation
ML = (C*(B*RL-A)^(-1)*B)^(-1);

% Calculation of different observebility matrix parameters

PsollL1 = [P_set_L1 - 0.0001];
PsollL2 = [P_set_L1 - 0.0005];
PsollL3 = [P_set_L1 - 0.001];

% Transponierte Werte
L1T = place(Ah',Ch',PsollL1);
L2T = place(Ah',Ch',PsollL2);
L3T = place(Ah',Ch',PsollL3);

% L = L Transponiert Transponiert
L1=L1T';
L2=L2T';
L3=L3T';

%%%%%%% PsollL2 gives the best possible result, the X - X^ value
% reaches to 0 around 30s. %%%%%%%%%

%%%%%%%%%%%%%% KALMANS BUCKY %%%%%%%%%%%%%%%%

Qk1 = 10^0*diag([10^0 10^0 10^0 10^0 10^0]);
Qk2 = 10^1*diag([10^0 10^0 10^0 10^0 10^0]);
Qk3 = 10^2*diag([10^0 10^0 10^0 10^0 10^0]);
Qk4 = 10^3*diag([10^0 10^0 10^0 10^0 10^0]);
Qk5 = 10^4*diag([10^0 10^0 10^0 10^0 10^0]);

Rk = 1;

%%%%%%%%% Kalman Filter Calculation %%%%%%%%%%
[L_Kalman S0 Eig_beo] = lqr(Ah', Ch', Qk5, Rk);

L_K = L_Kalman';
R_K = RL;
M_K = ML;

%%%%%%% Mit den gewählten Werten stabilisiert sich das System bei etwa 
% t=35 Sekunden. Die Grafik entspricht den Gütekriterien.

%%%%%%% Die X - X^ Grafik stabiliziert sich bei etwa t=10 Sekunden. Das 
% entspricht für einen schnellen Beobachter, der gut funtioniert.

%%%%%%%%%%%%%%%% REGLERVERHALTENSBEWERTUNG %%%%%%%%%%%%%%%%

% Definition der Parameter:

% Zeitkonstanten:
TRt = 1/250;    % einheit = s
Ty = 1/100;     % einheit = s
TD = 1/2;       % einheit = s
TW1 = 20;       % einheit = s
TW2 = 50;       % einheit = s

% Parameter:
KRt = 1;          % einheit = Ω/K
KR = 1/100;       % einheit = V/Ω
Ky = 3/40;        % einheit = bar/V
Kz = 2/100;       % einheit = m/bar
Kh = 1000;        % einheit = 1000 (m^3/min)/m
Kpv = 1;          % einheit = (m^3/min)/bar
KD = 20;          % einheit = K/(m^3/min)
KW = 1;           % einheit = K/(1/min)

% Definition des Arbeitspunkts:
Apu = 70+273; % einheit = K
Apo = 95+273; % einheit = K

% Definition der Arbeitspunktwerte:
TApu = 50+273;     % einheit = K
uTApu = 5;         % einheit = V
RTApu = 300;       % einheit = Ω
uYApu = 4;         % einheit = V
pYApu = 2.4;       % einheit = bar
hApu = 10;         % einheit = mm
mWApu = 30;        % einheit = l/min
pVApu = 1.6;       % einheit = bar
mDApu = 1.5;       % einheit = m^3/min

% Störungen:
ampl_pV = 1;           % einheit = bar
pdauer_pV = 100;       % einheit = s
ampl_m = 10;           % einheit = l/min
pdauer_m = 300;        % einheit = s

% Steuerkonstanten der Störgrößen:
z1str = 1;
z2str = 1;

% Steuerkonstante des Reglers:
Rstr = 1;

% Festlegung der Simulationsparameter:
t_start = 0;     % einheit = s
t_stop = 500;    % einheit = s
t_step = 1/500;  % einheit = s
t_jump = 10;     % einheit = s

% Festlegung der Qualitätskriterien:
dystep = Apo-Apu;
tol = 2;
xuemax = 5;
Tan = 30;
Taus = 50;

% Simulation in Simulink:

% Definiton der Systemmatrix:
A=[ -1/Ty           0               0                 0                 0;
   -(Kh*Kz)/TD    -1/TD             0                 0                 0;
      0             0               0                 1                 0;
      0        KD/(TW1*TW2)   -1/(TW1*TW2)  -(TW1+TW2)/(TW1*TW2)        0;
      0             0            KRt/TRt              0           -1/TRt];

B=[ Ky/Ty;
      0  ;
      0  ;
      0  ;
      0 ];

E=[  0            0      ;
   Kpv/TD         0      ;
     0            0      ;
     0      -KW/(TW1*TW2);
     0            0     ];

% Eingangsmatrix (Beide Führungs- und Störgrößen):
BE = [ Ky/Ty      0            0      ;
        0       Kpv/TD         0      ;
        0         0            0      ;
        0         0      -KW/(TW1*TW2); 
        0         0            0     ];

% Definition der Ausgangsmatrizen:
C = [0 0 0 0 KR]; %ORK

C_R = [0 0 1 0 0]; %GRK

% Definition der Durchgangsmatrix:
D = [0];

E = [0        0;
   Kpv/TD     0;
     0        0;
     0  -KW/(TW1*TW2);
     0        0];

F = [0        0];

%Eingenwerte ORORK:
eig_A =eig(A);

% Definition des Zustandsraummodells für den PI-Regler:
A_PI = [  A      zeros(5,1);
        -C_R         0    ];

B_PI = [ B ;
         0];

% Matrix Ricatti Regler:
QLQR = diag([10^0 10^0 10^0 10^0 10^0 10^13]);

RLQR = 1;

% LQR Parameter Berechnung:
[R_PI, S, E_LQR] = lqr(A_PI, B_PI, QLQR, RLQR);

% Reglerparameter für den PI-Regler:
Ri = -R_PI(1,6);
Rp = -(C_R*A^(-1)*B(:,1))^(-1);
Rx = R_PI(1,1:5)-Rp*C_R;

% Beobachtermodell (zero vectors were added to match the dimensions 
% of the matrix with the entrance values):
A_h = A;
B_h = [B,[0;0;0;0;0],[0;0;0;0;0]];
C_h = C;
D_h = D;

xs = [0;0;0;0;0];
xb = [10;10;10;10;10];

%%%%%%%%%%%%%FEINTUNING%%%%%%%%%%%

% Matrix Ricatti Regler:
QLQR_F = diag([10^0 10^0 10^0 10^0 10^0 10^0]);
RLQR_F = 1;

% LQR Parameter Berechnung:
[R_PI_F, S, E_LQR] = lqr(A_PI, B_PI, QLQR_F, RLQR_F);

% Reglerparameter für den PI-Regler:
Ri_F = -R_PI_F(1,6);
Rp_F = -(C_R*A^(-1)*B(:,1))^(-1);
Rx_F = R_PI_F(1,1:5)-Rp_F*C_R;

% Berechnung der Bewertungsmatrix für Beobachter:
Q_Kalman = 10^4*diag([10^0 10^0 10^5 10^0 10^5]);
R_Kalman = 1;

% Berechnung des Kalman-Filters:
[LKalmanT SKalman EKalman]= lqr(A_h',C_h',Q_Kalman,R_Kalman);

% Berechnung von step time:
t_step = 1/abs(min(real(EKalman)));

% Zuweisung zu den Parametern der Simulation:
LKalman = LKalmanT';
Rx = Rx;

%%%%%%%%%%%%%%%%DIAGRAMME AUS SIMULINK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simout = sim('LNB_ET_Gr3_Simulink.slx');

% Darstellung der Sprungantwort aus Simulink
figure(6);
plot(simout.ScopeData1{1}.Values.Time,simout.ScopeData1{1}.Values.Data);
title('GRK mit Zustandsraumregler');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
ylim([0, 50]);
xlim([0, 500]);
grid on

figure(7);
plot(simout.ScopeData2{1}.Values.Time,simout.ScopeData2{1}.Values.Data);
title('Differenz von Zustandsgrößen GRK mit Beobachter');
xlabel('Zeit (s)');
ylim([-10, 50]);
xlim([0, 500]);
legend('Py (bar)',"m'D (m^3/min)","ta (°C)","t'a (°C/s)",'Rt (ohm)')
grid on

figure(8);
plot(simout.ScopeData3{1}.Values.Time,simout.ScopeData3{1}.Values.Data);
title('Ausgangsgrößen GRK mit Beobachter');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
legend('GRK','Beobachter','Location','southeast')
ylim([0, 50]);
xlim([0, 500]);
grid on

figure(9);
plot(simout.ScopeData4{1}.Values.Time,simout.ScopeData4{1}.Values.Data);
title('Differenz der Zustandsgrößen GRK mit Störgrößen');
xlabel('Zeit (s)');
ylim([-10, 50]);
xlim([0, 500]);
legend('Py (bar)',"m'D (m^3/min)","ta (°C)","t'a (°C/s)",'Rt (ohm)')
grid on

figure(10);
plot(simout.ScopeData5{1}.Values.Time,simout.ScopeData5{1}.Values.Data);
title('Ausgangsgrößen GRK mit PI-Regler');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
legend('GRK','Beobachter','Location','southeast')
ylim([0, 50]);
xlim([0, 500]);
grid on

figure(11);
plot(simout.ScopeData6{1}.Values.Time,simout.ScopeData6{1}.Values.Data);
title('Polplatzierung mit verschiedenen D und T');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
ylim([34, 48]);
xlim([0, 500]);
xline(t_jump+30, '--r','HandleVisibility', 'off');
xline(t_jump+50, '--r','HandleVisibility', 'off');
yline(43, '--r','HandleVisibility', 'off');
yline(45, '--g','HandleVisibility', 'off');
yline(47, '--r','HandleVisibility', 'off');
legend('D=0,72 T=9,5s','D=0,72 T=5s','D=0,9 T=9,5s','D=0,9 T=5s',"Location","southeast")
grid on

figure(12);
plot(simout.ScopeData7{1}.Values.Time,simout.ScopeData7{1}.Values.Data);
title('LQR mit verschiedenen Q-Matrizen');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
ylim([34, 48]);
xlim([0, 500]);
xline(t_jump+30, '--r','HandleVisibility', 'off');
xline(t_jump+50, '--r','HandleVisibility', 'off');
yline(43, '--r','HandleVisibility', 'off');
yline(45, '--g','HandleVisibility', 'off');
yline(47, '--r','HandleVisibility', 'off');
legend('Q1','Q2','Q3','Q4',"Location","southeast")
grid on

figure(13);
plot(simout.ScopeData8{1}.Values.Time,simout.ScopeData8{1}.Values.Data);
title('Differenz zwischen der Zustandgrößen mit Luenberger-Beobachter');
xlabel('Zeit (s)');
ylim([-3, 10]);
xlim([0, 50]);
legend('Py (bar)',"m'D (m^3/min)","ta (°C)","t'a (°C/s)",'Rt (ohm)')
grid on

figure(14);
plot(simout.ScopeData9{1}.Values.Time,simout.ScopeData9{1}.Values.Data);
title('Wirkung von dem Luenberger-Beobachter');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
ylim([-0.025, 0.005]);
xlim([0, 50]);
grid on

figure(15);
plot(simout.ScopeData10{1}.Values.Time,simout.ScopeData10{1}.Values.Data);
title('Wirkung von dem Kalman-Bucy');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
ylim([-0.6, 1]);
xlim([0, 500]);
grid on

figure(16);
plot(simout.ScopeData11{1}.Values.Time,simout.ScopeData11{1}.Values.Data);
title('Regelgröße mit Kalman-Bucy-Filter');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
ylim([0, 50]);
xlim([0, 500]);
xline(t_jump+30, '--r','HandleVisibility', 'off');
xline(t_jump+50, '--r','HandleVisibility', 'off');
yline(43, '--r','HandleVisibility', 'off');
yline(45, '--g','HandleVisibility', 'off');
yline(47, '--r','HandleVisibility', 'off');
grid on

figure(17);
plot(simout.ScopeData12{1}.Values.Time,simout.ScopeData12{1}.Values.Data);
title('Verhalten von Zustandsgrößen Feintuning');
xlabel('Zeit (s)');
ylim([-10, 50]);
xlim([0, 500]);
xline(t_jump+30, '--r','HandleVisibility', 'off');
xline(t_jump+50, '--r','HandleVisibility', 'off');
yline(43, '--r','HandleVisibility', 'off');
yline(45, '--g','HandleVisibility', 'off');
yline(47, '--r','HandleVisibility', 'off');
legend('Py (bar)',"m'D (m^3/min)","ta (°C)","t'a (°C/s)",'Rt (ohm)','Location','east')
grid on

figure(18);
plot(simout.ScopeData13{1}.Values.Time,simout.ScopeData13{1}.Values.Data);
title('Verhalten von Regelgröße Feintuning');
xlabel('Zeit (s)');
ylabel('Temperatur (°C)');
ylim([40, 100]);
xlim([0, 500]);
xline(t_jump+30, '--r','HandleVisibility', 'off');
xline(t_jump+50, '--r','HandleVisibility', 'off');
yline(93, '--r','HandleVisibility', 'off');
yline(95, '--g','HandleVisibility', 'off');
yline(97, '--r','HandleVisibility', 'off');
grid on
