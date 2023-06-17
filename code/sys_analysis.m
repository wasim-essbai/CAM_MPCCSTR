clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%%  0. Parametri del reattore
q0 = 10;
V = 150;
k0 = 7.2*10^(10);
ER = 8750;
dHr = -50000;
UA = 50000;
ro = 1000;
Cp = 0.239;
tau_c = 1.5;
CAf = 1;
Tf = 370;

Ts = 1;

%%  1. Equazioni di stato
syms Ca T Tc Tr
syms dCa dT dTc
syms f1(Ca,T,Tc)
syms f2(Ca,T,Tc)
syms f3(Ca,T,Tc)

eq1 = (q0/V)*(CAf - Ca)-k0*exp(-ER/T)*Ca == 0;
eq2 = (q0/V)*(Tf - T) + (((-dHr)*k0)/(ro*Cp))*exp(-ER/T)*Ca + UA/(V*ro*Cp)*(Tc - T) == 0;
eq3 = (Tr - Tc)/tau_c == 0;

f1(Ca,T,Tc) = q0/V*(CAf - Ca)-k0*exp(-ER/T)*Ca;
f2(Ca,T,Tc) = q0/V*(Tf - T) + (((-dHr)*k0)/(ro*Cp))*exp(-ER/T)*Ca + UA/(V*ro*Cp)*(Tc - T);
f3(Ca,T,Tc) = (Tr - Tc)/tau_c;

%% 2. Sistema linearizzato
Ca_eq = 0.5054;
T_eq = 315.5491;
Tc_eq = 308;
Tr_eq = 308;

A = [-k0*exp(-ER/T_eq) - q0/V, -Ca_eq*ER*k0*exp(-ER/T_eq)/(T_eq^2), 0;
      ((-dHr)*k0/(ro*Cp))*exp(-ER/T_eq), -q0/V + Ca_eq*((-dHr)*k0/(ro*Cp))*exp(-ER/T_eq)*ER/(T_eq^2) - UA/(V*ro*Cp), UA/(V*ro*Cp);
      0, 0, -1/tau_c];
B = [0;0;1/tau_c];
C = [1, 0, 0];
D = 0;

continuous_time_ss = ss(A,B,C,D);

discrete_time_ss = c2d(continuous_time_ss, Ts);

% Matrice di raggiungibilità
Mr = ctrb(discrete_time_ss);
Mr_rank = rank(Mr);

% Autovalori del sistema
lambda0 = eig(discrete_time_ss.A);

%% 3. Vincoli sul sistema linearizzato
Hx = [1,0,0;
      -1,0,0;
      0,1,0;
      0,-1,0;
      0,0,1;
      0,0,-1];
hx = [0.954-Ca_eq;-0.38-Ca_eq; 10000; -10000; 10000; -10000];

Hu = [1;-1];
hu = [310 - Tr_eq; -280 - Tr_eq];

%% 4. Progettazione costo MPC

% Matrici del costo quadratico
Q = eye(3);
R = 1;

%   Computazione control invariant set
[CIS_H,CIS_h] = cis(A,B,[0; 0; 0],0,Hx,hx,Hu,hu,Q,R);
CIS = Polyhedron(CIS_H,CIS_h);
figure(1)
CIS.plot();
title('\textbf{Control invariant set sistema linearizzato}');
xlabel('$\theta$ [rad]');
ylabel('$\dot{\theta}$ [rad/s]');
xlim([-pi,pi]);
ylim([-2,2]);
