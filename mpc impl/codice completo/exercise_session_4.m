clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%%  0. Parametri del pendolo
g = 9.81;
l = 0.5;
m = 0.1;
b = 1e-2;

%%  1. Invariant set per il sistema controllato

%   Sampling time
Ts = 0.1;

%   Sistema linearizzato e discretizzato (Eulero in avanti)
A = [1 Ts; Ts*g/l -Ts*b/(m*l^2)];
B = [0; Ts/l];

%   Vincoli su stato e ingresso (NELLE COORDINATE DEL SISTEMA
%   LINEARIZZATO!)
Hx = [eye(2); -eye(2)];
hx = [pi/2;2;pi/2;2];
Hu = [1; -1];
hu = 5*ones(2,1);

%   Matrici del costo quadratico
Q = eye(2);
R = 1;

%   Computazione control invariant set
[CIS_H,CIS_h] = cis(A,B,[0; 0],0,Hx,hx,Hu,hu,Q,R);
CIS = Polyhedron(CIS_H,CIS_h);
figure(1)
CIS.plot();
title('\textbf{Control invariant set sistema linearizzato}');
xlabel('$\theta$ [rad]');
ylabel('$\dot{\theta}$ [rad/s]');
xlim([-pi,pi]);
ylim([-2,2]);

%%  2. N-step-controllable set dell'invariant set

%   Orizzonte di predizione
N = 5;

[Np_steps_H, Np_steps_h] = controllable_set(Hx,hx,Hu,hu,CIS_H,CIS_h,A,B,N);
Np_step_set = Polyhedron('A',Np_steps_H,'b',Np_steps_h);
figure(2)
Np_step_set.plot('Alpha',0);
title('\textbf{CIS e 10-step-controllable set}');
xlabel('$\theta$ [rad]');
ylabel('$\dot{\theta}$ [rad/s]');
hold on;
CIS.plot()
xlim([-pi,pi]);
ylim([-2,2]);

legend({'Control-invariant set', '10-step set'});

%%  3. Design MPC e simulazione

%   Numero di step simulati
T_sim = 60;

%   Riferimento
x_ref = [pi; 0];
u_ref = 0;

% N.B. RIFERIMENTO DEL SISTEMA LINEARIZZATO!!
x_ref_lin = [0; 0];
u_ref_lin = 0;

mpc = mpc_ingredients(A,B,Hx,hx,Hu,hu,x_ref_lin,u_ref_lin,Q,R,N);
mpc2 = mpc_ingredients(A,B,Hx,hx,Hu,hu,x_ref_lin,u_ref_lin,Q,R,10);

%   Log stati e ingresso sistema
x_log = zeros(2,T_sim+1);
u_log = zeros(1,T_sim);
flags = zeros(1,T_sim);

x_log(:,1) = [pi-0.4; -0.2];

for tt = 1:T_sim

    %   Stato sistema linearizzato
    x_lin = x_log(:,tt) - [pi; 0];
    x_lin_shifted = x_lin - x_ref_lin;

    %   Impostazioni MPC relative alla condizione iniziale
    f = mpc.f_base * x_lin_shifted;
    b_ineq = mpc.b_ineq_base - mpc.b_ineq_x0_factor * x_lin_shifted;
    
    % S = Polyhedron('A',mpc.A_ineq,'b',b_ineq);
    % isEmptySet(S)
    %   Risoluzione problema di ottimizzazione
    [delta_u_seq,~,exitflag] = quadprog(mpc.F,f,mpc.A_ineq,b_ineq);

    flags(tt) = exitflag;

    %   Risposta del sistema
    u_log(tt) = u_ref + delta_u_seq(1);
    
    dxdt = @(t,x) pendulum(t,x,u_log(tt),g,l,m,b);
    [~,xx] = ode45(dxdt,[0 Ts],x_log(:,tt));
    x_log(:,tt+1) = xx(end,:)';

end

%%  4. Plot dei risultati

%   Traslazione del CIS e dell'N-step set nelle coordinate originali
CIS_shifted = CIS + x_ref;
Np_step_set_shifted = Np_step_set + x_ref;

figure(3);
Np_step_set_shifted.plot('Alpha',0);
title('\textbf{Traiettoria del sistema}');
xlabel('$\theta$ [rad]');
ylabel('$\dot{\theta}$ [rad/s]');
hold on;
CIS_shifted.plot()
hold on
plot(x_log(1,:),x_log(2,:),'Color',[0 0 0.5])
scatter(x_log(1,:),x_log(2,:),'cyan');
xlim([0,2*pi]);
ylim([-2,2]);