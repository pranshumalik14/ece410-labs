%% ece410: linear control systems
%  lab3: State Feedback Stabilization of a Cart-Pendulum Robot
%  authors: Pranshu Malik and Varun Sampat
%  date: 24 November 2021

clc;
close all;
clear all;

%% parameters setup
% set up initial values and DE parameters
parameters = struct('M', 1.0731, 'm', 0.2300 , 'l', 0.3302, 'g', 9.8);

%% Syms setup
syms u M m l g
syms x1 x2 x3 x4
% x1 is y
% x2 is y_dot
% x3 is theta
% x4 is theta_dot
x    = [x1; x2; x3; x4];

x1_dot = x(2);
x2_dot = (- (m*l*sin(x(3))*(x(4))^2) + (m*g*sin(x(3))*cos(x(3))) + u)/(M + m*sin(x(3))^2);
x3_dot = x(4);
x4_dot = (- (m*l*sin(x(3))*cos(x(3))*(x(4))^2) + (m+M)*g*sin(x(3)) + u*cos(x(3)) )/ ( l * (M + m*sin(x(3))^2));

x_dot = [x1_dot; x2_dot; x3_dot; x4_dot]

% y = x1;
%% Linearize at (0,0,0,0)
x_star = [0; 0; 0; 0];
u_star = 0;

symA = jacobian(x_dot,x);
symB = jacobian(x_dot,u);


% substitute the symbolic variables with the equilibrium points 
A = subs(symA, {x1, x2, x3, x4, u}, {x_star(1), x_star(2), x_star(3), x_star(4), u})
B = subs(symB, {x1, x2, x3, x4, u}, {x_star(1), x_star(2), x_star(3), x_star(4), u})

%% parameters
Asub = subs(A, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l})
Bsub = subs(B, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l})


%% Test for controllability
Qc = ctrb(Asub, Bsub)

% Rank of Qc tells us if the system is controllable
disp(rank(Qc))

%% Generate K1
p = [-1 -2 -3 -4]
K1 = place(Asub, Bsub, p)