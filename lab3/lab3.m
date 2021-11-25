%% ece410: linear control systems
%  lab3: state feedback stabilization of a cart-pendulum robot
%  authors: Pranshu Malik and Varun Sampat
%  date: 24 November 2021

clc;
close all;
clear all;

%% parameters setup

% set up initial values and DE parameters
parameters = struct('M', 1.0731, 'm', 0.2300 , 'l', 0.3302, 'g', 9.8);

%% symbolic linearization of the control system

syms u M m l g
syms x1 x2 x3 x4

% x1 = y; horizontal position
% x2 = y_dot
% x3 = theta; angle from the vertically upwards axis
% x4 = theta_dot
x    = [x1; x2; x3; x4];

x1_dot = x(2);
x2_dot = (-(m*l*sin(x(3))*(x(4))^2) + (m*g*sin(x(3))*cos(x(3))) + u) / (M + m*sin(x(3))^2);
x3_dot = x(4);
x4_dot = (-(m*l*sin(x(3))*cos(x(3))*(x(4))^2) + (m+M)*g*sin(x(3)) + u*cos(x(3))) / (l * (M + m*sin(x(3))^2));

x_dot = [x1_dot; x2_dot; x3_dot; x4_dot];

% linearization equilibrium (for all comparisons going forward)
xstar = [0; 0; 0; 0];

symA = jacobian(x_dot,x);
symB = jacobian(x_dot,u);


% substitute the symbolic variables with the equilibrium points
A = subs(symA, {x1, x2, x3, x4, u}, {xstar(1), xstar(2), xstar(3), xstar(4), u});
B = subs(symB, {x1, x2, x3, x4, u}, {xstar(1), xstar(2), xstar(3), xstar(4), u});

%% converting to numeric state-space representation

numA = double(subs(A, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l}));
numB = double(subs(B, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l}));

%% test for controllability

Qc = ctrb(Asub, Bsub)

% Rank of Qc tells us if the system is controllable
disp(rank(Qc))

%% generate K1

p = [-1 -2 -3 -4]
K1 = place(Asub, Bsub, p)

%% Simulate at x(-0.5, 0, -pi/4, 0)
x0 = [-0.5;0;-pi/4; 0];

% setup Tspan
options    = odeset('RelTol',1e-7, 'AbsTol',1e-7);
Tspan      = linspace(0,20,2e3);

[t, xf] = ode45(@inverted_pendulum, Tspan, x0, options, K1, parameters);

x1 = xf(:,1);
x2 = xf(:,2);
x3 = xf(:,3);
x4 = xf(:,4);

plot(t, x1)

%% u over time
u_num = (K1 * xf')';
plot(t, u_num)
