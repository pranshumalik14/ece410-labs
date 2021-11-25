%% ece410: linear control systems
%  lab3: state feedback stabilization of a cart-pendulum robot
%  authors: Pranshu Malik and Varun Sampat
%  date: 24 November 2021

clc;
close all;
clear all;

%% system parameters

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

% output 1: symbolic forms of A and B
simplify(symA)
simplify(symB)

% substitute the symbolic variables with the equilibrium points
A = subs(symA, {x1, x2, x3, x4, u}, {xstar(1), xstar(2), xstar(3), xstar(4), u});
B = subs(symB, {x1, x2, x3, x4, u}, {xstar(1), xstar(2), xstar(3), xstar(4), u});

%% converting to numeric state-space representation

numA = double(subs(A, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l}));
numB = double(subs(B, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l}));

%% test for controllability

Qc = ctrb(numA, numB);

% system is controllable if Qc is full rank
assert(rank(Qc) == min(size(Qc)));

%% generate K1

p = [-1 -2 -3 -4];
K1 = -1*place(numA, numB, p) % place adopts the convention that the state feedback controller has the form of u = -Kx, but we want it to be u = Kx

%% simulate for x0 = (-0.5, 0, -pi/4, 0)

x0 = [-0.5; 0; -pi/4; 0];

% setup Tspan
options    = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
Tspan      = linspace(0,20,2e3);

[t, X1] = ode45(@inverted_pendulum, Tspan, x0, options, parameters, {K1});

x1_t = X1(:,1);
x2_t = X1(:,2);
x3_t = X1(:,3);
x4_t = X1(:,4);
u_t  = X1(:,5);

figure
plot(t, x1_t)

figure
plot(t, u_t)
