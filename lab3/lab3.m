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

% generate K1
p1 = [-1 -2 -3 -4];
K1 = -1*place(numA, numB, p1) % use u = Kx instead of u = -Kx

% create cls (closed loop system) for K1
A_cls1 = numA + numB*K1;
cls1   = @(t, x) A_cls1*x;

% simulate for x0 = (-0.5, 0, -pi/4, 0)
x0 = [-0.5; 0; -pi/4; 0];
u0 = K1*x0;

% setup Tspan
options    = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
Tspan      = linspace(0,10,1e3);

% [t1, X1] =  ode45(@inverted_pendulum, Tspan, [x0; u0], options, parameters, K1);
[t1, X1] = ode45(cls1, Tspan, x0);

% extract state evolutions for K1
x1_t1 = X1(:,1);
x2_t1 = X1(:,2);
x3_t1 = X1(:,3);
x4_t1 = X1(:,4);
% u_t1  = X1(:,5);

% extract control input, u, for K1
u_t1 = X1*K1'; % = (K1*X1')'

% generate K2
p2 = [-1 -2 -3 -20];
K2 = -1*place(numA, numB, p2) % use u = Kx instead of u = -Kx

% create cls system for K2
A_cls2 = numA + numB*K2;
cls2   = @(t, x) A_cls2*x;

%[t2, X2] = ode45(@inverted_pendulum, Tspan, [x0; u0], options, parameters, K2);
[t2, X2] = ode45(cls2, Tspan, x0);

% extract state evolutions for K2
x1_t2 = X2(:,1);
x2_t2 = X2(:,2);
x3_t2 = X2(:,3);
x4_t2 = X2(:,4);
% u_t2 = X2(:,5);

% extract control input, u, for K2
u_t2 = X2*K2'; % = (K2*X2')'

figure
plot(t1, x1_t1)
hold on
plot(t2, x1_t2)

figure
plot(t1, u_t1)
hold on
plot(t2, u_t2)

%% Linear Quadratic Optimal control
syms q1 q2

symQ = [q1 0 0 0; 0 0 0 0; 0 0 q2 0; 0 0 0 0];

%% Impact of changing q1
R = 0.5;
num_q2 = 5;

% exp 1
num_q1_1 = 0.1;

Q_num_1 = double(subs(symQ, {q1, q2}, {num_q1_1, num_q2}));
K1_1 = -1*lqr(numA, numB, Q_num_1, R)

% exp 2
num_q1_2 = 0.005;
Q_num_2 = double(subs(symQ, {q1, q2}, {num_q1_2, num_q2}));
K1_2 = -1*lqr(numA, numB, Q_num_2, R);

%% Impact of changing q2
num_q1 = 0.05;
R = 0.5;

% exp 1
num_q2_1 = 1;
Q_num_1 = double(subs(symQ, {q1, q2}, {num_q1_1, num_q2_1}));
K2_1 = -1*lqr(numA, numB, Q_num_1, R);

% exp 2
num_q2_2 = 2000;
Q_num_2 = double(subs(symQ, {q1, q2}, {num_q1_1, num_q2_2}));
K2_2 = -1*lqr(numA, numB, Q_num_2, R);

%% Impact of changing R
num_q1 = 0.05;
num_q2 = 5;
Q_num = double(subs(symQ, {q1, q2}, {num_q1_1, num_q2_2}));

% exp 1
num_R_1 = 0.0005;
K3_1 = -1*lqr(numA, numB, Q_num, num_R_1);

% exp 2
num_R_2 = 10;
K3_2 = -1*lqr(numA, numB, Q_num, num_R_1);

% need a function that takes in K, Tspan, A, B, x0, and outputs generated
% plots

%% nonlinear comparison of controller performance

x0 = [-1; 0; -pi/4; 0]; % new initial state

% get nonlinear evolution of X given our control law, u = Kx
[t, Xref] = ode45(@inverted_pendulum, Tspan, [x0; u0], options, parameters, K3_1);

x_non_lin = Xref(:, 1:4);
u_non_lin = Xref(:, 5);

% create cls system for K3_1
A_cls = numA + numB*K3_1;
cls   = @(t, x) A_cls*x;

[t2, Xlin] = ode45(cls, Tspan, x0);
u_lin = Xlin*(K3_1');

figure
plot(t, Xlin(:,1) , t, x_non_lin(:,1))

%% 
x0_bar = [-200; 0; -pi/4; 0]; 
Tspan      = linspace(0,100,1e4);

% get nonlinear evolution of X given our control law, u = Kx
[t, Xref] = ode45(@inverted_pendulum, Tspan, [x0_bar; u0], options, parameters, K3_1);

x_non_lin = Xref(:, 1:4);
u_non_lin = Xref(:, 5);

figure
plot(t, x_non_lin)
