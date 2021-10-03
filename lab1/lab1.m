%% ece410: linear control systems
%  lab1: numerical simulation of dynamical systems and symbolic linearization
%  authors: Pranshu Malik and Varun Sampat
%  date: 3 October 2021

clc;
close all;
clear all;

%% numerical simulation

% set up <>
parameters = struct('M', 0.2, 'g', 9.81, 'l', 0.15);
x_0        = [0                               0; 
              sqrt(parameters.g/parameters.l) 0.199*sqrt(parameters.g/parameters.l)];
options    = odeset('RelTol',1e-7,'AbsTol',1e-7);
Tspan      = linspace(0,10,1e3);

% numerical integration: calculating evolution of system states
[t_1,x_1t] = ode45(@pendulum, Tspan, x_0(:,1), options, parameters);
[t_2,x_2t] = ode45(@pendulum, Tspan, x_0(:,2), options, parameters);

% plotting (todo: fix labels, aspect ratios, export to files)
x_1t_1 = x_1t(:, 1); x_1t_2 = x_1t(:, 2);
x_2t_1 = x_2t(:, 1); x_2t_2 = x_2t(:, 2);

f_1_1 = figure('Name','Individual State Evolution (starting x_0_1)','NumberTitle','off');
f_1_2 = figure('Name','State Orbit (starting x_0_1)','NumberTitle','off');

figure(f_1_1);
subplot(2,1,1);
plot(t_1, x_1t_1)
xlabel('Time (seconds)')
ylabel('$\theta$', 'Interpreter','latex')
subplot(212);
plot(t_1, x_1t_2)
xlabel('Time (seconds)')
ylabel('$\dot{\theta}$', 'Interpreter','latex')

figure(f_1_2);
plot(x_1t_1, x_1t_2)
% todo: control aspect ratio of this plot
xlabel('x_1')
ylabel('x_2')

f_2_1 = figure('Name','Individual State Evolution (starting x_0_2)','NumberTitle','off');
f_2_2 = figure('Name','State Orbit (starting x_0_2)','NumberTitle','off');

figure(f_2_1);
subplot(2,1,1);
plot(t_2, x_2t_1)
xlabel('Time (seconds)')
ylabel('$\theta$', 'Interpreter','latex')
subplot(212);
plot(t_2, x_2t_2)
xlabel('Time (seconds)')
ylabel('$\dot{\theta}$', 'Interpreter','latex')

figure(f_2_2);
plot(x_2t_1, x_2t_2)
% todo: control aspect ratio of this plot
xlabel('x_1')
ylabel('x_2')

%% symbolic linearization

syms x1 x2 t u M l g

x    = [x1; x2];
xdot = [x(2); -g/l*sin(x(1)) - 1/(M*l)*cos(x(1))*u]; % = f(x, u)
y    = x1;

symA = jacobian(xdot,x);
symB = jacobian(xdot,u);
symC = jacobian(y,x);
symD = jacobian(y,u);

% linearization equilibrium
xstar = [0; 0];
ustar = 0;

% todo: fix this and use the xstar and ustar variables
A = subs(symA, {x1, x2, u}, {0, 0, 0});
B = subs(symB, {x1, x2, u}, {0, 0, 0});
C = subs(symC, {x1, x2, u}, {0, 0, 0});
D = subs(symD, {x1, x2, u}, {0, 0, 0});

%% symbolic expression to numerical integration

