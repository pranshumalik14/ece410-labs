%% ece410: linear control systems
%  lab1: numerical simulation of dynamical systems and symbolic linearization
%  authors: Pranshu Malik and Varun Sampat
%  date: 3 October 2021

clc;
close all;
clear all;

%% numerical simulation

% set up initial values and DE parameters
parameters = struct('M', 0.2, 'g', 9.81, 'l', 0.15);
x_0        = [0                               0; 
              sqrt(parameters.g/parameters.l) 0.199*sqrt(parameters.g/parameters.l)];
options    = odeset('RelTol',1e-7,'AbsTol',1e-7);
Tspan      = linspace(0,10,1e3);

% numerical integration: calculating evolution of system states
[t_1,x_1t] = ode45(@pendulum, Tspan, x_0(:,1), options, parameters);
[t_2,x_2t] = ode45(@pendulum, Tspan, x_0(:,2), options, parameters);

% plotting
x_1t_1 = x_1t(:, 1); x_1t_2 = x_1t(:, 2);
x_2t_1 = x_2t(:, 1); x_2t_2 = x_2t(:, 2);

f_1_1   = figure('Name','Individual State Evolution (starting x_0_1)','NumberTitle','off');
f_1_2_a = figure('Name','State Orbit (starting x_0_1; axis equal)','NumberTitle','off');
f_1_2_b = figure('Name','State Orbit (starting x_0_1)','NumberTitle','off');

figure(f_1_1);

subplot(2,1,1);
plot(t_1, x_1t_1);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');

subplot(212);
plot(t_1, x_1t_2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(f_1_2_a);
plot(x_1t_2, x_1t_1);
daspect([1 1 1])
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
xlabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(f_1_2_b);
plot(x_1t_1, x_1t_2);
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

f_2_1   = figure('Name','Individual State Evolution (starting x_0_2)','NumberTitle','off');
f_2_2_a = figure('Name','State Orbit (starting x_0_2; axis equal)','NumberTitle','off');
f_2_2_b = figure('Name','State Orbit (starting x_0_2)','NumberTitle','off');

figure(f_2_1);

subplot(2,1,1);
plot(t_2, x_2t_1);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');

subplot(212);
plot(t_2, x_2t_2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(f_2_2_a);
plot(x_2t_2, x_2t_1);
daspect([1 1 1])
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
xlabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(f_2_2_b);
plot(x_2t_1, x_2t_2);
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

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

syms z1 z2

z    = [z1; z2];
zdot = A*(z - xstar) + B*(u - ustar);

Xdot = [xdot; zdot]; % augmented state ODE
Xdot = subs(Xdot, {M, l, g, u}, {parameters.M, parameters.l, parameters.g, 0}); % u = 0

augmented_pend = matlabFunction(Xdot, 'vars', {t, [x;z]});

% numerical integration: calculating evolution of augmented system states
[t_1,X_1t] = ode45(augmented_pend, Tspan, [x_0(:,1); x_0(:,1)], options);
[t_2,X_2t] = ode45(augmented_pend, Tspan, [x_0(:,2); x_0(:,2)], options);

% plotting
x_1t_1 = X_1t(:, 1); x_1t_2 = X_1t(:, 2); z_1t_1 = X_1t(:, 3); z_1t_2 = X_1t(:, 4);
x_2t_1 = X_2t(:, 1); x_2t_2 = X_2t(:, 2); z_2t_1 = X_2t(:, 3); z_2t_2 = X_2t(:, 4);

F_1_1   = figure('Name','Individual State Evolution (starting X_0_1)','NumberTitle','off');
F_1_2_a = figure('Name','State Orbit (starting X_0_1; axis equal)','NumberTitle','off');
F_1_2_b = figure('Name','State Orbit (starting X_0_1)','NumberTitle','off');

figure(F_1_1);

subplot(2,1,1);
plot(t_1, x_1t_1);
hold on;
plot(t_1, z_1t_1);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_x$', '$\theta_z$'), 'Interpreter', 'latex');

subplot(212);
plot(t_1, x_1t_2);
hold on;
plot(t_1, z_1t_2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta_x}$', '$\dot{\theta_z}$'), 'Interpreter', 'latex');

figure(F_1_2_a);
plot(x_1t_2, x_1t_1);
hold on;
plot(z_1t_2, z_1t_1);
daspect([1 1 1])
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
xlabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\mathbf{x}(t)$', '$\mathbf{z}(t)$'), 'Interpreter', 'latex');

figure(F_1_2_b);
plot(x_1t_1, x_1t_2);
hold on;
plot(z_1t_1, z_1t_2);
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\mathbf{x}(t)$', '$\mathbf{z}(t)$'), 'Interpreter', 'latex');

F_2_1 = figure('Name','Individual State Evolution (starting X_0_2)','NumberTitle','off');
F_2_2_a = figure('Name','State Orbit (starting X_0_2; axis equal)','NumberTitle','off');
F_2_2_b = figure('Name','State Orbit (starting X_0_2)','NumberTitle','off');

figure(F_2_1);

subplot(2,1,1);
plot(t_2, x_2t_1);
hold on;
plot(t_2, z_2t_1);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_x$', '$\theta_z$'), 'Interpreter', 'latex');

subplot(212);
plot(t_2, x_2t_2);
hold on;
plot(t_2, z_2t_2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta_x}$', '$\dot{\theta_z}$'), 'Interpreter', 'latex');

figure(F_2_2_a);
plot(x_2t_2, x_2t_1);
hold on;
plot(z_2t_2, z_2t_1);
daspect([1 1 1])
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
xlabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\mathbf{x}(t)$', '$\mathbf{z}(t)$'), 'Interpreter', 'latex');

figure(F_2_2_b);
plot(x_2t_1, x_2t_2);
hold on;
plot(z_2t_1, z_2t_2);
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\mathbf{x}(t)$', '$\mathbf{z}(t)$'), 'Interpreter', 'latex');

%% lti representations

%% autoexport figures to (pdf) files

set(F_1_1,'Units','Inches');
pos = get(F_1_1,'Position');
set(F_1_1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(F_1_1,'/figs/test','-dpdf','-r0')
