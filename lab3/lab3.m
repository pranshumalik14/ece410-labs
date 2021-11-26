%% ece410: linear control systems
%  lab3: state feedback stabilization of a cart-pendulum robot
%  authors: Pranshu Malik and Varun Sampat
%  date: 26 November 2021

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

% output 1: symbolic forms of A and B
symA = simplify(jacobian(x_dot,x))
symB = simplify(jacobian(x_dot,u))

% substitute the symbolic variables with the equilibrium points
A = subs(symA, {x1, x2, x3, x4, u}, {xstar(1), xstar(2), xstar(3), xstar(4), u});
B = subs(symB, {x1, x2, x3, x4, u}, {xstar(1), xstar(2), xstar(3), xstar(4), u});

%% converting to numeric state-space representation

numA = double(subs(A, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l}));
numB = double(subs(B, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l}));

%% controllability and pole-assignment for stabilization

Qc = ctrb(numA, numB);

% system is controllable if Qc is full rank
assert(rank(Qc) == min(size(Qc)));

% generate K1
p1 = [-1 -2 -3 -4];
K1 = -1*place(numA, numB, p1) % use u = Kx instead of u = -Kx

% create linearized cls (closed loop system) for K1
A_cls1 = numA + numB*K1;
cls1   = @(t, x) A_cls1*x;

% simulate for x0 = (-0.5, 0, -pi/4, 0)
x0 = [-0.5; 0; -pi/4; 0];
u0 = K1*x0;

% setup integration error tol and Tspan
options    = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
Tspan      = linspace(0, 10, 1e3);

% [t1, X1] =  ode45(@inverted_pendulum, Tspan, [x0; u0], options, parameters, K1); % non-linear system sim
[t1, X1] = ode45(cls1, Tspan, x0);

% extract state evolutions for K1
x1_t1 = X1(:,1);
x2_t1 = X1(:,2);
x3_t1 = X1(:,3);
x4_t1 = X1(:,4);
% u_t1  = X1(:,5); % from sim for non-linear system

% extract control input, u, for K1
u_t1 = X1*(K1'); % = (K1*X1')'

% generate K2
p2 = [-1 -2 -3 -20];
K2 = -1*place(numA, numB, p2) % use u = Kx instead of u = -Kx

% create linearized cls system for K2
A_cls2 = numA + numB*K2;
cls2   = @(t, x) A_cls2*x;

%[t2, X2] = ode45(@inverted_pendulum, Tspan, [x0; u0], options, parameters, K2); % non-linear system sim
[t2, X2] = ode45(cls2, Tspan, x0);

% extract state evolutions for K2
x1_t2 = X2(:,1);
x2_t2 = X2(:,2);
x3_t2 = X2(:,3);
x4_t2 = X2(:,4);
% u_t2 = X2(:,5); % from sim for non-linear system

% extract control input, u, for K2
u_t2 = X2*(K2'); % = (K2*X2')'

% plot state and control evolution over time for both gains K1 and K2
fig_X_pp = figure('Name', 'State and Control Evolution with Pole Placement (starting x0)', 'NumberTitle', 'off');

figure(fig_X_pp);
subplot(5,1,1);
plot(t1, x1_t1)
hold on;
plot(t2, x1_t2)
ylabel('$y$ [m]', 'Interpreter', 'latex');
set(legend('$y_{K_1}$', '$y_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,2);
plot(t1, x2_t1)
hold on;
plot(t2, x2_t2)
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'latex');
set(legend('$\dot{y}_{K_1}$', '$\dot{y}_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,3);
plot(t1, x3_t1)
hold on;
plot(t2, x3_t2)
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_{K_1}$', '$\theta_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,4);
plot(t1, x4_t1)
hold on;
plot(t2, x4_t2)
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta}_{K_1}$', '$\dot{\theta}_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,5);
plot(t1, u_t1)
hold on;
plot(t2, u_t2)
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$u$ [N]', 'Interpreter', 'latex');
set(legend('$u_{K_1}$', '$u_{K_2}$'), 'Interpreter', 'latex');

%% linear quadratic optimal control for stabilization at equilibrium
% note: now to compute controlled state response (for linearized system, we 
% will use the function plot_linearized_cls_response()
% also note: all simulations are done for the same x0 = (-0.5, 0, -pi/4, 0)

syms q1 q2
symQ = [q1 0 0 0; 0 0 0 0; 0 0 q2 0; 0 0 0 0];

% impact of changing q1 on controlled state time response
R      = 0.5;
num_q2 = 5;

% K for q1 = 0.1
num_q1_1 = 0.1;
numQ_1   = double(subs(symQ, {q1, q2}, {num_q1_1, num_q2}));
K_q1_1   = -1*lqr(numA, numB, numQ_1, R);

% K for q1 = 0.005
num_q1_2 = 0.005;
numQ_2   = double(subs(symQ, {q1, q2}, {num_q1_2, num_q2}));
K_q1_2   = -1*lqr(numA, numB, numQ_2, R);

% time control responses for controllers u = K_q1_1*x and u = K_q1_2*x
fig_q1_1 = figure('Name', 'State and Control Evolution for K_q1_1 (starting x0)', 'NumberTitle', 'off');
fig_q1_2 = figure('Name', 'State and Control Evolution for K_q1_2 (starting x0)', 'NumberTitle', 'off');
cls_q1_1 = @(t, x) (numA + numB*K_q1_1)*x;
cls_q1_2 = @(t, x) (numA + numB*K_q1_2)*x;

[~, X_q1_1] = plot_linearized_cls_response(fig_q1_1, cls_q1_1, K_q1_1, Tspan, x0);
[~, X_q1_2] = plot_linearized_cls_response(fig_q1_2, cls_q1_2, K_q1_2, Tspan, x0);

% impact of changing q2 on controlled state time response
R      = 0.5;
num_q1 = 0.05;

% K for q2 = 1
num_q2_1 = 1;
numQ_1   = double(subs(symQ, {q1, q2}, {num_q1, num_q2_1}));
K_q2_1   = -1*lqr(numA, numB, numQ_1, R);

% K for q2 = 2000
num_q2_2 = 2000;
numQ_2   = double(subs(symQ, {q1, q2}, {num_q1, num_q2_2}));
K_q2_2   = -1*lqr(numA, numB, numQ_2, R);

% time control responses for controllers u = K_q2_1*x and u = K_q2_2*x
fig_q2_1 = figure('Name', 'State and Control Evolution for K_q2_1 (starting x0)', 'NumberTitle', 'off');
fig_q2_2 = figure('Name', 'State and Control Evolution for K_q2_2 (starting x0)', 'NumberTitle', 'off');
cls_q2_1 = @(t, x) (numA + numB*K_q2_1)*x;
cls_q2_2 = @(t, x) (numA + numB*K_q2_2)*x;

[~, X_q2_1] = plot_linearized_cls_response(fig_q2_1, cls_q2_1, K_q2_1, Tspan, x0);
[~, X_q2_2] = plot_linearized_cls_response(fig_q2_2, cls_q2_2, K_q2_2, Tspan, x0);

% impact of changing R
num_q1 = 0.05;
num_q2 = 5;
numQ   = double(subs(symQ, {q1, q2}, {num_q1, num_q2}));

% K for R = 0.005
numR_1 = 0.005;
K_R_1   = -1*lqr(numA, numB, numQ, numR_1);

% K for R = 10
numR_2 = 10;
K_R_2  = -1*lqr(numA, numB, numQ, numR_2);

% time control responses for controllers u = K_R_1*x and u = K_R_2*x
fig_R_1 = figure('Name', 'State and Control Evolution for K_R_1 (starting x0)', 'NumberTitle', 'off');
fig_R_2 = figure('Name', 'State and Control Evolution for K_R_2 (starting x0)', 'NumberTitle', 'off');
cls_R_1 = @(t, x) (numA + numB*K_R_1)*x;
cls_R_2 = @(t, x) (numA + numB*K_R_2)*x;

[~, X_R_1] = plot_linearized_cls_response(fig_R_1, cls_R_1, K_R_1, Tspan, x0);
[~, X_R_2] = plot_linearized_cls_response(fig_R_2, cls_R_2, K_R_2, Tspan, x0);

%% nonlinear comparison of lqr controller performance

x0 = [-1; 0; -pi/4; 0]; % new initial state

% get nonlinear evolution of X given our control law, u = Kx
[t, Xu_nl] = ode45(@inverted_pendulum, Tspan, [x0; u0], options, parameters, K_R_1);
X_nl = Xu_nl(:, 1:4);
u_nl = Xu_nl(:, 5);

% get linearized cls conntrolled state time response for K_R_1
A_cls = numA + numB*K_R_1;
cls   = @(t, x) A_cls*x;

[t2, X_l] = ode45(cls, Tspan, x0);
u_l = X_l*(K3_1');

figure
plot(t, X_l(:,1) , t, X_nl(:,1))

%%

x0_bar = [-175.7; 0; -pi/4; 0];
Tspan  = linspace(0,0.5,1e3);

% get nonlinear evolution of X given our control law, u = Kx
[t, X_nl] = ode45(@inverted_pendulum, Tspan, [x0_bar; u0], options, parameters, K3_1);

X_nl = X_nl(:, 1:4);
u_nl = X_nl(:, 5);

[t, X_l] = ode45(cls, Tspan, x0);
x_lin = X_l(:, 1:4);
u_l = X_l*(K3_1');

figure
plot(t, X_nl(:,1))
hold on
plot(t, x_lin(:,1))
figure
plot(t, X_nl(:,2))
hold on
plot(t, x_lin(:,2))
figure
plot(t, X_nl(:,3))
hold on
plot(t, x_lin(:,3))
figure
plot(t, X_nl(:,4))
hold on
plot(t, x_lin(:,4))
figure
plot(t, u_nl)
hold on
plot(t, u_l)

%% autoexport figures to (pdf) files
%  note: uncomment to save again

% savefig(f_1_1, './figs/section3_x0_1_state_evolution')