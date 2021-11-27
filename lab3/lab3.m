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
K1 = -1*place(numA, numB, p1) % we use u = Kx instead of u = -Kx

% create linearized cls (closed loop system) for K1
A_cls1 = numA + numB*K1;
cls1   = @(t, x) A_cls1*x;

% simulate for x0 = (-0.5, 0, -pi/4, 0)
x0 = [-0.5; 0; -pi/4; 0];

% setup integration error tol and Tspan
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
Tspan   = linspace(0, 10, 1e3);

% [t1, X1] =  ode45(@inverted_pendulum, Tspan, [x0; K1*x0], options, parameters, K1); % non-linear system sim
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

%[t2, X2] = ode45(@inverted_pendulum, Tspan, [x0; K2*x0], options, parameters, K2); % non-linear system sim
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
% will use the function linearized_cls_response()
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
cls_q1_1 = @(t, x) (numA + numB*K_q1_1)*x;
cls_q1_2 = @(t, x) (numA + numB*K_q1_2)*x;

[t, X_q1_1] = linearized_cls_response(cls_q1_1, K_q1_1, Tspan, x0);
[~, X_q1_2] = linearized_cls_response(cls_q1_2, K_q1_2, Tspan, x0);

% plot results
fig_X_q1 = figure('Name', 'State and Control Evolution for K_q1_1 and K_q1_2', 'NumberTitle', 'off');

figure(fig_X_q1);
subplot(5,1,1);
plot(t, X_q1_1(:, 1))
hold on;
plot(t, X_q1_2(:, 1))
ylabel('$y$ [m]', 'Interpreter', 'latex');
set(legend('$y_{K_1}$', '$y_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,2);
plot(t, X_q1_1(:, 2))
hold on;
plot(t, X_q1_2(:, 2))
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'latex');
set(legend('$\dot{y}_{K_1}$', '$\dot{y}_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,3);
plot(t, X_q1_1(:, 3))
hold on;
plot(t, X_q1_2(:, 3))
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_{K_1}$', '$\theta_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,4);
plot(t, X_q1_1(:, 4))
hold on;
plot(t, X_q1_2(:, 4))
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta}_{K_1}$', '$\dot{\theta}_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,5);
plot(t, X_q1_1(:, 5))
hold on;
plot(t, X_q1_2(:, 5))
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$u$ [N]', 'Interpreter', 'latex');
set(legend('$u_{K_1}$', '$u_{K_2}$'), 'Interpreter', 'latex');

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
cls_q2_1 = @(t, x) (numA + numB*K_q2_1)*x;
cls_q2_2 = @(t, x) (numA + numB*K_q2_2)*x;

[~, X_q2_1] = linearized_cls_response(cls_q2_1, K_q2_1, Tspan, x0);
[~, X_q2_2] = linearized_cls_response(cls_q2_2, K_q2_2, Tspan, x0);

% plot results
fig_X_q2 = figure('Name', 'State and Control Evolution for K_q2_1 and K_q2_2', 'NumberTitle', 'off');

figure(fig_X_q2);
subplot(5,1,1);
plot(t, X_q2_1(:, 1))
hold on;
plot(t, X_q2_2(:, 1))
ylabel('$y$ [m]', 'Interpreter', 'latex');
set(legend('$y_{K_1}$', '$y_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,2);
plot(t, X_q2_1(:, 2))
hold on;
plot(t, X_q2_2(:, 2))
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'latex');
set(legend('$\dot{y}_{K_1}$', '$\dot{y}_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,3);
plot(t, X_q2_1(:, 3))
hold on;
plot(t, X_q2_2(:, 3))
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_{K_1}$', '$\theta_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,4);
plot(t, X_q2_1(:, 4))
hold on;
plot(t, X_q2_2(:, 4))
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta}_{K_1}$', '$\dot{\theta}_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,5);
plot(t, X_q2_1(:, 5))
hold on;
plot(t, X_q2_2(:, 5))
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$u$ [N]', 'Interpreter', 'latex');
set(legend('$u_{K_1}$', '$u_{K_2}$'), 'Interpreter', 'latex');

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
cls_R_1 = @(t, x) (numA + numB*K_R_1)*x;
cls_R_2 = @(t, x) (numA + numB*K_R_2)*x;

[~, X_R_1] = linearized_cls_response(cls_R_1, K_R_1, Tspan, x0);
[~, X_R_2] = linearized_cls_response(cls_R_2, K_R_2, Tspan, x0);

% plot results
fig_X_R = figure('Name', 'State and Control Evolution for K_R_1 and K_R_2', 'NumberTitle', 'off');

figure(fig_X_R);
subplot(5,1,1);
plot(t, X_R_1(:, 1))
hold on;
plot(t, X_R_2(:, 1))
ylabel('$y$ [m]', 'Interpreter', 'latex');
set(legend('$y_{K_1}$', '$y_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,2);
plot(t, X_R_1(:, 2))
hold on;
plot(t, X_R_2(:, 2))
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'latex');
set(legend('$\dot{y}_{K_1}$', '$\dot{y}_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,3);
plot(t, X_R_1(:, 3))
hold on;
plot(t, X_R_2(:, 3))
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_{K_1}$', '$\theta_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,4);
plot(t, X_R_1(:, 4))
hold on;
plot(t, X_R_2(:, 4))
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta}_{K_1}$', '$\dot{\theta}_{K_2}$'), 'Interpreter', 'latex');
subplot(5,1,5);
plot(t, X_R_1(:, 5))
hold on;
plot(t, X_R_2(:, 5))
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$u$ [N]', 'Interpreter', 'latex');
set(legend('$u_{K_1}$', '$u_{K_2}$'), 'Interpreter', 'latex');

%% nonlinear comparison of lqr controller performance

x0 = [-1; 0; pi/4; 0]; % new initial state

% get nonlinear evolution of X given our control law, u = Kx
[t_nl, Xu_nl] = ode45(@inverted_pendulum, Tspan, [x0; K_R_1*x0], options, parameters, K_R_1);
X_nl = Xu_nl(:, 1:4);
u_nl = Xu_nl(:, 5);

% get linearized cls conntrolled state time response for K_R_1
A_cls = numA + numB*K_R_1;
cls   = @(t, x) A_cls*x;

[t_l, Xu_l] = linearized_cls_response(cls, K_R_1, Tspan, x0);
X_l = Xu_l(:, 1:4);
u_l = Xu_l(:, 5);

% plot state and control evolution over time for both linear and non-linear closed loop systems
fig_X_nl_l_comp = figure('Name', 'State and Control Evolution Non-linear Comparison', 'NumberTitle', 'off');

figure(fig_X_nl_l_comp);
subplot(5,1,1);
plot(t_l, X_l(:, 1))
hold on;
plot(t_nl, X_nl(:, 1))
ylabel('$y$ [m]', 'Interpreter', 'latex');
set(legend('$y_{l}$', '$y_{nl}$'), 'Interpreter', 'latex');
subplot(5,1,2);
plot(t_l, X_l(:, 2))
hold on;
plot(t_nl, X_nl(:, 2))
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'latex');
set(legend('$\dot{y}_{l}$', '$\dot{y}_{nl}$'), 'Interpreter', 'latex');
subplot(5,1,3);
plot(t_l, X_l(:, 3))
hold on;
plot(t_nl, X_nl(:, 3))
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_{l}$', '$\theta_{nl}$'), 'Interpreter', 'latex');
subplot(5,1,4);
plot(t_l, X_l(:, 4))
hold on;
plot(t_nl, X_nl(:, 4))
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta}_{l}$', '$\dot{\theta}_{nl}$'), 'Interpreter', 'latex');
subplot(5,1,5);
plot(t_l, u_l)
hold on;
plot(t_nl, u_nl)
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$u$ [N]', 'Interpreter', 'latex');
set(legend('$u_{l}$', '$u_{nl}$'), 'Interpreter', 'latex');

% non-linear closed loop system state evolution for shifted initial positions
x0_shiftpos1  = [-2; 0; pi/4; 0];
x0_shiftpos2  = [-6; 0; pi/4; 0];
x0_shiftpos3  = [-8; 0; pi/4; 0];
x0_shiftpos4  = [-11.8; 0; pi/4; 0];
[t_nl_shift, Xshift_nl1] = ode45(@inverted_pendulum, Tspan, [x0_shiftpos1; K_R_1*x0_shiftpos1], options, parameters, K_R_1);
[~, Xshift_nl2]          = ode45(@inverted_pendulum, Tspan, [x0_shiftpos2; K_R_1*x0_shiftpos2], options, parameters, K_R_1);
[~, Xshift_nl3]          = ode45(@inverted_pendulum, Tspan, [x0_shiftpos3; K_R_1*x0_shiftpos3], options, parameters, K_R_1);
[~, Xshift_nl4]          = ode45(@inverted_pendulum, Tspan, [x0_shiftpos4; K_R_1*x0_shiftpos4], options, parameters, K_R_1);

% plot results
fig_X_nl_shift = figure('Name', 'State and Control Evolution For Shifted Initial Position', 'NumberTitle', 'off');

figure(fig_X_nl_shift);
subplot(5,1,1);
plot(t_nl_shift, Xshift_nl1(:, 1))
hold on;
plot(t_nl_shift, Xshift_nl2(:, 1))
hold on;
plot(t_nl_shift, Xshift_nl3(:, 1))
hold on;
plot(t_nl_shift, Xshift_nl4(:, 1))
ylabel('$y$ [m]', 'Interpreter', 'latex');
set(legend('$y_0 = -2$', '$y_0 = -6$', '$y_0 = -8$', '$y_0 = -11.8$'), 'Interpreter', 'latex');
subplot(5,1,2);
plot(t_nl_shift, Xshift_nl1(:, 2))
hold on;
plot(t_nl_shift, Xshift_nl2(:, 2))
hold on;
plot(t_nl_shift, Xshift_nl3(:, 2))
hold on;
plot(t_nl_shift, Xshift_nl4(:, 2))
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'latex');
set(legend('$y_0 = -2$', '$y_0 = -6$', '$y_0 = -8$', '$y_0 = -11.8$'), 'Interpreter', 'latex');
subplot(5,1,3);
plot(t_nl_shift, Xshift_nl1(:, 3))
hold on;
plot(t_nl_shift, Xshift_nl2(:, 3))
hold on;
plot(t_nl_shift, Xshift_nl3(:, 3))
hold on;
plot(t_nl_shift, Xshift_nl4(:, 3))
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$y_0 = -2$', '$y_0 = -6$', '$y_0 = -8$', '$y_0 = -11.8$'), 'Interpreter', 'latex');
subplot(5,1,4);
plot(t_nl_shift, Xshift_nl1(:, 4))
hold on;
plot(t_nl_shift, Xshift_nl2(:, 4))
hold on;
plot(t_nl_shift, Xshift_nl3(:, 4))
hold on;
plot(t_nl_shift, Xshift_nl4(:, 4))
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$y_0 = -2$', '$y_0 = -6$', '$y_0 = -8$', '$y_0 = -11.8$'), 'Interpreter', 'latex');
subplot(5,1,5);
plot(t_nl_shift, Xshift_nl1(:, 5))
hold on;
plot(t_nl_shift, Xshift_nl2(:, 5))
hold on;
plot(t_nl_shift, Xshift_nl3(:, 5))
hold on;
plot(t_nl_shift, Xshift_nl4(:, 5))
set(legend('$y_0 = -2$', '$y_0 = -6$', '$y_0 = -8$', '$y_0 = -11.8$'), 'Interpreter', 'latex');
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$u$ [N]', 'Interpreter', 'latex');

%% failure in stabilization of nonlinear system by above lqr controller

x0_unstab = [-11.9; 0; pi/4; 0]; % unstable initial condition (near 12)
Tspan     = linspace(0, 4.5, 1e3);
[t_unstab, X_unstab] = ode45(@inverted_pendulum, Tspan, [x0_unstab; K_R_1*x0_unstab], options, parameters, K_R_1);

% plot results
fig_X_unstab = figure('Name', 'State and Control Evolution For Unstable Initial Condition', 'NumberTitle', 'off');

figure(fig_X_unstab);
subplot(5,1,1);
plot(t_unstab, X_unstab(:, 1))
ylabel('$y$ [m]', 'Interpreter', 'latex');
subplot(5,1,2);
plot(t_unstab, X_unstab(:, 2))
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'latex');
subplot(5,1,3);
plot(t_unstab, X_unstab(:, 3))
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
subplot(5,1,4);
plot(t_unstab, X_unstab(:, 4))
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
subplot(5,1,5);
plot(t_unstab, X_unstab(:, 5))
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$u$ [N]', 'Interpreter', 'latex');

%% autoexport figures to (pdf) files
%  note: uncomment to save again

% savefig(fig_X_pp, './figs/pole_place_K1_K2_state_evolutions')
% savefig(fig_X_q1, './figs/q1_K1_K2_state_evolutions')
% savefig(fig_X_q2, './figs/q2_K1_K2_state_evolutions')
% savefig(fig_X_R, './figs/R_K1_K2_state_evolutions')
% savefig(fig_X_nl_l_comp, './figs/lin_nonlin_state_evolutions')
% savefig(fig_X_nl_shift, './figs/nonlin_shiftpos_state_evolutions')
% savefig(fig_X_unstab, './figs/nonlin_unstable_state_evolution')
