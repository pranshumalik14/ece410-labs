%% ece410: linear control systems
%  lab4: output feedback stabilization of a cart-pendulum robot
%  authors: Pranshu Malik and Varun Sampat
%  date: 5 December 2021

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
symA = simplify(jacobian(x_dot,x));
symB = simplify(jacobian(x_dot,u));

% substitute the symbolic variables with the equilibrium points
A = subs(symA, {x1, x2, x3, x4, u}, {xstar(1), xstar(2), xstar(3), xstar(4), u});
B = subs(symB, {x1, x2, x3, x4, u}, {xstar(1), xstar(2), xstar(3), xstar(4), u});

%% converting to numeric state-space representation (with output)

numA = double(subs(A, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l}));
numB = double(subs(B, {g, m, M, l}, {parameters.g, parameters.m, parameters.M, parameters.l}));
numC = [1 0 0 0; 0 0 1 0]; % only linear and angular positions are known by the observer (output)

%% observability and state estimation

Qo = obsv(numA, numC);

% system is observable iff Qo is full rank
assert(rank(Qo) == min(size(Qo)));

% generate L1 (state correction term in observer 1)
p_obs_1 = [-10 -11 -12 -13];
L1      = -1*place(numA', numC', p_obs_1)'

% generate L2 (state correction term in observer 2)
p_obs_2 = [-40 -41 -42 -43];
L2      = -1*place(numA', numC', p_obs_2)'

% generate K (state feedback controller)
p_ctrl = [-5.1 -5.2 -5.3 -5.4];
K      = -1*place(numA, numB, p_ctrl); % we use u = Kx instead of u = -Kx

% simulate cls for Z0 = (-0.5, 0, -pi/4, 0)
Z0 = [-0.5; 0; -pi/4; 0];
X0 = [0; 0; 0; 0];

% setup integration error tol and Tspan
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
Tspan   = linspace(0, 10, 1e3);

% get evolutions for cls (linearized) system state, Z(t), and observer 1 state, X1(t)
[t_lin, ZX1_lin] = ode45(@lin_cls_observer_sep_response, Tspan, [Z0; X0], options, {numA, numB, numC, K, L1});
Z1_lin_t = ZX1_lin(:, 1:4);
X1_lin_t = ZX1_lin(:, 5:8);

% get evolutions for cls (linearized) system state, Z(t), and observer 2 state, X2(t)
[~, ZX2_lin] = ode45(@lin_cls_observer_sep_response, Tspan, [Z0; X0], options, {numA, numB, numC, K, L2});
Z2_lin_t = ZX2_lin(:, 1:4);
X2_lin_t = ZX2_lin(:, 5:8);

% get evolutions for cls (nonlinear) system state, Z(t), and observer 1 state, X1(t)
[t_nlin, ZX1_nlin] = ode45(@nonlin_cls_observer_sep_response, Tspan, [Z0; X0], options, parameters, {numA, numB, numC, K, L1});
Z1_nlin_t = ZX1_nlin(:, 1:4);
X1_nlin_t = ZX1_nlin(:, 5:8);

% get evolutions for cls (nonlinear) system state, Z(t), and observer 2 state, X2(t)
[~, ZX2_nlin] = ode45(@nonlin_cls_observer_sep_response, Tspan, [Z0; X0], options, parameters, {numA, numB, numC, K, L2});
Z2_nlin_t = ZX2_nlin(:, 1:4);
X2_nlin_t = ZX2_nlin(:, 5:8);

%% plot results (Linearized system response)
fig_ZX_lin = figure('Name', 'Linear System State Estimation Error Evolution', 'NumberTitle', 'off');
figure(fig_ZX_lin);

% Plot tilde(z)
subplot(4,1,1);
plot(t_lin, X1_lin_t(:, 1)-Z1_lin_t(:, 1))
hold on;
plot(t_lin, X2_lin_t(:, 1)-Z2_lin_t(:, 1))
ylabel('$\tilde{z}$ [m]', 'Interpreter', 'latex');
xlabel('Time $t$ [s]', 'Interpreter','latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Don't think this legend is visible. Should just keep L_1 and L_2
% set(legend('$\tilde{z}_{L_1}$', '$\tilde{z}_{L_2}$'), 'Interpreter', 'latex');

% Plot tilde(zdot)
subplot(4,1,2);
plot(t_lin, X1_lin_t(:, 2)-Z1_lin_t(:, 2))
hold on;
plot(t_lin, X2_lin_t(:, 2)-Z2_lin_t(:, 2))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{z}}$ [m/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');
% set(legend('$\dot{\tilde{z}}_{L_1}$', '$\dot{\tilde{z}}_{L_2}$'), 'Interpreter', 'latex');

% Plot tilde(theta)
subplot(4,1,3);
plot(t_lin, X1_lin_t(:, 3)-Z1_lin_t(:, 3))
hold on;
plot(t_lin, X2_lin_t(:, 3)-Z2_lin_t(:, 3))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\tilde{\theta}$ [rad]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');
% set(legend('$\tilde{\theta}_{L_1}$', '$\tilde{\theta}_{L_2}$'), 'Interpreter', 'latex');

% Plot tilde(thetadot)
subplot(4,1,4);
plot(t_lin, X1_lin_t(:, 4)-Z1_lin_t(:, 4))
hold on;
plot(t_lin, X2_lin_t(:, 4)-Z2_lin_t(:, 4))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{\theta}}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');
% set(legend('$\dot{\tilde{\theta}}_{L_1}$', '$\dot{\tilde{\theta}}_{L_2}$'), 'Interpreter', 'latex');

%% Plot results (non-linear system response) 
% todo: follow above and just change signals to nonlin
fig_ZX_nlin = figure('Name', 'Nonlinear System State Estimation Error Evolution', 'NumberTitle', 'off');
figure(fig_ZX_nlin);

% Plot tilde(z)
subplot(4,1,1);
plot(t_nlin, X1_nlin_t(:, 1) - Z1_nlin_t(:, 1))
hold on;
plot(t_nlin, X2_nlin_t(:, 1) - Z2_nlin_t(:, 1))
ylabel('$\tilde{z}$ [m]', 'Interpreter', 'latex');
xlabel('Time $t$ [s]', 'Interpreter','latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(zdot)
subplot(4,1,2);
plot(t_nlin, X1_nlin_t(:, 2) - Z1_nlin_t(:, 2))
hold on;
plot(t_nlin, X2_nlin_t(:, 2) - Z2_nlin_t(:, 2))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{z}}$ [m/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(theta)
subplot(4,1,3);
plot(t_nlin, X1_nlin_t(:, 3) - Z1_nlin_t(:, 3))
hold on;
plot(t_nlin, X2_nlin_t(:, 3) - Z2_nlin_t(:, 3))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\tilde{\theta}$ [rad]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(thetadot)
subplot(4,1,4);
plot(t_nlin, X1_nlin_t(:, 4) - Z1_nlin_t(:, 4))
hold on;
plot(t_nlin, X2_nlin_t(:, 4) - Z2_nlin_t(:, 4))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{\theta}}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');
%% state estimation with measurement noise

% create gaussian white noise process, W(t), with mean 'mu' and covariance 'covar'
mu = [0; 0]; covar = [0.005 0; 0 0.001];
L  = chol(covar, 'lower'); % transformation matrix to get covariance of covar from std normal
Ws = (L*randn(2, 1000))';  % samples of the WGN
W  = @(t) interp1(Tspan, Ws, t)'; % interpolate W(t) between samples for t in [0, 10] seconds

% get evolutions for cls (linearized) system state, Z(t), and noisy observer 1 state, XN1(t)
[~, ZXN1_lin] = ode45(@lin_cls_noisy_observer_sep_response, Tspan, [Z0; X0], options, {numA, numB, numC, K, L1}, W);
XN1_lin_t = ZXN1_lin(:, 5:8);
ZN1_lin_t = ZXN1_lin(:, 1:4); % ######################### PRANSHU ######################## why was this not computed??

% get evolutions for cls (linearized) system state, Z(t), and noisy observer 2 state, XN2(t)
[~, ZXN2_lin] = ode45(@lin_cls_noisy_observer_sep_response, Tspan, [Z0; X0], options, {numA, numB, numC, K, L2}, W);
XN2_lin_t = ZXN2_lin(:, 5:8);

%% plot results

fig_ZXN_lin = figure('Name', 'Linear System State Estimation Error Evolution with Noise', 'NumberTitle', 'off');
figure(fig_ZXN_lin);

% plot tilde(z)
subplot(4,1,1);
% plot(t_lin, XN1_lin_t(:, 1)-Z1_lin_t(:, 1))
plot(t_lin, XN1_lin_t(:, 1) - ZN1_lin_t(:, 1)) % ######################### PRANSHU ######################## why was this not used??
hold on;
plot(t_lin, XN2_lin_t(:, 1)-Z2_lin_t(:, 1))
ylabel('$\tilde{y}$ [m]', 'Interpreter', 'latex');
% set(legend('$\tilde{y}_{L_1}$', '$\tilde{y}_{L_2}$'), 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% plot tilde(zdot)
subplot(4,1,2);
plot(t_lin, XN1_lin_t(:, 2)-Z1_lin_t(:, 2))
hold on;
plot(t_lin, XN2_lin_t(:, 2)-Z2_lin_t(:, 2))
ylabel('$\dot{\tilde{y}}$ [m/s]', 'Interpreter', 'latex');
% set(legend('$\dot{\tilde{y}}_{L_1}$', '$\dot{\tilde{y}}_{L_2}$'), 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% plot tilde(theta)
subplot(4,1,3);
plot(t_lin, XN1_lin_t(:, 3)-Z1_lin_t(:, 3))
hold on;
plot(t_lin, XN2_lin_t(:, 3)-Z2_lin_t(:, 3))
ylabel('$\tilde{\theta}$ [rad]', 'Interpreter', 'latex');
% set(legend('$\tilde{\theta}_{L_1}$', '$\tilde{\theta}_{L_2}$'), 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% plot tilde(theta dot)
subplot(4,1,4);
plot(t_lin, XN1_lin_t(:, 4)-Z1_lin_t(:, 4))
hold on;
plot(t_lin, XN2_lin_t(:, 4)-Z2_lin_t(:, 4))
ylabel('$\dot{\tilde{\theta}}$ [rad/s]', 'Interpreter', 'latex');
% set(legend('$\dot{\tilde{\theta}}_{L_1}$', '$\dot{\tilde{\theta}}_{L_2}$'), 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

%% calculate mean squared error (MSE) over last half of Tspan
Xtilde1 = XN1_lin_t-ZXN1_lin(:, 1:4);
Xtilde2 = XN2_lin_t-ZXN2_lin(:, 1:4);
n    = ceil(length(Tspan)/2);
mse1 = 1/n*sum(Xtilde1(n+1:end, :).^2);
mse2 = 1/n*sum(Xtilde2(n+1:end, :).^2);

%%  Noiseless output feedback control

% Calculations
% get evolutions for cls (linearized) system state, Z(t), and observer 1 state, X1(t)
[t_lin_feed, ZX1_lin_feed] = ode45(@lin_cls_observer_feedback_response, Tspan, [Z0; X0], options, {numA, numB, numC, K, L1});
Z1_lin_feed_t = ZX1_lin_feed(:, 1:4);
X1_lin_feed_t = ZX1_lin_feed(:, 5:8);

% get evolutions for cls (linearized) system state, Z(t), and observer 2 state, X2(t)
[~, ZX2_lin_feed] = ode45(@lin_cls_observer_feedback_response, Tspan, [Z0; X0], options, {numA, numB, numC, K, L2});
Z2_lin_feed_t = ZX2_lin_feed(:, 1:4);
X2_lin_feed_t = ZX2_lin_feed(:, 5:8);

% get evolutions for cls (nonlinear) system state, Z(t), and observer 1 state, X1(t)
[t_nlin_feed, ZX1_nlin_feed] = ode45(@nonlin_cls_observer_feedback_response, Tspan, [Z0; X0], options, parameters, {numA, numB, numC, K, L1});
Z1_nlin_feed_t = ZX1_nlin_feed(:, 1:4);
X1_nlin_feed_t = ZX1_nlin_feed(:, 5:8);

% get evolutions for cls (nonlinear) system state, Z(t), and observer 2 state, X2(t)
[~, ZX2_nlin_feed] = ode45(@nonlin_cls_observer_feedback_response, Tspan, [Z0; X0], options, parameters, {numA, numB, numC, K, L2});
Z2_nlin_feed_t = ZX2_nlin_feed(:, 1:4);
X2_nlin_feed_t = ZX2_nlin_feed(:, 5:8);

%% Plot results (linear system)
fig_ZX_lin_feed = figure('Name', 'Linear System State Evolution with feedback control', 'NumberTitle', 'off');
figure(fig_ZX_lin_feed);

% Plot tilde(z)
subplot(4,1,1);
plot(t_lin_feed, X1_lin_feed_t(:, 1)-Z1_lin_feed_t(:, 1))
hold on;
plot(t_lin_feed, X2_lin_feed_t(:, 1)-Z2_lin_feed_t(:, 1))
ylabel('$\tilde{z}$ [m]', 'Interpreter', 'latex');
xlabel('Time $t$ [s]', 'Interpreter','latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(zdot)
subplot(4,1,2);
plot(t_lin_feed, X1_lin_feed_t(:, 2)-Z1_lin_feed_t(:, 2))
hold on;
plot(t_lin_feed, X2_lin_feed_t(:, 2)-Z2_lin_feed_t(:, 2))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{z}}$ [m/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(theta)
subplot(4,1,3);
plot(t_lin_feed, X1_lin_feed_t(:, 3)-Z1_lin_feed_t(:, 3))
hold on;
plot(t_lin_feed, X2_lin_feed_t(:, 3)-Z2_lin_feed_t(:, 3))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\tilde{\theta}$ [rad]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(thetadot)
subplot(4,1,4);
plot(t_lin_feed, X1_lin_feed_t(:, 4)-Z1_lin_feed_t(:, 4))
hold on;
plot(t_lin_feed, X2_lin_feed_t(:, 4)-Z2_lin_feed_t(:, 4))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{\theta}}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

%% Plot results (nlinear system)
fig_ZX_nlin_feed = figure('Name', 'Nonlinear System State Evolution with feedback control', 'NumberTitle', 'off');
figure(fig_ZX_nlin_feed);

% Plot tilde(z)
subplot(4,1,1);
plot(t_nlin_feed, X1_nlin_feed_t(:, 1)-Z1_nlin_feed_t(:, 1))
hold on;
plot(t_nlin_feed, X2_nlin_feed_t(:, 1)-Z2_nlin_feed_t(:, 1))
ylabel('$\tilde{z}$ [m]', 'Interpreter', 'latex');
xlabel('Time $t$ [s]', 'Interpreter','latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(zdot)
subplot(4,1,2);
plot(t_nlin_feed, X1_nlin_feed_t(:, 2)-Z1_nlin_feed_t(:, 2))
hold on;
plot(t_nlin_feed, X2_nlin_feed_t(:, 2)-Z2_nlin_feed_t(:, 2))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{z}}$ [m/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(theta)
subplot(4,1,3);
plot(t_nlin_feed, X1_nlin_feed_t(:, 3)-Z1_nlin_feed_t(:, 3))
hold on;
plot(t_nlin_feed, X2_nlin_feed_t(:, 3)-Z2_nlin_feed_t(:, 3))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\tilde{\theta}$ [rad]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(thetadot)
subplot(4,1,4);
plot(t_nlin_feed, X1_nlin_feed_t(:, 4)-Z1_nlin_feed_t(:, 4))
hold on;
plot(t_nlin_feed, X2_nlin_feed_t(:, 4)-Z2_nlin_feed_t(:, 4))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{\theta}}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

%% Output feedback control with measurement noise (non-linear system)
% get evolutions for cls (linearized) system state, Z(t), and noisy observer 1 state, XN1(t)
[t_nlin_noise_feed, ZXN1_nlin] = ode45(@nonlin_cls_noisy_observer_feedback_response, Tspan, [Z0; X0], options, parameters, {numA, numB, numC, K, L1}, W);
XN1_nlin_t = ZXN1_nlin(:, 5:8);
ZN1_nlin_t = ZXN1_nlin(:, 1:4);

% get evolutions for cls (linearized) system state, Z(t), and noisy observer 2 state, XN2(t)
[~, ZXN2_nlin] = ode45(@nonlin_cls_noisy_observer_feedback_response, Tspan, [Z0; X0], options, parameters, {numA, numB, numC, K, L2}, W);
XN2_nlin_t = ZXN2_nlin(:, 5:8);
ZN2_nlin_t = ZXN2_nlin(:, 1:4);

%% Plot results (nlinear system w noise)
fig_ZX_nlin_noise_feed = figure('Name', 'Noisy nonlinear System State Evolution with feedback control', 'NumberTitle', 'off');
figure(fig_ZX_nlin_noise_feed);

% Plot tilde(z)
subplot(4,1,1);
plot(t_nlin_noise_feed, XN1_nlin_t(:, 1)-ZN1_nlin_t(:, 1))
hold on;
plot(t_nlin_noise_feed, XN2_nlin_t(:, 1)-ZN2_nlin_t(:, 1))
ylabel('$\tilde{z}$ [m]', 'Interpreter', 'latex');
xlabel('Time $t$ [s]', 'Interpreter','latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(zdot)
subplot(4,1,2);
plot(t_nlin_noise_feed, XN1_nlin_t(:, 2)-ZN1_nlin_t(:, 2))
hold on;
plot(t_nlin_noise_feed, XN2_nlin_t(:, 2)-ZN2_nlin_t(:, 2))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{z}}$ [m/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(theta)
subplot(4,1,3);
plot(t_nlin_feed, XN1_nlin_t(:, 3)-ZN1_nlin_t(:, 3))
hold on;
plot(t_nlin_feed, XN2_nlin_t(:, 3)-ZN2_nlin_t(:, 3))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\tilde{\theta}$ [rad]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

% Plot tilde(thetadot)
subplot(4,1,4);
plot(t_nlin_noise_feed, XN1_nlin_t(:, 4)-ZN1_nlin_t(:, 4))
hold on;
plot(t_nlin_noise_feed, XN2_nlin_t(:, 4)-ZN2_nlin_t(:, 4))
xlabel('Time $t$ [s]', 'Interpreter','latex');
ylabel('$\dot{\tilde{\theta}}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$L_1$', '$L_2$'), 'Interpreter', 'latex');

%% autoexport figures to (pdf) files
%  note: uncomment to save again

% savefig(fig_ZX_lin, './figs/lin_noiseless_state_est_error')
% savefig(fig_ZX_nlin, './figs/nlin_noiseless_state_est_error')
% savefig(fig_ZXN_lin, './figs/lin_noisy_state_est_error')
% savefig(fig_ZX_lin_feed, './figs/lin_noiseless_state_est_error_feedback')
% savefig(fig_ZX_nlin_feed, './figs/nlin_noiseless_state_est_error_feedback')
% savefig(fig_ZX_nlin_noise_feed, './figs/nlin_noisy_state_est_error_feedback')