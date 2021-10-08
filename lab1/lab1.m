%% ece410: linear control systems
%  lab1: numerical simulation of dynamical systems and symbolic linearization
%  authors: Pranshu Malik and Varun Sampat
%  date: 3 October 2021

clc;
close all;
clear all;

%% parameters setup

% set up initial values and DE parameters
parameters = struct('M', 0.2, 'g', 9.81, 'l', 0.15);
x_0        = [0                               0; 
              sqrt(parameters.g/parameters.l) 1.99*sqrt(parameters.g/parameters.l)];
% note x_0 here is set as: [x_0_1 | x_0_2] where x_0_1, for instance, is the
% first intial condition

options    = odeset('RelTol',1e-7, 'AbsTol',1e-7);
Tspan      = linspace(0,10,1e3);

%% numerical integration

% numerical integration: calculating evolution of system states
[t_1,x_1t] = ode45(@pendulum, Tspan, x_0(:,1), options, parameters);
[t_2,x_2t] = ode45(@pendulum, Tspan, x_0(:,2), options, parameters);

% plotting vars
x_1t_1 = x_1t(:, 1); x_1t_2 = x_1t(:, 2);
x_2t_1 = x_2t(:, 1); x_2t_2 = x_2t(:, 2);

f_1_1   = figure('Name', 'Individual State Evolution (starting x_0_1)', 'NumberTitle', 'off');
f_1_2_a = figure('Name', 'State Orbit (starting x_0_1; axis equal)', 'NumberTitle', 'off');
f_1_2_b = figure('Name', 'State Orbit (starting x_0_1)', 'NumberTitle', 'off');

figure(f_1_1);

% plotting evolution of x1 over time with sys starting from i.c. x_0_1
subplot(2,1,1);
plot(t_1, x_1t_1);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');

% plotting evolution of x2 over time with sys starting from i.c. x_0_1
subplot(2,1,2);
plot(t_1, x_1t_2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

% plotting state orbit, i.e. x2 vs x1 over time, with sys starting from i.c. x_0_1
figure(f_1_2_a);
plot(x_1t_2, x_1t_1); % flip axes for a manageable plot aspect ratio
daspect([1 1 1])
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
xlabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(f_1_2_b);
plot(x_1t_1, x_1t_2);
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

f_2_1   = figure('Name', 'Individual State Evolution (starting x_0_2)', 'NumberTitle', 'off');
f_2_2_a = figure('Name', 'State Orbit (starting x_0_2; axis equal)', 'NumberTitle', 'off');
f_2_2_b = figure('Name', 'State Orbit (starting x_0_2)', 'NumberTitle', 'off');

figure(f_2_1);

% plotting evolution of x1 over time with sys starting from i.c. x_0_2
subplot(2,1,1);
plot(t_2, x_2t_1);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');

% plotting evolution of x2 over time with sys starting from i.c. x_0_2
subplot(2,1,2);
plot(t_2, x_2t_2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

% plotting state orbit with sys starting from i.c. x_0_2
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
xdot = [x(2); -g/l*sin(x(1)) - 1/(M*l)*cos(x(1))*u]; % xdot = f(x, u)
y    = x1;

symA = jacobian(xdot,x);
symB = jacobian(xdot,u);
symC = jacobian(y,x);
symD = jacobian(y,u);

% linearization equilibrium (for all comparisons going forward)
xstar = [0; 0];
ustar = 0;

% substitute the symbolic variables with the equilibrium points 
A = subs(symA, {x1, x2, u}, {xstar(1), xstar(2), ustar});
B = subs(symB, {x1, x2, u}, {xstar(1), xstar(2), ustar});
C = subs(symC, {x1, x2, u}, {xstar(1), xstar(2), ustar});
D = subs(symD, {x1, x2, u}, {xstar(1), xstar(2), ustar});

% linearization equilibrium for discussion in lab report (not used in any
% subsequent plot comparisons or simulations)
syms theta_star
xstar_theta = [theta_star; 0];
ustar_theta = -M*g*tan(theta_star);
A_thetastar = subs(symA, {x1, x2, u}, {xstar_theta(1), xstar_theta(2), ustar_theta});
B_thetastar = subs(symB, {x1, x2, u}, {xstar_theta(1), xstar_theta(2), ustar_theta});
C_thetastar = subs(symC, {x1, x2, u}, {xstar_theta(1), xstar_theta(2), ustar_theta});
D_thetastar = subs(symD, {x1, x2, u}, {xstar_theta(1), xstar_theta(2), ustar_theta});

%% symbolic expression to numerical integration

syms z1 z2

z    = [z1; z2];
zdot = A*(z - xstar) + B*(u - ustar);

Xdot = [xdot; zdot]; % Xdot = f(X, u), the augmented state ODE
Xdot = subs(Xdot, {M, l, g, u}, {parameters.M, parameters.l, parameters.g, 0}); % u = 0

augmented_pend = matlabFunction(Xdot, 'vars', {t, [x;z]});

% numerical integration: calculating evolution of augmented system states
% starting from the same initial condition x_0 = z_0
[t_1,X_1t] = ode45(augmented_pend, Tspan, [x_0(:,1); x_0(:,1)], options);
[t_2,X_2t] = ode45(augmented_pend, Tspan, [x_0(:,2); x_0(:,2)], options);

% plotting vars
x_1t_1 = X_1t(:, 1); x_1t_2 = X_1t(:, 2); z_1t_1 = X_1t(:, 3); z_1t_2 = X_1t(:, 4);
x_2t_1 = X_2t(:, 1); x_2t_2 = X_2t(:, 2); z_2t_1 = X_2t(:, 3); z_2t_2 = X_2t(:, 4);

F_1_1   = figure('Name', 'Individual State Evolution (starting X_0_1)', 'NumberTitle', 'off');
F_1_2_a = figure('Name', 'State Orbit (starting X_0_1; axis equal)', 'NumberTitle', 'off');
F_1_2_b = figure('Name', 'State Orbit (starting X_0_1)', 'NumberTitle', 'off');

figure(F_1_1);

% plotting evolution of theta over time with sys starting from i.c. X_0_1
subplot(2,1,1);
plot(t_1, x_1t_1);
hold on;
plot(t_1, z_1t_1);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_x$', '$\theta_z$'), 'Interpreter', 'latex');

% plotting evolution of theta dot over time with sys starting from i.c. X_0_1
subplot(2,1,2);
plot(t_1, x_1t_2);
hold on;
plot(t_1, z_1t_2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta_x}$', '$\dot{\theta_z}$'), 'Interpreter', 'latex');

% plotting state orbit, i.e. theta vs theta dot over time, with sys starting from i.c. X_0_1
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

F_2_1 = figure('Name', 'Individual State Evolution (starting X_0_2)', 'NumberTitle', 'off');
F_2_2_a = figure('Name', 'State Orbit (starting X_0_2; axis equal)', 'NumberTitle', 'off');
F_2_2_b = figure('Name', 'State Orbit (starting X_0_2)', 'NumberTitle', 'off');

figure(F_2_1);

% plotting evolution of theta over time with sys starting from i.c. X_0_2
subplot(2,1,1);
plot(t_2, x_2t_1);
hold on;
plot(t_2, z_2t_1);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
set(legend('$\theta_x$', '$\theta_z$'), 'Interpreter', 'latex');

% plotting evolution of theta dot over time with sys starting from i.c. X_0_1
subplot(2,1,2);
plot(t_2, x_2t_2);
hold on;
plot(t_2, z_2t_2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
set(legend('$\dot{\theta_x}$', '$\dot{\theta_z}$'), 'Interpreter', 'latex');

% plotting state orbit, i.e. theta vs theta dot over time, with sys starting from i.c. X_0_2
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

% converting state-space representation matrices into numeric form
numA = double(subs(A, {g, l}, {parameters.g, parameters.l}));
numB = double(subs(B, {M, l}, {parameters.M, parameters.l}));
numC = double(C);
numD = double(D);

sys         = ss(numA, numB, numC, numD); % state-space model of the system
Gs          = tf(sys);     % transfer function of the system
[V,lambda]  = eig(numA);   % eigen- vectors and values of the system matrix A
[Gz,Gp,Gk]  = zpkdata(Gs); % zero-pole data of the plant (pendcart) 

% pole-zero plot of the pendcart system (plant)
f_Gs = pole_zero_plot(sys, 'G(s) pzplot');

%% pendulum stabilization

% controller tf: lead controller with break frequencies at 10 and 1000 and
% gain of -30
Cs   = tf([-3 -30], [1e-5 1]);
f_cs = pole_zero_plot(Cs, 'C(s) pzplot');

% testing stabililty of the closed-loop system using nyquist stability thm.
Ls     = Cs * Gs;
cls_tf = zpk(1 + Ls);

% pole-zero plot of 1 + L(s) for nyquist stability test
f_cls = pole_zero_plot(cls_tf, 'CLS pzplot');

% get controller in ss form for realizing it as (sugmented) system state DE
[F, G, H, L] = ssdata(Cs);

% simulate state evolution of the controlled pendulum (closed loop control)
Tc_span      = linspace(0,2,1e3);
[tc_1,Xc_t1] = ode45(@controlled_pendulum, Tc_span, [x_0(:,1); 0; 0], options, parameters, {F, G, H, L});
[tc_2,Xc_t2] = ode45(@controlled_pendulum, Tc_span, [x_0(:,2); 0; 0], options, parameters, {F, G, H, L});

% plot state evolution diagrams for both initial conditions with z=0 and u=0
F_Xc_1_1 = figure('Name', 'Controlled State Evolution (starting x_0_1)', 'NumberTitle', 'off');
F_Xc_1_2 = figure('Name', 'Controlled State Orbit (starting x_0_1)', 'NumberTitle', 'off');
F_Xc_2_1 = figure('Name', 'Controlled State Evolution (starting x_0_2)', 'NumberTitle', 'off');
F_Xc_2_2 = figure('Name', 'Controlled State Orbit (starting x_0_2)', 'NumberTitle', 'off');

figure(F_Xc_1_1);

% plotting evolution of theta over time with sys starting from i.c. x_0_1; z=0; u=0
subplot(2,1,1);
plot(tc_1, Xc_t1(:,1));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');

% plotting evolution of theta dot over time with sys starting from i.c. x_0_1; z=0; u=0
subplot(2,1,2);
plot(tc_1, Xc_t1(:,2));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

% plotting state orbit, i.e. theta vs theta dot over time, with sys starting from i.c. x_0_1; z=0; u=0
figure(F_Xc_1_2);
plot(Xc_t1(:,1), Xc_t1(:,2));
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(F_Xc_2_1);

% plotting evolution of theta over time with sys starting from i.c. x_0_2; z=0; u=0
subplot(2,1,1);
plot(tc_2, Xc_t2(:,1));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');

% plotting evolution of theta dot over time with sys starting from i.c. x_0_2; z=0; u=0
subplot(2,1,2);
plot(tc_2, Xc_t2(:,2));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

% plotting state orbit, i.e. theta vs theta dot over time, with sys starting from i.c. x_0_2; z=0; u=0
figure(F_Xc_2_2);
plot(Xc_t2(:,1), Xc_t2(:,2));
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

% simulations for initial conditions away from the linearization equilibrium
x_0_tilde = [pi/4 pi/2 pi; 
             0    0    0];
x_0       = [x_0 x_0_tilde]; % append the new initial conditions to the old matrix

[tc_3,Xc_t3] = ode45(@controlled_pendulum, Tc_span, [x_0(:,3); 0; 0], options, parameters, {F, G, H, L});
[tc_4,Xc_t4] = ode45(@controlled_pendulum, Tc_span, [x_0(:,4); 0; 0], options, parameters, {F, G, H, L});
[tc_5,Xc_t5] = ode45(@controlled_pendulum, Tc_span, [x_0(:,5); 0; 0], options, parameters, {F, G, H, L});

% plot solutions
F_Xc_3_1 = figure('Name', 'Controlled State Evolution (starting x_0_3)', 'NumberTitle', 'off');
F_Xc_3_2 = figure('Name', 'Controlled State Orbit (starting x_0_3)', 'NumberTitle', 'off');
F_Xc_4_1 = figure('Name', 'Controlled State Evolution (starting x_0_4)', 'NumberTitle', 'off');
F_Xc_4_2 = figure('Name', 'Controlled State Orbit (starting x_0_4)', 'NumberTitle', 'off');
F_Xc_5_1 = figure('Name', 'Controlled State Evolution (starting x_0_5)', 'NumberTitle', 'off');
F_Xc_5_2 = figure('Name', 'Controlled State Orbit (starting x_0_5)', 'NumberTitle', 'off');

figure(F_Xc_3_1);
subplot(2,1,1);
plot(tc_3, Xc_t3(:,1));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
subplot(2,1,2);
plot(tc_3, Xc_t3(:,2));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(F_Xc_1_2);
plot(Xc_t3(:,1), Xc_t3(:,2));
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(F_Xc_4_1);
subplot(2,1,1);
plot(tc_4, Xc_t4(:,1));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
subplot(2,1,2);
plot(tc_4, Xc_t4(:,2));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(F_Xc_4_2);
plot(Xc_t4(:,1), Xc_t4(:,2));
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(F_Xc_5_1);
subplot(2,1,1);
plot(tc_5, Xc_t5(:,1));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
subplot(2,1,2);
plot(tc_5, Xc_t5(:,2));
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

figure(F_Xc_5_2);
plot(Xc_t5(:,1), Xc_t5(:,2));
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

%% autoexport figures to (pdf) files
%  note: uncomment to save again

% savefig(f_1_1, './figs/section3_x0_1_state_evolution')
% savefig(f_1_2_a, './figs/section3_x0_1_state_orbit_axis_equal')
% savefig(f_1_2_b, './figs/section3_x0_1_state_orbit')
% 
% savefig(f_2_1, './figs/section3_x0_2_state_evolution')
% savefig(f_2_2_a, './figs/section3_x0_2_state_orbit_axis_equal')
% savefig(f_2_2_b, './figs/section3_x0_2_state_orbit')
% 
% savefig(F_1_1, './figs/section5_X0_1_state_evolution')
% savefig(F_1_2_a, './figs/section5_X0_1_state_orbit_axis_equal')
% savefig(F_1_2_b, './figs/section5_X0_1_state_orbit')
% 
% savefig(F_2_1, './figs/section5_X0_2_state_evolution')
% savefig(F_2_2_a, './figs/section5_X0_2_state_orbit_axis_equal')
% savefig(F_2_2_b, './figs/section5_X0_2_state_orbit')

% savefig(f_Gs, './figs/section6_pole_zero_gs')
% savefig(f_cs, './figs/section7_pole_zero_cs')
% savefig(f_cls, './figs/section7_pole_zero_cls')

savefig(F_Xc_1_1, './figs/section7_controlled_state_evolution_x_0_1')
savefig(F_Xc_1_2, './figs/section7_controlled_state_orbit_x_0_1')
savefig(F_Xc_2_1, './figs/section7_controlled_state_evolution_x_0_2')
savefig(F_Xc_2_2, './figs/section7_controlled_state_orbit_x_0_2')
savefig(F_Xc_3_1, './figs/section7_controlled_state_evolution_x_0_3')
savefig(F_Xc_3_2, './figs/section7_controlled_state_orbit_x_0_3')
savefig(F_Xc_4_1, './figs/section7_controlled_state_evolution_x_0_4')
savefig(F_Xc_4_2, './figs/section7_controlled_state_orbit_x_0_4')
savefig(F_Xc_5_1, './figs/section7_controlled_state_evolution_x_0_5')
savefig(F_Xc_5_2, './figs/section7_controlled_state_orbit_x_0_5')
