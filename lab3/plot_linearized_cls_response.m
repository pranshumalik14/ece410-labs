% plots autonomous system (closed loop state feedback) time response for inverted pendulum
function [t, Xu] = plot_linearized_cls_response(fig, clsfun, K, Tspan, x0)

% numerically integrate linear sys ODE
[t, X] = ode45(clsfun, Tspan, x0);

% extract state responses
x1_t = X(:,1);
x2_t = X(:,2);
x3_t = X(:,3);
x4_t = X(:,4);

% calculate control input
u_t  = X*(K');

% create figure
figure(fig);
subplot(5,1,1);
plot(t, x1_t)
ylabel('$y$ [m]', 'Interpreter', 'latex');
subplot(5,1,2);
plot(t, x2_t)
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'latex');
subplot(5,1,3);
plot(t, x3_t)
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
subplot(5,1,4);
plot(t, x4_t)
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
subplot(5,1,5);
plot(t, u_t)
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Control Input $u$ [N]', 'Interpreter', 'latex');

Xu = [X, u_t];

end

