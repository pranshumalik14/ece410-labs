% pendulum state evolution ODE appended with control signal u
function Xdot = inverted_pendulum(t, x, parameters, K)

% extract parameters
M = parameters.M;
g = parameters.g;
l = parameters.l;
m = parameters.m;

% extract states
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

% set input
u = K*x(1:4);

% set output: state evolution DE
x1_dot = x(2);
x2_dot = (-(m*l*sin(x(3))*(x(4))^2) + (m*g*sin(x(3))*cos(x(3))) + u) / (M + m*sin(x(3))^2);
x3_dot = x(4);
x4_dot = (-(m*l*sin(x(3))*cos(x(3))*(x(4))^2) + (m+M)*g*sin(x(3)) + u*cos(x(3))) / (l * (M + m*sin(x(3))^2));

xdot = [x1_dot; x2_dot; x3_dot; x4_dot];
Xdot = [xdot; u];

end
