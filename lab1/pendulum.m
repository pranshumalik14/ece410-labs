% pendulum state evolution ODE
function xdot = pendulum(t, x, parameters)

% extract parameters
M = parameters.M;
g = parameters.g;
l = parameters.l;

% extract states
x1 = x(1); x2=x(2);

% set input
u = 0;

% set output: state evolution DE
xdot = [x2; -g/l*sin(x1) - 1/(M*l)*cos(x1)*u];

