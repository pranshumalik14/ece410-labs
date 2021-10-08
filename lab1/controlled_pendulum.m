% pendulum state evolution ODE
function Xdot = controlled_pendulum(t, X, {parameters, [F G H L]})

% extract parameters
M = parameters.M;
g = parameters.g;
l = parameters.l;

% extract system states
x1 = X(1); x2=X(2);

% extract controller state
z = X(3);

% set input
u = H*z - L*x1;

% set output: state evolution DE
xdot = [x2; -g/l*sin(x1) - 1/(M*l)*cos(x1)*u];
zdot = F*z - G*x1;
Xdot = [xdot; zdot; u];
