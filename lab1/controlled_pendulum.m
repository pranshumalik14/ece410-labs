% pendulum state evolution ODE
function Xdot = controlled_pendulum(t, X, parameters, controller_matrices)

% extract parameters
M = parameters.M;
g = parameters.g;
l = parameters.l;

% controller matrices
F = controller_matrices.F;
G = controller_matrices.G;
H = controller_matrices.H;
L = controller_matrices.L;

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
