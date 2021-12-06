% pendulum state evolution ODE
function Xdot = controlled_pendulum(t, X, parameters, ssmatrices)

% extract parameters
M = parameters.M;
g = parameters.g;
l = parameters.l;

% controller matrices
F = ssmatrices{1};
G = ssmatrices{2};
H = ssmatrices{3};
L = ssmatrices{4};

% extract system states
x1 = X(1); % = y (observed output, theta)
x2 = X(2);
z  = X(3); % = controller internal state

% calculate control input
u = H*z - L*x1; % ustar + Hz - L(y-thetastar); thetastar = ustar = 0 while linearizing

% set output: state evolution DE
xdot = [x2; -g/l*sin(x1) - 1/(M*l)*cos(x1)*u];
zdot = F*z - G*x1;         % zdot = F*z ? G*(y-thetastar); thetastar = 0 while linearizing
udot = H*zdot - L*x2;      % udot = d/dt(u)
Xdot = [xdot; zdot; udot]; % augmented controller-system state
