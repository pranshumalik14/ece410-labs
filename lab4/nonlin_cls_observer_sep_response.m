% returns pendulum state evolution ODE appended with control signal u
function ZXdot = inverted_pendulum(t, ZX, parameters, ssmatrices)

% extract parameters
M = parameters.M;
g = parameters.g;
l = parameters.l;
m = parameters.m;

% controller matrices
A = ssmatrices{1};
B = ssmatrices{2};
C = ssmatrices{3};
K = ssmatrices{4};
L = ssmatrices{5};

% cls and observer states
Z = ZX(1:4);
X = ZX(5:8);

% set input (u = Kx)
u  = K*Z;

% set output: state evolution DE
z1_dot = Z(2);
z2_dot = (-(m*l*sin(Z(3))*(Z(4))^2) + (m*g*sin(Z(3))*cos(Z(3))) + u) / (M + m*sin(Z(3))^2);
z3_dot = Z(4);
z4_dot = (-(m*l*sin(Z(3))*cos(Z(3))*(Z(4))^2) + (m+M)*g*sin(Z(3)) + u*cos(Z(3))) / (l * (M + m*sin(Z(3))^2));

Zdot = [z1_dot; z2_dot; z3_dot; z4_dot];
Y    = C*Z;
Xdot = (A + L*C)*X + B*K*Z - L*Y;

ZXdot = [Zdot; Xdot]; % augmented controller-observer-system state

end
