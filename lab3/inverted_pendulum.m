% pendulum state evolution ODE
function xdot = inverted_pendulum(t, x, K, parameters)


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
disp(K)
u = K*x;

% set output: state evolution DE
x1_dot = x2;
x2_dot = (- (m*l*sin(x3)*(x(4))^2) + (m*g*sin(x3)*cos(x3)) + u)/(M + m*sin(x3)^2);
x3_dot = x4;
x4_dot = (- (m*l*sin(x3)*cos(x3)*(x(4))^2) + (m+M)*g*sin(x3) + u*cos(x3) )/ ( l * (M + m*sin(x3)^2));

xdot = [x1_dot; x2_dot; x3_dot; x4_dot];

