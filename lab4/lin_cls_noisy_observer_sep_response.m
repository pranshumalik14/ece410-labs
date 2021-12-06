% closed-loop state feedback control and observer state estimate evolution stochastic ODE
function ZXdot = lin_cls_noisy_observer_sep_response(t, ZX, ssmatrices, N)

% controller matrices
A = ssmatrices{1};
B = ssmatrices{2};
C = ssmatrices{3};
K = ssmatrices{4};
L = ssmatrices{5};

% cls and observer states
Z = ZX(1:4);
X = ZX(5:8);

Zdot = (A + B*K)*Z;
Y    = C*Z + N(t); % added noise
Xdot = (A + L*C)*X + B*K*Z - L*Y;

ZXdot = [Zdot; Xdot]; % augmented controller-observer-system state

end
