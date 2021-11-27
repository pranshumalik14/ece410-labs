% returns autonomous system (closed loop state feedback) time response
%  appended with u
function [t, Xu] = linearized_cls_response(clsfun, K, Tspan, x0)

% numerically integrate linear sys ODE
[t, X] = ode45(clsfun, Tspan, x0);

% calculate control input
u_t  = X*(K');

Xu = [X, u_t];

end
