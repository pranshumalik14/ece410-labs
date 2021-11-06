function [A_hat, B_hat] = kalman_decomposition(A, B)

Qc = ctrb(A, B);

V = orth(Qc);
% V is the image of Qc, hence orth() can be used to compute V

W = null(V');
% W is the complement subspace of V

P = [V W];

% AP = P*A_hat
A_hat = P\(A*P);

% B_hat = inv(P)*B
B_hat = P\B;

end