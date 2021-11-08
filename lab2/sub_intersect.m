function C = sub_intersect(A, B)

[row_A, col_A] = size(A); % size: nxl (ell); basis for R^l

[row_B, col_B] = size(B); % size: nxm; basis for R^m

% x = V * lambda = W * mu
% To be V intersect W
% Where lamda belongs to R_l
% mu belongs to R_m

X = [A -B];

nullspace = null(X); % nullspace of X will return values of lambda and mu 

% First l rows correspond to lambda
lambda = nullspace(1:col_A, :);

% Next m rows correspond to mu
mu = nullspace(col_A + 1: col_B + col_A, :);

C = orth(A*lambda);
% C = orth(B * mu); % will generate the same output

end