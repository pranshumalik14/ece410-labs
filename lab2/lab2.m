%%  Output 1 (Basic operations on subspaces)

%% Output 1.1
V = [1 1 3; -1 1 1; 0 0 0; 1 0 1];
W = [1 1; 0 0; 2 -2; 1 0];

basis_V = orth(V)

basis_W = orth(W)

% W was already lin independent since dim(W) = rank(W) = 2
% V was not lin independent since dim(V) = 3 but rank(V) = 2 

% what has happened to W? nothing bro it's fine only <Ask TA>
%%  Output 1.2
span_sum = sub_sum(V, W)
span_intersect = sub_intersect(V, W)

%%  Output 3 (Linear transformations and linear equations)

%%  Output 3.3 
x1 = [1;1];
x2 = [2;1];

P = [x1 x2];
% x = Pz hence z = inv(P)*x

x = [2;1];
z = P\x

z1 = z(1);
z2 = z(2);

numeric_x = x1*z1 + x2*z2
% this should equal x = (2,1)
%%  Output 3.4 <Pranshu>
A = [1 2 0 -1; 0 1 -1 -2; 1 0 3 4; 0 -1 2 3; 0 0 2 2];




%%  Output 3.5 
% The linear transformation by A maps R_4 to R_5
% dim(A) will be the number of points needed to describe any point in A
% hence dim(A) can be # of columns (WLOG)
dim_A = size(A,2)

dim_img_A = rank(A)

% rank-nullity states dimA = rank(A) + null(A)
% hence using the theorem
null_A = dim_img_A - dim_A

% This can be verified by directly computing the nullity (dim of nullspace/kernel)
null_space_dim = size(null(A),2) %need columns

fprintf('\nThe transformation is not surjective as the output value is in R_5\nbut the rank is 3 => A cannot span all of R_5 with just 3 vectors')
fprintf('\nSimilarly, the transformation is not injective as the rank(A) is 3.\n The row space of A cannot span all of R_4 as there are only 3 linearly independent vectors')
%%  Output 3.6  
b1 = [1;0;0;0;0]
b2 = [1;1;-2;-2;-2]


%%  Output 4.7 (A-invariance and Representation Theorem)
A = [1 2 2; 1 -3 -5; -1 2 4]
V = [0 1; 1 -2; -1 1]

% Must check Rank([AV | V]) = Rank(V)
LHS = rank([A*V V])
RHS = rank(V)

W = null(V')
P = [V W]

A_hat = P^-1 * A * P
%%  Output 5.8 (Controllability and Kalman decomposition)
