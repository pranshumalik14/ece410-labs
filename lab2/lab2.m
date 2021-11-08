%%  Section 3.1 (Basic operations on subspaces)

clear all;

%% Output 1

fprintf('######################### OUTPUT 1 #########################\n');
V = [1 1 3; -1 1 1; 0 0 0; 1 0 1];
W = [1 1; 0 0; 2 -2; 1 0];

basis_V = orth(V)

basis_W = orth(W)

% W was already lin independent since dim(W) = rank(W) = 2
% V was not lin independent since dim(V) = 3 but rank(V) = 2 

% The column space in V was not linearly independent; the rank of V was 2 < the columns of V. 
% This led to the orth basis of V containing only 2 vectors instead of 3.
% The column space in W was linearly independent; the orth basis contained 
% the same number of column vectors as itself.
% One interesting point to note is that W contained [1 0 -2 0]T but in the orth basis, 
% the last 0 did not show up. This is because orth returns vectors that are not only 
% linearly independent but also orthogonal and normalized.

%%  Output 2

fprintf('######################### OUTPUT 2 #########################\n');
fprintf('V + W spans:\n')
span_sum = sub_sum(V, W)

fprintf('V \cap W spans:\n')
span_intersect = sub_intersect(V, W)

fprintf(['To ensure the given calculations are right, we could check whether', ...
    'span_intersect is in span_sum, which should be true for any subspaces,', ...
    'i.e., V \cap W must be in the image of V+W.\n']);
rank_sum = rank(span_sum)
rank_intersect_sum = rank([span_sum span_intersect])

%%  Section 3.2 (Linear transformations and linear equations)

clear all;

%%  Output 3 

fprintf('######################### OUTPUT 3 #########################\n');
x1 = [1;1];
x2 = [2;1];

% The basis covered by x1 x2 will just be a matrix P that is a
% concatenation of both x1, x2
P = [x1 x2];
% x = Pz hence z = inv(P)*x, where z are the coor

x = [2;1];
fprintf(['The new coordinates (z) for (x) are determined using the equation ', ...
    'x = Pz.\nIn this case, P and x are known, so solving for z is straightforward.\n']);
z = P\x

fprintf(['The input x is represented as x = 0 * x1 + 1 * x2 in the new coordinate system.', ...
    '\nThis makes sense as x = x2.\n']);

z1 = z(1);
z2 = z(2);
fprintf(['Further multiplying z1*x1 + z2*x2 should transform z back into x', ...
    '(in the original coordinates),\nwhere x1 and x2 are the basis vectors', ...
    'established at the start of the question.\n']);

numeric_x = x1*z1 + x2*z2

%%  Output 4 

fprintf('######################### OUTPUT 4 #########################\n')
A = [1 2 0 -1; 0 1 -1 -2; 1 0 3 4; 0 -1 2 3; 0 0 2 2];

P = [1 1 1 1; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Q = [1 1 0 0 1; 1 -1 0 0 0; 0 0 1 1 0; 0 0 1 -1 0; 0 0 0 0 1];

% Each column of AP represents the transformed P
% Basically, we have Qj = A_hat * (A*xj),
% where xj is the jth column of P
% and Qj is the jth column of Q
% We can compute all columns of A_hat at once
% by first transforming all of P into AP (A * P)
% Then solve Q = A_hat * AP (Q and AP are known)

AP = A*P
%this returns the all A*xj where xj is the jth column of P and j goes from 1 to n

A_hat = Q\AP

% To test A_hat, can perform the conceptual procedure provided above for
% one column of P and verify that it matches the corresponding column in
% A_hat computed above
p1 = P(:, 1)

Ap1 = A*p1 % This is the first column of P transformed by A

% This is the jth column of A_hat == this is expressing the Ap1 in the coordinate 
% basis of Q (hence A_hat1). This matches the first column of A_hat computed above
A_hat1 = Q\Ap1

%%  Output 5

fprintf('######################### OUTPUT 5 #########################\n')
% The linear transformation by A maps R_4 to R_5
% dim(A) will be the number of points needed to describe any point in A
% hence dim(A) can be # of columns (WLOG)
dim_A = size(A,2)

dim_img_A = rank(A) 
% recall rank = dim(img(A))

% rank-nullity theorem states dimA = rank(A) + null(A)
% hence using the theorem
% recall nullility(A) = dim(nullspace(A)) or dim(kernel(A))
dim_null_A = dim_A - dim_img_A

% This can be verified by directly computing the nullity (dim of nullspace/kernel)
null_space_dim = size(null(A),2) % need columns

% The transformation is not surjective because the output value is in R_5
% but the rank is 3 => A cannot span all of R_5 with just 3 LI column vectors
% Similarly the transformation is not injective as the rank(A) is 3. The row 
% space of A cannot span all of R_4 as there are only 3 LI row vectors.

% use the generic matrix mapping characterization function
[injective, surjective] = characterize_matrix_mapping(A);

% print result
fprintf('A is ');
if (injective && surjective)
    fprintf('a bijective map\n');
elseif (injective && ~surjective)
    fprintf('an injective but not a surjective map\n');
elseif (~injective && surjective)
    fprintf('not an injective but a surjective map\n');
else
    fprintf('neither an injective nor a surjective map\n');
end

%%  Output 6  

fprintf('######################### OUTPUT 6 #########################\n')
rank_A = rank(A)

fprintf(['To determine whether a vector b has a solution, multiple solutions', ...
    'or no solutions to Ax = b,\n we must determine the rank of [A|b].\n', ...
    'If the rank is unaffected then b is in the image of A, hence rank([A|b]) = rank(A)\n']);
% Must check rank of [A | b] 
% If rank(A) = rank([A | b]) => b will is in the image of the column space of A
% Essentially, an unaffected rank will imply b is in the column space
% But if the rank increases, then Ax = b does not have a solution

b1 = [1; 0; 0; 0; 0];
rank_A_b1 = rank([A b1])
% ASK PRANSHU: how can we determine whether there's one or multiple
% solutions?
% if rank_A_b1 == rank_A:
%     fprintf('\nRank [A | b] == rank(A). This implies there is at least one solution')
% rank increases, so no solution

b2 = [1; 1; -2; -2; -2];
rank_A_b2 = rank([A b2])
% rank unaffected here -> there exists a solution for A x = b2

x = A\b2
% Now, the nullity(A) is not 0 so to any solution for A*x = b2
% x1 can be defined as x1 = x + k*y, where y is a vector in the nullspace
% and k is any real number
% null(A) indicates there is 1 vector in the nullspace of A (nullity = 1)

% Consider x' = x + null(A)
x1 = x + null(A)

% Can check if Ax1 = b2
A*x1
% This is equal to b2. Hence, x1 is also a valid solution

%% Section 4 (A-invariance and Representation Theorem)

clear all;

%% Output 7

fprintf('######################### OUTPUT 7 #########################')
A = [1 2 2; 1 -3 -5; -1 2 4]
V = [0 1; 1 -2; -1 1]

% Must check Rank([AV | V]) = Rank(V)
fprintf('\nRank of [AV | V]\n')
LHS = rank([A*V V])

fprintf('\nRank of V\n')
RHS = rank(V)

if RHS == LHS
    fprintf('Both ranks are equal. Hence, V is A-invariant\n')
else
    fprintf('Ranks are unequal. Hence, V is not A-invariant\n')
end
% LHS == RHS => V is A-invariant

% To find P, need to first determine W, that is the complement subspace of V
W = null(V');

fprintf('\nMatrix P of the representation theorem:\n')
P = [V W]

% A_hat = inv(P) * A * P
fprintf('\nA_hat = inv(P) * A * P =\n')
A_hat = P\(A*P)

% Both yield the same result, the second option is computationally faster
% A_hat is block upper triangular.

%% Section 5 (Controllability and Kalman Decomposition)

clear all;

%% Output 8

fprintf('######################### OUTPUT 8 #########################')
A = [5 -11 5; 0 -6 0; -5 5 -1];
B = [1 -2; 0 0; 1 2];

% To check for controllability, must determine if rank(Qc) = n (3 in this case)
n = size(A, 2);
k = rank(ctrb(A, B));
% Rank here is 2 (which means system is not controllable!)  

% Check kalman_decomposition.m for more calculations and comments
[A_hat, B_hat] = kalman_decomposition(A, B);

% Rounding off values in A_hat, B_hat
A_hat = round(A_hat, 4)
B_hat = round(B_hat, 4)

% correspond to the controllable subsystem
sympref('FloatingPointOutput',true);
syms z1 z2 z3 u1 u2;
syms z1_dot z2_dot z3_dot

controllable_system = A_hat(1:k, 1:k )* [z1;z2] + A_hat(1:k, k+1: n) * z3 + B_hat(1:k, :)*[u1;u2];
fprintf('\nThe following are differential equations for the controlloable systems (z1_dot, z2_dot):')
z1_dot = controllable_system(1)
z2_dot = controllable_system(2)

% The last row corresponds to the uncontrollable subsystem
fprintf('\nThe following are differential equations for the uncontrollable system (z3_dot)\n')
z3_dot = A_hat(k+1:n, 1:k)* [z1; z2] + A_hat(k+1:n, k+1:n)*z3
