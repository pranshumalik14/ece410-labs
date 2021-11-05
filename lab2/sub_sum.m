function C = sub_sum(A, B)

% To find the image of A + B
% Just concatenate A with B and take its image
% Orthonomal basis (X) 
X = [A B];

C = orth(X);

end