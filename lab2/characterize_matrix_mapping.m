% returns whether the matrix A is injective, surjective, neither, or both.
% note that the return variables are boolean variables
function [injective, surjective] = characterize_matrix_mapping(A)

[m, n] = size(A);
rankA  = rank(A);

% if the A does not have full rank then it is neither injective nor surjective
if rankA < min(m, n)
    injective  = false;
    surjective = false;
elseif (m > n) && (n == rankA)  % note that in this case dim(ker(A)) = 0
    injective   = true;
    surjective  = false;
elseif (n > m) && (m == rankA)  % full row rank ==> surjectivity
    injective   = false;
    surjective  = true;
elseif (m == n) && (n == rankA) % both surjective and injective (= bijective mapping)
    injective  = true;
    surjective = true;
else                            % deafult case: should not enter
    injective  = false;
    surjective = false;
end

end

