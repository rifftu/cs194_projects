% [Motion, Shape, T] = factorization(Wo, iterMax1, iterMax2, stopError1, stopError2)
% 
% This function computes the 3D object shape from missing and degenerate data. 
% 
% 
% Input arguments:
% 
%     Wo - The data matrix is defined as:
%                 Wo = [ u_1^1  ...  u_P^1
%                        v_1^1  ...  v_P^1
%                          .    .      .
%                          .     .     .
%                          .      .    .
%                        u_1^F  ...  u_P^F
%                        v_1^F  ...  v_P^F ]
%                 The missing data entries must be identified with NaN.
%                 
%     iterMax1 [Typical values: 3000-15000] - Maximum number of iterations for 
%                                             the missing data algorithm
%     iterMax2 [100-500] - Maximum number of iterations for Rigid Factorization
%     stopError1 [10^(-4) - 10^(-7)] - The missing data algorithm stops if the
%                                      error between two consecutive iterations
%                                      is lower than this threshold
%     stopError2 [10^(-1) - 10^(-3)] - Rigid Factorization algorithm stops if 
%                                      the error between two consecutive iterations 
%                                      is lower than this threshold
% 
% Note that Rigid Factorization is one step of the global missing data algorithm.
% 
% 
% Output arguments:
% 
%     Motion - The motion matrix stacks the camera matrices in all frames. The
%              camera matrix in frame f is a Stiefel matrix and is composed by lines 
%              2*f-1 and 2*f.
%     Shape - 3D object shape
%     T - Translation vector
%     
% 
% For details, see:
% 
%    Marques, M., and Costeira, J.. "Estimating 3D shape from degenerate sequences with missing data", 
%    Computer Vision and Image Understanding, 113(2):261-272, February 2009.
             
function [Motion, T] = getMotionTranslation(Wo, Shape, iterMax1)

W = Wo;

iter1 = 0;

num_kps = size(W, 2);


T = mean(W,2);
offset = T*ones(1, num_kps);

W_centered = W - offset;

R = [1 0 0;0 1 0];

try
while iter1 < iterMax1
    
    
    A_f = projStiefel(R');

    Motion = A_f';

    R = W_centered*pinv(Shape);
% 
    iter1 = iter1 + 1;
% 
    Motionret=Motion;
    Tret=T;

end

catch
    Motion=Motionret;
    T=Tret;
end


end

function W = projStiefel(Wo)
[U,D,V] = svd(Wo,'econ');
c = mean(diag(D));
W = c*U*V';
end