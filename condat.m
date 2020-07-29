function [u, relnorm] = condat(form, model, algo)
% CONDAT A Generic Proximal Algorithm for Convex Optimization applied to the
% task of general audio reconstruction.
%
% The algorithm solves the problem of minimizing l1 norm of signal
% coefficients while both the signal and the coefficients lie in box-type
% sets.
%
% Although it uses the notation of frame analysis and synthesis, the
% operators may in fact be arbitrary, as long as model.ana is the adjoint
% operator of model.syn and vice versa.
%
% Input arguments
%       form     'analysis' or 'synthesis'
%       model    structure specifying the model
%                - model.ana             frame analysis
%                - model.syn             frame synthesis
%                - model.sparse          proximal operator of the
%                                        sparsity-inducing penalty
%                - model.projT           projection onto the T domain
%                                        feasible set
%                - model.projTF          projection onto the TF domain
%                                        feasible set
%                - model.dim = [ M N ]   dimensions of T and TF domains
%       algo     structure specifying the parameters of the algorithm
%                - algo.tau
%                - algo.sigma
%                - algo.rho
%                - algo.maxit
%                - algo.tol
%
% Output arguments
%       u        output signal / coefficients
%       relnorm  vector maxit x 1 containing the relative norm of the
%                solution during iterations
%
% Date: 16/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% initialization
if strcmp(form, 'analysis')
    analysisform = true;
else
    analysisform = false;
end

M = model.dim(1); % dimension of the T domain
N = model.dim(2); % dimension of the TF domain

if analysisform
    u = zeros(M, 1);
else
    u = zeros(N, 1);
end
v1 = zeros(N, 1);
v2 = zeros(M, 1);
v3 = zeros(N, 1);

relnorm = NaN(algo.maxit, 1);

%% main iteration for the analysis formulation
if analysisform
    for n = 1:algo.maxit
       
        % precomputation
        anau   = model.ana(u);
        
        % update corresponding to the sparse penalty
        v1tilde = v1 + algo.sigma * anau...
                   - algo.sigma * model.sparse(v1/algo.sigma + anau);
        v1new   = algo.rho * v1tilde + (1-algo.rho) * v1;

        % projection in the T domain
        v2tilde = v2 + algo.sigma * u...
                   - algo.sigma * model.projT(v2/algo.sigma + u);
        v2new   = algo.rho * v2tilde + (1-algo.rho) * v2;

        % projection in the TF domain
        v3tilde = v3 + algo.sigma * anau...
                   - algo.sigma * model.projTF(v3/algo.sigma + anau);
        v3new   = algo.rho * v3tilde + (1-algo.rho) * v3;

        % update u
        unew   = u - algo.rho * algo.tau * model.syn(2*v1tilde + 2*v3tilde - v1 - v3)...
                   - algo.rho * algo.tau * (2*v2tilde - v2);              
               
        % compute the relative norm
        relnorm(n) = norm(unew-u)/norm(u);
        
        % rename variables for the subsequent iteration
        u  = unew;
        v1 = v1new;
        v2 = v2new;
        v3 = v3new;
        
        % check the relative norm
        if relnorm(n) < algo.tol
            break
        end  
    end
end

%% main iteration for the synthesis formulation
if ~analysisform
    for n = 1:algo.maxit
       
        % precomputation
        synu   = model.syn(u);

        % update corresponding to the sparse penalty
        v1tilde = v1 + algo.sigma * u...
                   - algo.sigma * model.sparse(v1/algo.sigma + u);
        v1new   = algo.rho * v1tilde + (1-algo.rho) * v1; 

        % projection in the T domain
        v2tilde = v2 + algo.sigma * synu...
                   - algo.sigma * model.projT(v2/algo.sigma + synu);
        v2new   = algo.rho * v2tilde + (1-algo.rho) * v2;

        % projection in the TF domain
        v3tilde = v3 + algo.sigma * u...
                   - algo.sigma * model.projTF(v3/algo.sigma + u);
        v3new   = algo.rho * v3tilde + (1-algo.rho) * v3;

        % update u
        unew   = u - algo.rho * algo.tau * (2*v1tilde + 2*v3tilde - v1 - v3)...
                   - algo.rho * algo.tau * model.ana(2*v2tilde - v2);

        % compute the relative norm
        relnorm(n) = norm(unew-u)/norm(u);
        
        % rename variables for the subsequent iteration
        u  = unew;
        v1 = v1new;
        v2 = v2new;
        v3 = v3new;
        
        % check the relative norm
        if relnorm(n) < algo.tol
            break
        end   
    end
end
    
%% output
relnorm = relnorm(1:n);

end