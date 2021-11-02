function [phihat, rhohat] = OLSEstPhiRho(y)
%% DESCRIPTION: estimation and prediction in linear trend model
%---INPUT VARIABLE(S)---
%   (1) y: (TxN) data matrix stacking series columnwise
%---OUTPUT VARIABLE(S)---
%   (1) phihat: OLS estimate for phi
%   (2) rhohat: OLS estimate for rho

    % Dimensions
    [T, N] = size(y);

    % Auxiliary matrices
    mZ = [ones(T-2,1) (1:T-2)'];
    mM = eye(T-2) - mZ*((mZ'*mZ)\mZ');

    % Initialize vectors
    pi = NaN(1,N);
    qi = NaN(1,N);
    ri = NaN(1,N);
    vistar = NaN(1,N);
    wistar = NaN(1,N);

    % OLS estimator
    for i = 1:N
        diffyilag = y(2:T-1,i)-y(1:T-2,i);
        pi(i) = y(2:T-1,i)'*mM*y(2:T-1,i);
        qi(i) = y(2:T-1,i)'*mM*diffyilag;
        ri(i) = diffyilag'*mM*diffyilag;
        vistar(i) = y(2:T-1, i)'*mM*y(3:T,i);
        wistar(i) = diffyilag'*mM*y(3:T,i);
    end
    OLSEst = [mean(pi) mean(qi); mean(qi) mean(ri)]\[mean(vistar); mean(wistar)];
    
    % Return output
    phihat = OLSEst(1);
    rhohat = OLSEst(2);
end

