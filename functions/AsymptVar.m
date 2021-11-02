function [OmegaTilde] = AsymptVar(y,phiTilde, rhoTilde, BiasPhiTilde, BiasRhoTilde)
%% DESCRIPTION: Asymptotic variance computation
%---INPUT VARIABLE(S)---
%   (1) y: (TxN) data matrix stacking series columnwise
%   (2) phiTilde: 
%   (3) rhoTilde:
%   (4) BiasPhiTilde:
%   (5) BiasRhoTilde:
%---OUTPUT VARIABLE(S)---
%   (1) OmegaTilde: consistent estimate of the asymptotic covariance matrix
%   of the bias-corrected estimates

    % Dimensions
    [T, N] = size(y);

    % Initialize vectors/matrices
    pi = NaN(1,N);
    qi = NaN(1,N);
    ri = NaN(1,N);
    vistar = NaN(1,N);
    wistar = NaN(1,N);
    epsTilde = NaN(T-2,N);
    xiTilde = NaN(2,N);
    sumxidot = zeros(2);

    % Auxiliary matrices
    mZ = [ones(T-2,1) (1:T-2)'];
    mM = eye(T-2) - mZ*((mZ'*mZ)\mZ');

    % Recursions
    for i = 1:N
        yi = y(3:T,i);
        yilag = y(2:T-1,i);
        yitwo = y(1:T-2,i);
        diffyilag = yilag - yitwo;

        pi(i) = y(2:T-1,i)'*mM*y(2:T-1,i);
        qi(i) = y(2:T-1,i)'*mM*diffyilag;
        ri(i) = diffyilag'*mM*diffyilag;

        epsTilde(:,i) = mM*yi - phiTilde*mM*yilag-rhoTilde*mM*diffyilag;
        xiTilde(:,i) = [yilag'*epsTilde(:,i); diffyilag'*epsTilde(:,i)] - [pi(i) qi(i); qi(i) ri(i)]*[BiasPhiTilde; BiasRhoTilde];
        sumxidot = sumxidot + xiTilde(:,i)*xiTilde(:,i)';
    end
    
    % Calculation of asymptotic covariance matrix
    Temp = [mean(pi) mean(qi); mean(qi) mean(ri)];                                      % Matrix of averages
    SigmaTilde = (Temp)\((1/N)*sumxidot)/(Temp);                                        % Asymptotic variance of OLS estimator
    [dbPhiTilde, dbRhoTilde] = FiniteDifDerivativeApprox(rhoTilde, T, 1E-4);            % Finite difference approximation to bias derivatives
    G = [1 -dbPhiTilde/(1+dbRhoTilde); 0 1/(1+dbRhoTilde)];                             % Jacobian of transformation

    % Asymptotic covariance matrix of bias-corrected estimators
    OmegaTilde = G*SigmaTilde*G';
end

%----- FUNCTIONS -----%
function [dbphi, dbrho] = FiniteDifDerivativeApprox(rho, T, dx)
    [bphi1, brho1] = AsymptBiasPhiRho(rho, T);
    [bphi2, brho2] = AsymptBiasPhiRho(rho + dx, T);
    dbphi = (bphi2 - bphi1)/dx;
    dbrho = (brho2 - brho1)/dx;
end

