function [BiasPhi, BiasRho] = AsymptBiasPhiRho(rho, T)
%% DESCRIPTION: Asymptotic bias, cf. Theorem 3.1 in De Blander and Dhaene (2012)
%---INPUT VARIABLE(S)---
%   (1) rho: autoregressive parameter for error process
%   (2) T: sample size
%---OUTPUT VARIABLE(S)---
%   (1) BiasPhi: asymptotic bias in phi
%   (2) BiasRho: asymptotic bias in rho

    % Matrix definitions
    L = tril(ones(T-2));
    F = zeros(T-2); F(1:end-1,2:end) = eye(T-3);
    Z = [ones(T-2,1) (1:T-2)'];
    M = eye(T-2) - Z*((Z'*Z)\Z');
    A = Amatrix(rho, T);
    D = Dmatrix(rho, T);
    
    % Computations
    MLA = M*L*A;
    MFDtran = M*F*D';
    VectorBias = [trace(MLA*L') trace(MLA); trace(MLA) trace(M*A)]\[trace(MFDtran*L'); trace(MFDtran)];
    BiasPhi = VectorBias(1);
    BiasRho = VectorBias(2);
end

%----- FUNCTIONS -----%
function [mD] = Dmatrix(rho, T)
    MatrixSize = T-2;
    mD = eye(MatrixSize);
    for diagiter = 1:T-3
        % rho^i to fill off-diagonals
        rhoi = rho^diagiter;
        for coliter = 1:(MatrixSize-diagiter)
            mD(coliter+diagiter, coliter) = rhoi;
        end
    end
end

function [mA] = Amatrix(rho, T)
    MatrixSize = T-2;
    mD = Dmatrix(rho,T);
    mA = (mD+mD'-eye(MatrixSize))/(1-rho^2);
end

