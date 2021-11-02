function [testStat] = fBlanderDhaeneModel2(y)
%% DESCRIPTION: Unit root test statistic, cf. Theorem 4.1 in De Blander and Dhaene (2012)
%---INPUT VARIABLE(S)---
%   (1) y: (TxN) data matrix stacking series columnwise
%---OUTPUT VARIABLE(S)---
%   (1) testStat: the test statistic for the unit root test
    
    % Dimensions
    [T, N] = size(y);

    if T<5
        error('The time dimension is too small to estimate all parameters')
    end

    %% Compute functions for bias correction
    rholist = -0.95:0.05:0.95;
    plimRholistT5 = NaN(length(rholist),1);
    plimRholistTrueT = NaN(length(rholist),1);
    for iter = 1:length(rholist)
        rhoj = rholist(iter);
        % Bias for T = 5
        [~, bRhoT5] = AsymptBiasPhiRho(rhoj, 5);
        plimRholistT5(iter) =  rhoj+bRhoT5;
        % Bias for T
        [~, bRho] = AsymptBiasPhiRho(rhoj, T);
        plimRholistTrueT(iter) = rhoj+bRho;
    end

    %% Compute bias-corrected rho based on first 5 observations
    yT5 = y(1:5,1:N);
    [~,rhoHatT5] = OLSEstPhiRho(yT5);
    % Case 1: rhoHatT5 exceeds plimRholistT5
    %   => select smallest rho value
    if rhoHatT5 > max(plimRholistT5)
        rhoTildeT5 = min(rholist);
    % Case 2: rhoHatT5 less than plimRholistT5
    %   => select largest rho value
    elseif rhoHatT5 < min(plimRholistT5)
        rhoTildeT5 = max(rholist);
    % Case 3: rhoHatT5 falls within the range of plimRholistT5
    else
        rhoTildeT5 = interp1(plimRholistT5, rholist, rhoHatT5);
    end

    %% Full sample estimation
    [phiHat, rhoHat] = OLSEstPhiRho(y);
    % Case 1: rhohat does not exceed values is plimRholist
    %   => continue with cases
    if rhoHat < max(plimRholistTrueT)

        % Case 1a: plimRholist is monotonic
        %   => decide between boundary solution or interior
        if checkMonotonic(plimRholistTrueT)
            
            % rhoHat is less than the minimum value in plimRholist
            %   => check which boundary provides minimum            
            if rhoHat <= min(plimRholistTrueT)
                if T==5
                    rhoTilde = max(rholist);
                else
                    rhoTilde = min(rholist);
                end
            % rhoHat is within the range of the (monotonic) plimRholist
            %   => interpolate
            else
                rhoTilde = interp1(plimRholistTrueT, rholist, rhoHat);
            end
        % Case 1b: plimRholist is not monotonic
        %   => compute rhoTildeA and rhoTildeB and select solution closest
        %   to rhoTildeT5
        else
            [rhoTildeA, rhoTildeB] = rhoTildeAandB(plimRholistTrueT, rholist, rhoHat);
            if abs(rhoTildeA - rhoTildeT5) < abs(rhoTildeB - rhoTildeT5)
                rhoTilde = rhoTildeA;
            else
                rhoTilde = rhoTildeB;
            end
        end
    % Case 2: rhohat exceeds values is plimRholist
    %   => use rhoTilde as computed with T = 5
    else
        rhoTilde = rhoTildeT5;
    end

    % Perform test
    [bPhiTilde, bRhoTilde] = AsymptBiasPhiRho(rhoTilde, T);                     % Asymptotic bias at selected rho
    phiTilde = phiHat - bPhiTilde;                                              % Bias-corrected phi
    OmegaTilde = AsymptVar(y, phiTilde, rhoTilde, bPhiTilde, bRhoTilde);        % Asymptotic covariance matrix of bias-corrected estimator
    testStat = sqrt(N)*(phiTilde-1)/sqrt(OmegaTilde(1,1));                      % Test statistic
end

%----- FUNCTIONS -----%
function [rhoTildeA, rhoTildeB] = rhoTildeAandB(plimRholist, rholist, rhoHat)
    % Find index with maximum element in plimRholist
    [~, maxindex] = max(plimRholist);

    % Case 1: rhoHat less than smallest value in plimRholist 
    %   => select minima on both sides
    if rhoHat < min(plimRholist)  
        rhoTildeA = min(rholist(1:maxindex));
        rhoTildeB = max(rholist(maxindex:length(plimRholist)));
    % Case 2: rhoHat is included in right segment only
    %   => select rhoTildeB solution
    elseif rhoHat < min(plimRholist(1:maxindex)) && rhoHat >= min(plimRholist(maxindex:length(plimRholist)))
        rhoTildeB = interp1(plimRholist(maxindex:length(plimRholist)), rholist(maxindex:length(plimRholist)), rhoHat);
        rhoTildeA = rhoTildeB;
    % Case 3: rhoHat is included in left segment only
    %   => select rhoTildeA solution
    elseif rhoHat < min(plimRholist(maxindex:length(plimRholist))) && rhoHat >= min(plimRholist(1:maxindex))
        rhoTildeA = interp1(plimRholist(1:maxindex), rholist(1:maxindex), rhoHat);
        rhoTildeB = rhoTildeA;
    % Case 4: rhoHat is included in both the right and the left segment
    %   => report both solutions
    else
        rhoTildeA = interp1(plimRholist(1:maxindex), rholist(1:maxindex), rhoHat);
        rhoTildeB = interp1(plimRholist(maxindex:length(plimRholist)), rholist(maxindex:length(plimRholist)), rhoHat);
    end
end

function [monotonic] = checkMonotonic(list)
    monotonic = prod(diff(list) > 0) + prod(diff(list) < 0);
end

