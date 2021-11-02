clear variables; clc; close all;
addpath("./functions")

%% Description: Replication of entries in Table 1
% Current code replicates N=100, T=20, rho=0 and phi=1. The reported
% rejection rate in Table 1 of the paper is 0.064 with our result being
% 0.0642.


%% Simulation Settings
N = 100;
T = 20;
alphai = 0;
betai = 0;
rho = 0;
phi = 1;
sgnlevel = 0.05;
sigma = 1;
Nsim = 5000;
rng(12345)

%% Initializing lists
testStat = NaN(Nsim,1);
reject = NaN(Nsim, 1);

%% Start simulations
fprintf('\nStarting simulations\n');
for simiter = 1:Nsim

    % Report progress
    if mod(simiter, 1E3) == 0
        fprintf('\tIteration %5d out of %5d \n', simiter, Nsim);
    end

    % Generate data
    y = NaN(T,N);
    for i = 1:N
        z_i = zeros(T+2,1);
        z_i(1) = normrnd(0,sigma);
        z_i(2) = phi*z_i(1) + normrnd(0,sigma);
        for titer = 3:(T+2)
            z_i(titer) = phi*z_i(titer-1) + rho*( z_i(titer-1)-z_i(titer-2) ) + normrnd(0,sigma);
        end
        z_i(1:2) = []; % drop first two observations

        % Compute y_it
        y(:,i) = alphai + betai*(1:T)' + z_i;       % Eq. (2.3) in paper 
    end

    % Compute test statistic
    testStat(simiter) = fBlanderDhaeneModel2(y);

    % Hypothesis decision
    reject(simiter) = (testStat(simiter) <= norminv(sgnlevel));
end
fprintf('Simulations finished...\n\n');

%% Print output to screen
fprintf('Simulated rejection rate: %5.4f\n', mean(reject));
