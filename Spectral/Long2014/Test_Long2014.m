% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with the Long2014 method

% Set the size of the phantom
ObjectSize = 32;
NbIters = 10;
NbSubsets = 4;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Simulate noisy preps
noise = true;
[y, ~, A, M, S, T, ~] = SimulatePreps(ObjectSize, noise, false, 'none');
S = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature

% Set regularization parameters
lambda = [100, 100, 100];
delta_hyperbola = [0.0001, 0.0001, 0.01];

% - the LongAndFessler2014 method
runOnGPU = false;
tic
[Long2014iterates, costs] = Long2014(y, ObjectSize, NbPixelsPerProj, A, M, S, T, NbIters, lambda, NbSubsets, delta_hyperbola, runOnGPU);
toc

% Plot the resulting costs, in loglog scale
loglog(costs - min(costs(:)))

PlayIterates(Long2014iterates, [0 0 0], [0.015 0.015 1.5])