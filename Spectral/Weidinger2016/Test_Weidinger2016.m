% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with the Weidinger2016 method, internally using a basis of
% synthetic materials to speed up convergence

% Set the size of the phantom
ObjectSize = 64;
NbIters = 200;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Simulate noisy preps
noise = true;
[y, ~, A, M, S, T, ~] = SimulatePreps(ObjectSize, noise, false, 'none');
S = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature

% Set regularization parameters
lambda = [1, 1, 1];

% Perform reconstruction
runOnGPU = false;
[Weidinger2016_Regul_iterates, costs ]= Weidinger2016(y, ObjectSize, A, M, S, T, NbIters, lambda, runOnGPU);

% Show the resulting sequence of iterates
PlayIterates(Weidinger2016_Regul_iterates, [0 0 0], [0.015 0.015 1.5])

% Plot the resulting costs, in loglog scale
plot(costs)

% % Save the iterates, if necessary
% save('../Results/Weidinger2016Regul.mat', 'Weidinger2016_Regul_iterates');
