% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with the Mechlem2017 method, internally using a basis of
% synthetic materials to speed up convergence

% Set the size of the phantom
ObjectSize = 32;
NbIters = 50;
NbSubsets = 4;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Simulate noisy preps
noise = true;
[y, ~, R, M, S, T, ~] = SimulatePreps(ObjectSize, noise, false, 'none');
S = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature

% Set regularization parameters
lambda = [10, 10, 100];
delta_huber = [0.001, 0.001, 0.1];

% With regularization
runOnGPU = true;
[Mechlem2017_iterates, costs]= Mechlem2017(y, ObjectSize, NbPixelsPerProj,R, M, S, T, NbIters, lambda, NbSubsets, delta_huber, runOnGPU);

% Show the resulting sequence of iterates
PlayIterates(Mechlem2017_iterates, [0 0 0], [0.015 0.015 1.5])

% Plot the resulting costs, in loglog scale
plot(costs)

% % Save the iterates, if necessary
% save('../Results/Mechlem2017Regul.mat', 'Mechlem2017_Regul_iterates');
