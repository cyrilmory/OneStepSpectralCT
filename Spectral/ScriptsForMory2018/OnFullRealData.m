% Import photon counts
pc = read_mhd('/home/mory/data/tests/rtksparseprojectionmatrixexport/full/photonCounts.mhd');
measuredCounts = pc.data;
[pix_per_proj, rows, projs, bins] = size(measuredCounts);
measuredCounts = reshape(measuredCounts, [pix_per_proj * projs, bins]); % Fails if rows ~= 1

% Import detector response and spectrum model
det = read_mhd('/home/mory/data/Hamburg/realScannerModel/detector_response.mhd');
detectorResponse = det.data';
[detectedEnergies, incidentEnergies] = size(detectorResponse);

sp = read_mhd('/home/mory/data/Hamburg/realScannerModel/spectrumForMatlab.mhd');
incidentSpectrum = squeeze(sp.data).';
% IMPORTANT: rtkprojectionmatrix does not accept projections with non-identity transform. 
% Since PRH projections come in with negative spacing on x, they had to be clitkAffineTransform-ed
% and re-ordered. The spectrum therefore needs to be reordered too
% So this is only for data whose projection matrix is exported via
% rtkprojectionmatrix
incidentSpectrum = fliplr(incidentSpectrum);

mat = read_mhd('/home/mory/data/Hamburg/realScannerModel/material_attenuations.mhd');
M = mat.data';

S_unbinned = detectorResponse .* reshape(incidentSpectrum, [1 incidentEnergies pix_per_proj]);

% Define thresholds and get S binned
S = cat(1, sum(S_unbinned(30:50, :, :), 1), ...
           sum(S_unbinned(51:61, :, :), 1), ...
           sum(S_unbinned(62:71, :, :), 1), ...
           sum(S_unbinned(72:82, :, :), 1), ...
           sum(S_unbinned(83:end, :, :), 1));
       
% No change of basis, use T=identity
T = eye(size(M, 2));
       
% Import sparse geometry matrix
% fileID = fopen('/home/mory/data/tests/rtksparseprojectionmatrixexport/full/fullprojectionmatrix.bin');
% A = fread(fileID,'double');
% A = reshape(A, [3,numel(A)/3]);
% A(1,:) = A(1,:)+1;
% A(2,:) = A(2,:)+1;
% A = spconvert(A');
% fclose(fileID);
load('/home/mory/data/tests/rtksparseprojectionmatrixexport/full/matlabRewritten.mat', 'B');

% Set reconstruction parameters with Mechlem
ObjectSize = 380;
NbIters = 100;
runOnGPU = false;
NbSubsets = 4; % More than that speeds up computation, but makes it unstable, and ultimately, divergent
regulWeights = [10000, 10000, 10]; % Huber is scaled by delta, so where there is a low delta, there should be a high lambda to compensate
delta_huber = [0.001, 0.001, 0.01];

% Run Mechlem2017
[Mechlem2017_iterates, Mechlem2017_costs ]= Mechlem2017(measuredCounts, ObjectSize, pix_per_proj, B, M, S, T, NbIters, regulWeights, NbSubsets, delta_huber, runOnGPU);

% View results
PlayIteratesTransposed(Mechlem2017_iterates, [0 0 0], [0.0015 0.0015 1.5]);
result = Mechlem2017_iterates(:,:,NbIters);
result = reshape(result, [ObjectSize ObjectSize 3]);

writeImageResult(result(:,:,1)', ['/home/mory/Documents/Presentations/20180226_ConferenceCallWithPRH/figures/iodine_Mechlem2017.png'], 0, 0.0015);
writeImageResult(result(:,:,2)', ['/home/mory/Documents/Presentations/20180226_ConferenceCallWithPRH/figures/gadolinium_Mechlem2017.png'], 0, 0.0015);
writeImageResult(result(:,:,3)', ['/home/mory/Documents/Presentations/20180226_ConferenceCallWithPRH/figures/water_Mechlem2017.png'], 0, 1.5);

% Compare with PRH-reconstructions
load('/home/mory/data/Hamburg/2017_Jan_17.7905.1_materialDecomp/PhilipsReconstruction/2017_Jan17.7905.1_materialImages.mat', 'materialImages');
PRH_result = cat(3, materialImages.iodine(:,:,4), materialImages.gadolinium(:,:,4), materialImages.water(:,:,4));
PlayIteratesTransposed(reshape(PRH_result, [512*512 3 1]), [0 0 0], [0.0015 0.0015 1.5]);

writeImageResult(PRH_result(:,:,1)', ['/home/mory/Documents/Presentations/20180226_ConferenceCallWithPRH/figures/iodine_PRH.png'], 0, 0.0015);
writeImageResult(PRH_result(:,:,2)', ['/home/mory/Documents/Presentations/20180226_ConferenceCallWithPRH/figures/gadolinium_PRH.png'], 0, 0.0015);
writeImageResult(PRH_result(:,:,3)', ['/home/mory/Documents/Presentations/20180226_ConferenceCallWithPRH/figures/water_PRH.png'], 0, 1.5);

