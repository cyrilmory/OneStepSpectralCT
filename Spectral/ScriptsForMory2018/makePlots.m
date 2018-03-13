% Load all results
load('/home/mory/data/MatlabOneStep/Cai/Cai2013_fessler_100000_100000_100.mat', 'Cai2013_iterates', 'Cai2013_costs');
load('/home/mory/data/MatlabOneStep/Long/Long2014_100000_100000_100.mat', 'Long2014_iterates', 'Long2014_costs');
load('/home/mory/data/MatlabOneStep/Weidinger/Weidinger2016_10000_10000_10.mat', 'Weidinger2016_iterates', 'Weidinger2016_costs');
load('/home/mory/data/MatlabOneStep/Barber/Barber2016_100_100_10000.mat', 'Barber2016_iterates', 'Barber2016_gaps', 'Barber2016_costs');
load('/home/mory/data/MatlabOneStep/Mechlem/Mechlem2017_100000_100000_100.mat', 'Mechlem2017_iterates', 'Mechlem2017_costs');

% Write the final iterates for each method
results=cell(5,2);
results{1,1}=Cai2013_iterates;          results{1,2}='Cai2013';
results{2,1}=Long2014_iterates;         results{2,2}='Long2014';
results{3,1}=Weidinger2016_iterates;    results{3,2}='Weidinger2016';
results{4,1}=Barber2016_iterates;       results{4,2}='Barber2016';
results{5,1}=Mechlem2017_iterates;      results{5,2}='Mechlem2017';

% Generate ground truth
ObjectSize = sqrt(size(results{1,1}, 1));
unit = ObjectSize/8;
water=zeros(ObjectSize);
water(unit+1:7*unit, unit+1:7*unit) = ones(6*unit, 6*unit);
iodine=zeros(ObjectSize);
iodine(2*unit+1:3*unit, 2*unit+1:3*unit) = ones(unit, unit) * 0.01;
gadolinium=zeros(ObjectSize);
gadolinium(4*unit+1:5*unit, 5*unit+1:6*unit) = ones(unit, unit) * 0.01;

% Initialize the Relative MSE table
RMSE_table = zeros(5,3);

for m=1:5
    result = results{m,1};
    result = result(:,:,end);
    
    % Compute MSE with the ground truth
    RMSE_table(m, 1) = sum((result(:,1) - iodine(:)).^2)     / sum(iodine(:).^2);
    RMSE_table(m, 2) = sum((result(:,2) - gadolinium(:)).^2) / sum(gadolinium(:).^2);
    RMSE_table(m, 3) = sum((result(:,3) - water(:)).^2)      / sum(water(:).^2);
    
    % Write one image per material
    result = reshape(result, [ObjectSize ObjectSize 3]);
%     writeImageResult(result(:,:,1), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/iodine_',     results{m,2}, '.png'], 0.0075, 0.0125);
%     writeImageResult(result(:,:,2), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/gadolinium_', results{m,2}, '.png'], 0.0075, 0.0125);
%     writeImageResult(result(:,:,3), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/water_',      results{m,2}, '.png'], 0.75, 1.25);
    writeImageResult(result(:,:,1), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/iodine_',     results{m,2}, '.png'], 0, 0.015);
    writeImageResult(result(:,:,2), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/gadolinium_', results{m,2}, '.png'], 0, 0.015);
    writeImageResult(result(:,:,3), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/water_',      results{m,2}, '.png'], 0, 1.5);
end

% To determine the concentration in ROIs, crop by a few pixels
% in order to avoid partial volume effect on the borders
cropROI=2;

%% Get iodine concentration curve for each, and plot it
figure(1);
hold off

[y_Cai,x_Cai,~] = getMaterialConcentration('iodine', Cai2013_iterates, 10, 0, cropROI);
semilogx(x_Cai,y_Cai, '+','linewidth', 2)
hold on

[y_Long,x_Long,~] = getMaterialConcentration('iodine', Long2014_iterates, 10, 0, cropROI);
semilogx(x_Long,y_Long, '--', 'linewidth', 2)

[y_Weidinger,x_Weidinger,~] = getMaterialConcentration('iodine', Weidinger2016_iterates, 10, 0, cropROI);
semilogx(x_Weidinger,y_Weidinger, '-.', 'linewidth', 2)

[y_Barber,x_Barber,~] = getMaterialConcentration('iodine', Barber2016_iterates, 10, 0, cropROI);
semilogx(x_Barber,y_Barber, 'o', 'linewidth', 2)

[y_Mechlem,x_Mechlem,~] = getMaterialConcentration('iodine', Mechlem2017_iterates, 1, 0, cropROI);
semilogx(x_Mechlem,y_Mechlem, 'x', 'linewidth', 2)

target_x = 1:max([x_Cai(end), x_Long(end), x_Weidinger(end), x_Barber(end), x_Mechlem(end)]);
target_y = ones(size(target_x))*0.01;
semilogx(target_x, target_y,'linewidth', 2)

lgd = legend('Cai2013','Long2014','Weidinger2016','Barber2016','Mechlem2017', 'Ground truth', 'Location', 'southeast');
lgd.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Iodine concentration in ROI','fontsize',14)
set(gcf, 'Position', [500, 500, 800, 600])
hold off

% Rescale the graph manually before running this line, otherwise the legend
% appears on top of some curves
saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/iodineConcentrations.png');

%% Get gadolinium concentration curve for each, and plot it
figure(2);
hold off

[y_Cai,x_Cai,~] = getMaterialConcentration('gadolinium', Cai2013_iterates, 10, 0, cropROI);
semilogx(x_Cai,y_Cai, '+','linewidth', 2)
hold on

[y_Long,x_Long,~] = getMaterialConcentration('gadolinium', Long2014_iterates, 10, 0, cropROI);
semilogx(x_Long,y_Long, '--', 'linewidth', 2)

[y_Weidinger,x_Weidinger,~] = getMaterialConcentration('gadolinium', Weidinger2016_iterates, 10, 0, cropROI);
semilogx(x_Weidinger,y_Weidinger, '-.', 'linewidth', 2)

[y_Barber,x_Barber,~] = getMaterialConcentration('gadolinium', Barber2016_iterates, 10, 0, cropROI);
semilogx(x_Barber,y_Barber, 'o', 'linewidth', 2)

[y_Mechlem,x_Mechlem,~] = getMaterialConcentration('gadolinium', Mechlem2017_iterates, 1, 0, cropROI);
semilogx(x_Mechlem,y_Mechlem, 'x', 'linewidth', 2)

target_x = 1:max([x_Cai(end), x_Long(end), x_Weidinger(end), x_Barber(end), x_Mechlem(end)]);
target_y = ones(size(target_x))*0.01;
semilogx(target_x, target_y,'linewidth', 2)

lgd = legend('Cai2013','Long2014','Weidinger2016','Barber2016','Mechlem2017', 'Ground truth', 'Location', 'southeast');
lgd.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Gadolinium concentration in ROI','fontsize',14)
set(gcf, 'Position', [500, 500, 800, 600])
hold off

% Rescale the graph manually before running this line, otherwise the legend
% appears on top of some curves
saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/gadoliniumConcentrations.png');

%% Plot the costs functions
Cai2013_costs_toplot = Cai2013_costs(10:10:end);
Cai2013_costs_toplot = Cai2013_costs_toplot - min(Cai2013_costs_toplot);

Long2014_costs_toplot = Long2014_costs(10:10:end);
Long2014_costs_toplot = Long2014_costs_toplot - min(Long2014_costs_toplot);

Weidinger2016_costs_toplot = Weidinger2016_costs(10:10:end);
Weidinger2016_costs_toplot = Weidinger2016_costs_toplot - min(Weidinger2016_costs_toplot);

Barber2016_costs_toplot = Barber2016_costs;
Barber2016_costs_toplot = Barber2016_costs_toplot - min(Barber2016_costs_toplot);

Mechlem2017_costs_toplot = Mechlem2017_costs;
Mechlem2017_costs_toplot = Mechlem2017_costs_toplot - min(Mechlem2017_costs_toplot);

figure(3);
hold off
loglog(x_Cai,Cai2013_costs_toplot, '+','linewidth', 2)
hold on 
loglog(x_Long, Long2014_costs_toplot, '--','linewidth', 2)
loglog(x_Weidinger, Weidinger2016_costs_toplot, '-.','linewidth', 2)
load('/home/mory/data/MatlabOneStep/Barber/preps.mat', 'y', 'A', 'M', 'S', 'T');
% Barber2016_costs = zeros(size(Barber2016_iterates, 3), 1);
% for it=1:size(Barber2016_iterates, 3)
%     Barber2016_costs(it) = BarberComputeCosts(Barber2016_iterates(:,:,it), y, A, M, S);
% end
loglog(1:10000, Barber2016_costs_toplot, 'o','linewidth', 2)
loglog(x_Mechlem, Mechlem2017_costs_toplot, 'x','linewidth', 2)

lgd2 = legend('Cai2013','Long2014','Weidinger2016','Barber2016','Mechlem2017', 'Location', 'southwest');
lgd2.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Costs - min(Costs)','fontsize',14)
set(gcf, 'Position', [500, 500, 600, 450])
hold off

% Rescale the graph manually before running this line, otherwise the legend
% appears on top of some curves
saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/CostFunctions.png');

%% Thresholds in the concentration of iodine
TargetConcentration = 0.01;
[~,~,Cai_50_iodine] = getMaterialConcentration('iodine', Cai2013_iterates, 10, TargetConcentration*0.5, cropROI);
[~,~,Long_50_iodine] = getMaterialConcentration('iodine', Long2014_iterates, 10, TargetConcentration*0.5, cropROI);
[~,~,Weidinger_50_iodine] = getMaterialConcentration('iodine',Weidinger2016_iterates, 10, TargetConcentration*0.5, cropROI);
[~,~,Barber_50_iodine] = getMaterialConcentration('iodine', Barber2016_iterates, 10, TargetConcentration*0.5, cropROI);
[~,~,Mechlem_50_iodine] = getMaterialConcentration('iodine', Mechlem2017_iterates, 1, TargetConcentration*0.5, cropROI);

[~,~,Cai_80_iodine] = getMaterialConcentration('iodine', Cai2013_iterates, 10, TargetConcentration*0.8, cropROI);
[~,~,Long_80_iodine] = getMaterialConcentration('iodine', Long2014_iterates, 10, TargetConcentration*0.8, cropROI);
[~,~,Weidinger_80_iodine] = getMaterialConcentration('iodine', Weidinger2016_iterates, 10, TargetConcentration*0.8, cropROI);
[~,~,Barber_80_iodine] = getMaterialConcentration('iodine', Barber2016_iterates, 10, TargetConcentration*0.8, cropROI);
[~,~,Mechlem_80_iodine] = getMaterialConcentration('iodine', Mechlem2017_iterates, 1, TargetConcentration*0.8, cropROI);

[~,~,Cai_50_gadolinium] = getMaterialConcentration('gadolinium', Cai2013_iterates, 10, TargetConcentration*0.5, cropROI);
[~,~,Long_50_gadolinium] = getMaterialConcentration('gadolinium', Long2014_iterates, 10, TargetConcentration*0.5, cropROI);
[~,~,Weidinger_50_gadolinium] = getMaterialConcentration('gadolinium', Weidinger2016_iterates, 10, TargetConcentration*0.5, cropROI);
[~,~,Barber_50_gadolinium] = getMaterialConcentration('gadolinium', Barber2016_iterates, 10, TargetConcentration*0.5, cropROI);
[~,~,Mechlem_50_gadolinium] = getMaterialConcentration('gadolinium', Mechlem2017_iterates, 1, TargetConcentration*0.5, cropROI);

[~,~,Cai_80_gadolinium] = getMaterialConcentration('gadolinium', Cai2013_iterates, 10, TargetConcentration*0.8, cropROI);
[~,~,Long_80_gadolinium] = getMaterialConcentration('gadolinium', Long2014_iterates, 10, TargetConcentration*0.8, cropROI);
[~,~,Weidinger_80_gadolinium] = getMaterialConcentration('gadolinium', Weidinger2016_iterates, 10, TargetConcentration*0.8, cropROI);
[~,~,Barber_80_gadolinium] = getMaterialConcentration('gadolinium', Barber2016_iterates, 10, TargetConcentration*0.8, cropROI);
[~,~,Mechlem_80_gadolinium] = getMaterialConcentration('gadolinium', Mechlem2017_iterates, 1, TargetConcentration*0.8, cropROI);

thresholdsTable = [Cai_50_iodine,       Cai_80_iodine,          Cai_50_gadolinium,       Cai_80_gadolinium; ...
                   Long_50_iodine,      Long_80_iodine,         Long_50_gadolinium,      Long_80_gadolinium; ... 
                   Weidinger_50_iodine, Weidinger_80_iodine,    Weidinger_50_gadolinium, Weidinger_80_gadolinium; ...
                   Barber_50_iodine,    Barber_80_iodine,       Barber_50_gadolinium,    Barber_80_gadolinium; ...
                   Mechlem_50_iodine,   Mechlem_80_iodine,      Mechlem_50_gadolinium,   Mechlem_80_gadolinium];

timePerIter=[5.2; 265/200; 95/200; 204/200; 124/200];

% Where the threshold hasn't been reached, put value 10000 instead of 0
thresholdsTable(thresholdsTable==0)=10000;

% Get the times to reach the thresholds
timeTable=thresholdsTable .* timePerIter;

%% Compare Cai2013 with different mu-preconditioning

Cai2013_none = load('/home/mory/data/MatlabOneStep/Cai/Cai2013_none_100000_100000_100.mat');
Cai2013_none_iterates = Cai2013_none.Cai2013_iterates;
Cai2013_none_costs = Cai2013_none.Cai2013_costs;
Cai2013_normalize = load('/home/mory/data/MatlabOneStep/Cai/Cai2013_normalize_100000_100000_100.mat');
Cai2013_normalize_iterates = Cai2013_normalize.Cai2013_iterates;
Cai2013_normalize_costs = Cai2013_normalize.Cai2013_costs;
Cai2013_orthonormalize = load('/home/mory/data/MatlabOneStep/Cai/Cai2013_orthonormalize_100000_100000_100.mat');
Cai2013_orthonormalize_iterates = Cai2013_orthonormalize.Cai2013_iterates;
Cai2013_orthonormalize_costs = Cai2013_orthonormalize.Cai2013_costs;
Cai2013_fessler = load('/home/mory/data/MatlabOneStep/Cai/Cai2013_fessler_100000_100000_100.mat');
Cai2013_fessler_iterates = Cai2013_fessler.Cai2013_iterates;
Cai2013_fessler_costs = Cai2013_fessler.Cai2013_costs;

Cai2013_none_costs_toplot = Cai2013_none_costs(10:10:end);
Cai2013_none_costs_toplot = Cai2013_none_costs_toplot - min(Cai2013_none_costs_toplot);

Cai2013_normalize_costs_toplot = Cai2013_normalize_costs(10:10:end);
Cai2013_normalize_costs_toplot = Cai2013_normalize_costs_toplot - min(Cai2013_normalize_costs_toplot);

Cai2013_orthonormalize_costs_toplot = Cai2013_orthonormalize_costs(10:10:end);
Cai2013_orthonormalize_costs_toplot = Cai2013_orthonormalize_costs_toplot - min(Cai2013_orthonormalize_costs_toplot);

Cai2013_fessler_costs_toplot = Cai2013_fessler_costs(10:10:end);
Cai2013_fessler_costs_toplot = Cai2013_fessler_costs_toplot - min(Cai2013_fessler_costs_toplot);

x_Cai = (1:size(Cai2013_none_iterates, 3)) * 10;

figure(4);
hold off
loglog(x_Cai,Cai2013_none_costs_toplot, '+','linewidth', 2)
hold on 
loglog(x_Cai, Cai2013_normalize_costs_toplot, '--','linewidth', 2)
loglog(x_Cai, Cai2013_orthonormalize_costs_toplot, '-.','linewidth', 2)
loglog(x_Cai, Cai2013_fessler_costs_toplot, 'x','linewidth', 2)

lgd3 = legend('none','normalize','orthonormalize','Fessler', 'Location', 'southwest');
lgd3.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Costs - min(Costs)','fontsize',14)
% set(gcf, 'Position', [500, 500, 600, 450])
hold off

saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/mupreconditionings_costs.png');

%% Get iodine concentration curve for each method of mu-preconditioning for Cai2013
figure(5);
hold off
cropROI = 2; 

[y_none,x_Cai,~] = getMaterialConcentration('gadolinium', Cai2013_none_iterates, 10, 0, cropROI);
semilogx(x_Cai,y_none, '+','linewidth', 2)
hold on

[y_normalize,x_Cai,~] = getMaterialConcentration('gadolinium', Cai2013_normalize_iterates, 10, 0, cropROI);
semilogx(x_Cai,y_normalize, '--', 'linewidth', 2)

[y_orthonormalize,x_Cai,~] = getMaterialConcentration('gadolinium', Cai2013_orthonormalize_iterates, 10, 0, cropROI);
semilogx(x_Cai,y_orthonormalize, '-.', 'linewidth', 2)

[y_fessler,x_Cai,~] = getMaterialConcentration('gadolinium', Cai2013_fessler_iterates, 10, 0, cropROI);
semilogx(x_Cai,y_fessler, 'o', 'linewidth', 2)

target_y = ones(size(x_Cai))*0.01;
semilogx(x_Cai, target_y,'linewidth', 2)

lgd4 = legend('none','normalize','orthonormalize','Fessler', 'Location', 'northwest');
lgd4.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Iodine concentration in ROI','fontsize',14)
set(gcf, 'Position', [500, 500, 600, 500])
hold off

saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/mupreconditionings_gadolinium.png');

