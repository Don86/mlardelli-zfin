addpath code;
initCobraToolbox;

model = readCbModel('models/Ec_core_flux1','fileType','SBML');

% run baseline FBA
fba = optimizeCbModel(model);

% sample, uses by default 10 chains of 200 samples with thinning = 50
sol = bfba(model, fba);

%% visualise

% plot first 15 flux distributions with fba
figure; plotfluxes(model, sol, 1:15, fba);

% plot 7x7 flux pair grid
figure; plotfluxpair(model, sol, [10 2 40 42 45 53 55], fba);

% plot 8 example flux pairs
figure; plotfluxpair2(model, sol, [30 40; 30 41; 40 47; 15 17; 4 5; 1 6; 10 18; 42 45], fba);
