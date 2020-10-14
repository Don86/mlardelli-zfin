%fbamodel = sbml_to_mat('iJO1366_DC_glucose.xml');

fbamod = readCbModel('zebragem_2_0_updated_20200228.xml');

% Annoyingly, this runs silently.
% Output are new vars created in the workspace: 
% test_expr.mat expr expr_cols
% Modified to take in new data
expr_txt_to_mat;

% run SPOT
%uses gene expression data and a genome-scale metabolic model
[spot_flux, correl_m] = spot(1);