function [spot_flux, correl_m] = spot(measured_condition)
% SPOT produces a condition-specific metabolic flux distribution from transcriptomic data, which can be used when objective function is unknown 

% by Min Kyung Kim Oct, 21, 2015

%load gene expression data
load test_expr;

%load a genome-scale metabolic model
load ijo1366_dc_test;
% load ijo1366_ac_test;

if (~all(fbamodel.vmin == 0 | fbamodel.vmin == -Inf) ...
    && ~all(fbamodel.vmax == 0 | fbamodel.vmax == Inf))
    disp('FBA model bounds must be 0 or +/- Infinity');
    return;
end

% process gene expression data based on GPR association
pts_expr = zeros(length(fbamodel.pts), length(expr_cols)); 

for i = 1:length(fbamodel.pts) 
    pts = fbamodel.pts{i};

   
    level = 1; 
    j = 1; 
    op = [];  
    val = []; 
    
    while j <= length(pts) 
        
        if pts(j) == '('  
            level = level + 1; 
            j = j + 1; 
        
        elseif pts(j) == ')'
            newval = val(level, :);
            level = level - 1;                       
           
            if  length(op) >= level && op(level) > 0 
            
                if  op(level) == 1  
                    val(level, :) = val(level, :) + newval; 
                elseif op(level) == 2 
                    val(level, :) = min([val(level, :); newval]);
                end 
                op(level) = 0; 
            else 
                val(level, :) = newval; 
            end 
            j = j + 1; 
        
        elseif strcmp(pts(j:(j + 3)), ' or ') 
       
            op(level) = 1; 
            j = j + 4;
        
        elseif strcmp(pts(j:(j + 4)), ' and ') 
            op(level) = 2; 
            j = j + 5; 
        
        else 
            gene = pts(j:(j + 4)); 
            iexpr = find(strcmp(gene, expr{1})); 
            if ~isempty(iexpr) 
                newval = mean(expr{2}(iexpr, :), 1); 
            else 
                newval = Inf(1, length(expr_cols));
            end     
     
            if  length(op) >= level && op(level) > 0          
                if  op(level) == 1 
                    val(level, :) = val(level, :) + newval;  
                elseif op(level) == 2 
                    val(level, :) = min([val(level, :); newval]); 
                end 
                op(level) = 0;
            else 
                val(level, :) = newval; 
            end 
            
            j = j + 5; 
       
        end 
      
    end
     pts_expr(i, :) = val(1, :);  
      
end 

all_expr = fbamodel.G' * pts_expr;
all_expr(all(fbamodel.G == 0)', :) = Inf;  

gene_data = all_expr(:,measured_condition); 

% max vg, s.t. Sv = 0 , 0/-Inf <= v <= 0/+Inf
brev = fbamodel.vmin < 0 & fbamodel.vmax > 0; 
birr = ~brev;
irev = find(brev);
iirr = find(birr);

fbamodel_rev = sparse(0, fbamodel.nrxn);
for i = 1:length(irev)
    fbamodel_rev(i, irev(i)) = 1;
end

fbamodel_irr = sparse(0, fbamodel.nrxn);
for i = 1:length(iirr)
    fbamodel_irr(i, iirr(i)) = 1;
end

gene_rev = fbamodel_rev*gene_data;
gene_rev(gene_rev == inf) = 0;
gene_irr = fbamodel_irr*gene_data;
gene_irr(gene_irr == inf) = 0;

S_rev = fbamodel.S*fbamodel_rev';
S_irr = fbamodel.S*fbamodel_irr';

vmin_rev =  fbamodel_rev*fbamodel.vmin;
vmax_rev =  fbamodel_rev*fbamodel.vmax;
vmin_irr =  fbamodel_irr*fbamodel.vmin;
vmax_irr =  fbamodel_irr*fbamodel.vmax;

model1.A = [S_irr S_rev -S_rev] ; 
model1.sense = char('='* ones(size(model1.A,1), 1)); 
model1.rhs = [zeros(fbamodel.nmetab,1)] ; 
model1.modelsense = 'max'; 
model1.lb = [vmin_irr; zeros(length(vmin_rev),1); zeros(length(vmin_rev),1)];  
model1.ub = [vmax_irr; vmax_rev; -vmin_rev]; 
model1.obj = [gene_irr; gene_rev ; gene_rev]; 

% % |v|^2 <= 1 => quadratic constraint
model1.quadcon(1).Qc = speye(length(model1.lb)); 
model1.quadcon(1).q =  zeros(length(model1.lb), 1);  
model1.quadcon(1).rhs = 1; 

params1.outputflag = 0; %in order to shut off Gurobi output
params1.BarHomogeneous = 1; %setting homogeneous barrier algorithm to 1 forces it on.
result1 = gurobi(model1, params1);

maxvg_v = result1.x; % resulting flux distribution (decomposed network)

v_irr = maxvg_v(1:length(vmin_irr));
f_v_rev = maxvg_v(length(vmin_irr)+1: length(vmin_irr) + length(vmin_rev));
b_v_rev = maxvg_v(length(vmin_irr) + length(vmin_rev) + 1: length(maxvg_v));
v_rev = f_v_rev - b_v_rev;

spot_flux = fbamodel_irr'*v_irr + fbamodel_rev'*v_rev;

%%
% import txt file shows the predicted flux - the measured flux association 
FLUX_LIST = 'fluxmatched.txt'; %%%%%%%%%%%

fid = fopen(FLUX_LIST); 
match_list = textscan(fid, '%s');
fclose(fid);

% process the predicted fluxes according to the meausred flux - the predicted flux association
modified_correl = zeros(size(match_list{1}));

for i = 1:length(match_list{1})
    silico_vitro = match_list{1}{i};
    
    lv = 1; 
    j = 1; 
    op = [];  
    val = []; 
     
    while j <= length(silico_vitro)
        
        if silico_vitro(j) == '(' 
            lv = lv + 1; 
            j = j + 1; 
        
        elseif silico_vitro(j) == ')' 
            newval = val(lv);          
            lv = lv - 1;                         
           
            if  length(op) >= lv && op(lv) > 0 
            
                if  op(lv) == 1 
                    val(lv) = val(lv) + newval;
                    
                elseif op(lv) == 2 
                    val(lv) = min([val(lv); newval]);

                end 
                op(lv) = 0; 
            else 
                val(lv) = newval; 

            end 
            j = j + 1; 
        
        elseif strcmp(silico_vitro(j:(j + 3)), '_OR_') 
       
            op(lv) = 1; 
            j = j + 4; 
        
        elseif strcmp(silico_vitro(j:(j + 4)), '_AND_') 
            op(lv) = 2; 
            j = j + 5; 
        
        else 
            % section the names of fluxes
            k = strfind(silico_vitro(j:end),')');
            if ~isempty(k)
                flux_string = silico_vitro(j:(j + k(1) - 2));
            else 
                flux_string = silico_vitro(j:end);
            end
            
            % find and assign eflux value of corresponding flux name
            if strncmp('-', flux_string, 1) == 1 
               str_to_compare = ['R_' flux_string(2:end)];
               iflux = find(strcmp(str_to_compare, fbamodel.rxns)); 
               if ~isempty(iflux) 
                   newval = -spot_flux(iflux); 
                                      
               else 
                   newval = inf; 
                   
               end
            else
               iflux = find(strcmp(['R_' flux_string], fbamodel.rxns));
               if ~isempty(iflux) 
                   newval = spot_flux(iflux); 
                  
               else 
                   newval = inf;
                 
               end                
            end
            
            if  length(op) >= lv && op(lv) > 0
                if  op(lv) == 1 
                    val(lv) = val(lv) + newval;
                    
                elseif op(lv) == 2 
                    val(lv) = min([val(lv); newval]); 

                end 
                op(lv) = 0; 
                
            else 
                val(lv) = newval; 
            end 
            
            j = j + length(flux_string);       
                           
        end 
      
    end
    modified_correl(i) = val(1);  

end
                
    
%%
insilico_flux = modified_correl;

% import measured fluxes 
condition = {'_rf', '_pgm', '_pgi', '_gapC', '_zwf', '_rpe', '_wt5', '_wt7'};
MEASURED_FLUX = ['measured_flux' condition{measured_condition} '.txt'];
fid = fopen(MEASURED_FLUX); 
measured_flux = fscanf(fid, '%f');
fclose(fid);

% calculate correl b/w measured flux & predicted flux - except for NaN elements in the measured fluxes
w = measured_flux(~isnan(measured_flux))/100;
v = insilico_flux(~isnan(measured_flux));
correl_m = dot(v,w)/(norm(v)*norm(w));