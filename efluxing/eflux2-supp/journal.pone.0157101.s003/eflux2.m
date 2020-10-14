function [eflux2_flux, correl_m] = eflux2(measured_condition)
% E-Flux2 produces a condition-specific metabolic flux distribution from transcriptomic data, which can be used when objective function is known 

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

for i = 1:length(fbamodel.pts) % length(fbamodel.pts): the number of GPR associations in the model
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
            
                if  op(level) == 1  % if op(level) = OR
                    val(level, :) = val(level, :) + newval; % OR relationship among genes => sum
                elseif op(level) == 2 %if op(level) = AND
                    val(level, :) = min([val(level, :); newval]); % AND relationship among genes => minimum
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
        
        else % other than '(', ')', '_OR_', '_AND_' -> that is, gene name (b number)
       
            gene = pts(j:(j + 4)); 
            
%             % if gene symbols have varied lengths
%             k = strfind(pts(j:end),')');
%             if ~isempty(k)
%                 gene_string = pts(j:(j + k(1) - 2));
%             else
%                 gene_string = pts(j:end);
%             end
            
            iexpr = find(strcmp(gene, expr{1})); 
            if ~isempty(iexpr)
                newval = mean(expr{2}(iexpr, :), 1); 
            else % if there's no corresponding measured gene exprssion data for the variable, "gene"
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
all_expr(all(fbamodel.G == 0)', :) = Inf; % for reactions without corresponding GPR association
gene_data = all_expr(:,measured_condition); 

%Step 1 of E-Flux2: peform E-Flux
%change the lower and upper bounds of each reaction using the gene expression data 
fbamodel.vmax = gene_data ;  
fbamodel.vmin(fbamodel.vmin ~= 0) = -gene_data(fbamodel.vmin ~= 0, 1);
eflux = flux_balance(fbamodel, true);

%Step 2 of E-Flux2: minimize l2 norm
model2.A = [fbamodel.S;fbamodel.f'] ; 
model2.sense = char('='* ones(size(model2.A,1), 1)); 
model2.rhs = [zeros(fbamodel.nmetab,1); fbamodel.f'*eflux] ; % maintain maximized biomass flux value calculated in step 1
model2.modelsense = 'min'; % max = -1, min = 1 
model2.lb = fbamodel.vmin;
model2.ub = fbamodel.vmax;
model2.Q = sparse(diag(ones(fbamodel.nrxn,1)));
model2.obj = zeros(fbamodel.nrxn,1);

params2.outputflag = 0; %in order to shut off Gurobi output
params2.BarHomogeneous = 1; %setting homogeneous barrier algorithm to 1 forces it on.

result2 = gurobi(model2, params2);
eflux2_flux = result2.x;

%%
% import txt file shows the predicted flux - the measured flux association 
FLUX_LIST = 'fluxmatched.txt'; 
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
           
            k = strfind(silico_vitro(j:end),')');
            if ~isempty(k)
                flux_string = silico_vitro(j:(j + k(1) - 2));
            else 
                flux_string = silico_vitro(j:end);
            end
            
            % find and assign eflux value of corresponding flux name
            if strncmp('-', flux_string, 1) == 1 % flux name has a negative value in front of itself
               str_to_compare = ['R_' flux_string(2:end)];
               iflux = find(strcmp(str_to_compare, fbamodel.rxns)); 
               if ~isempty(iflux) % matched with the flux name
                   newval = -eflux2_flux(iflux); % assign 'negative' value of eflux
                                      
               else
                   newval = inf; 
                   
               end
            else
               iflux = find(strcmp(['R_' flux_string], fbamodel.rxns));
               if ~isempty(iflux) 
                   newval = eflux2_flux(iflux); 
                  
               else 
                   newval = inf;
                 
               end                
            end
            
            if  length(op) >= lv && op(lv) > 0
                if  op(lv) == 1 
                    val(lv) = val(lv) + newval; % OR relationship among fluxes => sum
                    
                elseif op(lv) == 2 
                    val(lv) = min([val(lv); newval]); % AND relationhsip among fluxes => minimum

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
m = measured_flux(~isnan(measured_flux))/100;
p = insilico_flux(~isnan(measured_flux));
correl_m = dot(p,m)/(norm(p)*norm(m));

end

%flux balance analysis
% by Desmond S. Lun Feb, 26, 2013

function [v, fmax, fmin] = flux_balance(fbamodel, quiet)

if nargin < 2
    quiet = false;
end

nrxn   = fbamodel.nrxn;
nmetab = fbamodel.nmetab;

c = fbamodel.f;
lb = fbamodel.vmin;
ub = fbamodel.vmax;
   
Aeq = [ fbamodel.S                                                                                                               
        sparse(1:nnz(~fbamodel.present), find(~fbamodel.present), ones(nnz(~fbamodel.present), 1), nnz(~fbamodel.present), nrxn) ];
beq = [ zeros(nmetab, 1);
        zeros(nnz(~fbamodel.present), 1); ];

vartype = char('C' * ones(nrxn, 1));
            
options = optimset('Display', 'off');
            
[x, vbiomass, exitflag] = solve_milp(c, [], [], Aeq, beq, lb, ub, vartype, -1, options);

fmin = NaN;
fmax = NaN;
max_vsynth = NaN;
v = [];

if exitflag > 0
    v = x(1:nrxn);

    if nnz(fbamodel.g) > 0
        c = fbamodel.g;

        lb_nobiomass = lb;
        ub_nobiomass = ub;
        lb_nobiomass(fbamodel.f > 0) = 0;
        ub_nobiomass(fbamodel.f < 0) = 0;
        [x, max_vsynth] = solve_milp(c, [], [], Aeq, beq, lb_nobiomass, ub_nobiomass, vartype, -1, options);         

        Aeq = [ Aeq;
                fbamodel.f' ];
        beq = [ beq;
                vbiomass ];

        [x, fmin] = solve_milp(c, [], [], Aeq, beq, lb, ub, vartype, 1, options);
        [x, fmax, exitflag] = solve_milp(c, [], [], Aeq, beq, lb, ub, vartype, -1, options);
        if exitflag > 0
            v = x(1:nrxn);
        end
    end
end

if ~quiet
    fprintf('Biomass flux:    %f\n', vbiomass);
    if ~isnan(fmin)
        fprintf('Synthetic flux:  [%f, %f] of %f\n', fmin, fmax, max_vsynth);
    end
end

end

function [x, fval, exitflag] = solve_milp(f, A, b, Aeq, beq, lb, ub, vartype, sense, options)

SOLVER = 'Gurobi';

if nargin < 10
    options = optimset('Display', 'final');
end
if nargin < 9
    sense = 1;
end
if nargin < 8
    vartype = char('C' * ones(size(f)));
end
if nargin < 6
    lb = [];
    ub = [];
end
if nargin < 4
    Aeq = [];
    beq = [];
end

switch SOLVER
   
    case 'Gurobi'
        params = struct();
        params.FeasibilityTol = 1e-9;
        params.IntFeasTol = 1e-9;
        switch options.Display
            case 'off'
                params.OutputFlag = 0;
            case 'iter'
                params.DisplayInterval = 5;
            case 'final'
                params.OutputFlag = 0;
            case 'notify'
                params.OutputFlag = 0;
        end
        
        model = struct();
        model.A = [A; Aeq];
        model.obj = f;
        model.sense = [ char('<' * ones(size(A, 1), 1)); char ('=' * ones(size(Aeq, 1), 1)) ];
        model.rhs = [b; beq];
        model.lb = lb;
        model.ub = ub;     
        model.vtype = vartype;
        if sense < 0
            model.modelsense = 'max';
        else
            model.modelsense = 'min';
        end
        
        result = gurobi(model, params);
        
        if strcmp(result.status, 'OPTIMAL')
            exitflag = 1;
            x = result.x;
            fval = result.objval;
        else
            exitflag = -1;
            x = [];
            fval = NaN;
        end
end

end

