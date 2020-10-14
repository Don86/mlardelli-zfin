%    FBAMODEL = READ_SBML(SBMLFNAME) reads the SBML file SBMLFNAME.txt
%    and returns the FBA model contained in it as FBAMODEL

% by Desmond S. Lun on Oct, 9, 2012
function fbamodel = sbml_to_mat(sbmlfname) 

model = TranslateSBML(sbmlfname, 1); 

fbamodel.nrxn = length(model.reaction); 
fbamodel.nmetab = length(model.species); 
fbamodel.S = sparse(fbamodel.nmetab, fbamodel.nrxn); 
fbamodel.rxns = cell(fbamodel.nrxn, 1); 
fbamodel.metabs = cell(fbamodel.nmetab, 1); 
fbamodel.f = zeros(fbamodel.nrxn, 1); 
fbamodel.g = zeros(fbamodel.nrxn, 1); 
fbamodel.vmin = zeros(fbamodel.nrxn, 1); 
fbamodel.vmax = zeros(fbamodel.nrxn, 1); 
fbamodel.present = true(fbamodel.nrxn, 1); 

fbamodel.pts = {};
fbamodel.G = sparse(0, fbamodel.nrxn);
ipt = 1;

boundary = zeros(fbamodel.nmetab, 1);

for imetab = 1:fbamodel.nmetab 
    fbamodel.metabs{imetab} = model.species(imetab).id; 
    fbamodel.metabname{imetab} = model.species(imetab).name; 
    
    boundary(imetab) = model.species(imetab).boundaryCondition; 
end


for irxn = 1:fbamodel.nrxn 
    fbamodel.rxns{irxn} = model.reaction(irxn).id; 
    fbamodel.rxnname{irxn} = model.reaction(irxn).name; 
    
    for ireactant = 1:length(model.reaction(irxn).reactant) 
        fbamodel.S(strcmp(model.reaction(irxn).reactant(ireactant).species, fbamodel.metabs), irxn) ...
            = -model.reaction(irxn).reactant(ireactant).stoichiometry; 
    end
    
    for iproduct = 1:length(model.reaction(irxn).product)
        fbamodel.S(strcmp(model.reaction(irxn).product(iproduct).species, fbamodel.metabs), irxn) ...
            = model.reaction(irxn).product(iproduct).stoichiometry; 
    end
    
    for iparam = 1:length(model.reaction(irxn).kineticLaw.parameter)
        param = model.reaction(irxn).kineticLaw.parameter(iparam);
        if strcmp('LOWER_BOUND', param.id)
            fbamodel.vmin(irxn) = param.value;
        elseif strcmp('UPPER_BOUND', param.id)
            fbamodel.vmax(irxn) = param.value;
        elseif strcmp('OBJECTIVE_COEFFICIENT', param.id)
            fbamodel.f(irxn) = param.value; 
        end
    end
    
    % by Min Kyung Kim on May, 30, 2013
    for inotes = 1:length(model.reaction(irxn).notes) 
        istart = strfind(model.reaction(irxn).notes, 'GENE_ASSOCIATION:') + 18; % site where gene association starts, after':'
        iend = strfind(model.reaction(irxn).notes, '</p>') -1;  % site where gene association ends , before '<'
        gene_association{irxn} = model.reaction(irxn).notes(istart:iend);% sectioning string from istart to iend
        if ~isempty(gene_association{irxn})
            pt = gene_association{irxn};
            tf = strcmp(pt,fbamodel.pts);
            if any(tf)
                fbamodel.G(tf, irxn) = 1;
            else
                fbamodel.G(ipt, irxn) = 1;
                fbamodel.pts{ipt} = pt;
                ipt = ipt + 1;
            end
        end
    end
           
end

fbamodel.metabs = fbamodel.metabs(~boundary); 
fbamodel.S = fbamodel.S(~boundary, :);
fbamodel.nmetab = length(fbamodel.metabs);
fbamodel.metabname = fbamodel.metabname(~boundary);
