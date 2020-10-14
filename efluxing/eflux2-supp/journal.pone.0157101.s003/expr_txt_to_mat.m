% This simiple matlab script reads a gene expression data (.txt) and returns MAT-file

% by Min Kyung Kim, Oct, 21, 2015

EXPR_DATA_FNAME = 'test_expr.txt'; % txt file name where gene expresison data is stored

fid = fopen(EXPR_DATA_FNAME);
C = textscan(fid, '%s%f%f%f%f%f%f%f%f', 'Headerlines', 1, 'Delimiter', '\t'); % first column : gene names used in the model, the other columns : measured gene expression data
fclose(fid);
gene_names = C{1};
expr_data = cell2mat(C(2:end));

expr = cell(1, 2);
expr{1} = gene_names;
expr{2} = expr_data;
expr_cols = { 'RF', 'pgm', 'pgi', 'gapC', 'zwf', 'rpe', 'WT_0_5', 'WT_0_7' }; % experimental conditions

save test_expr.mat expr expr_cols;
