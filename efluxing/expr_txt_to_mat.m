% This simiple matlab script reads a gene expression data (.txt) and returns MAT-file

% by Min Kyung Kim, Oct, 21, 2015

% modified for mlardelli

EXPR_DATA_FNAME = 'expr_ncbi_ids_n21927.tsv'; % txt file name where gene expresison data is stored

fid = fopen(EXPR_DATA_FNAME);
% first column : gene names used in the model, the other columns : measured gene expression data
C = textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Headerlines', 1, 'Delimiter', '\t');
fclose(fid);
gene_names = C{1};
expr_data = cell2mat(C(2:end));

expr = cell(1, 2);
expr{1} = gene_names;
expr{2} = expr_data;
mlardelli_expr_cols = { 'wt_6_0_1', 'wt_6_0_2', 'wt_6_0_3', 'wt_6_0_4', 'wt_6_1_2',...
       'wt_6_1_3', 'wt_6_1_4', 'wt_6_1_1', 'q96_6_0_1', 'q96_6_0_2',...
       'q96_6_0_3', 'q96_6_0_4', 'q96_6_1_1', 'q96_6_1_2', 'q96_6_1_3',...
       'q96_6_1_4', 'wt_24_0_1', 'wt_24_0_2', 'wt_24_0_3', 'wt_24_0_4',...
       'wt_24_1_1', 'wt_24_1_2', 'wt_24_1_3', 'wt_24_1_4', 'q96_24_0_1',...
       'q96_24_0_2', 'q96_24_0_3', 'q96_24_0_4', 'q96_24_1_1', 'q96_24_1_2',...
       'q96_24_1_3', 'q96_24_1_4' }; % experimental conditions

save mlardelli_expr.mat expr mlardelli_expr_cols;
