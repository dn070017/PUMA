%% Using PUMA to infer gene regulatory network.
%%
%% 1. Reading in input data (expression data, motif prior, TF PPI data)
%% 2. Computing coexpression network
%% 3. Normalizing networks
%% 4. Running PUMA algorithm
%% 5. Writing out PUMA network (optional)
%%
%% Authors: cychen, marieke, kglass
%% Script adapted from Marieke's pLIONESSpcss.m, modified to run PUMA only.

disp(datestr(now));

%% ============================================================================
%% Set Program Parameters and Path
%% ============================================================================
% Run configuration to set parameter first (e.g., run('puma_config.m');)
fprintf('Input expression file: %s\n', exp_file);
fprintf('Input motif file: %s\n', motif_file);
fprintf('Input PPI file: %s\n', ppi_file);
fprintf('Output PUMA network: %s\n', puma_out);
fprintf('Output temp data: %s\n', save_temp);
fprintf('Alpha: %.2f\n', alpha);
addpath(lib_path);

%% ============================================================================
%% Read in Data
%% ============================================================================
disp('Reading in expression data!');
tic
    fid = fopen(exp_file, 'r');
    headings = fgetl(fid);
    n = length(regexp(headings, '\t'));
    frewind(fid);
    %Exp = textscan(fid, ['%s', repmat('%f', 1, n)], 'delimiter', '\t', 'CommentStyle', '#');
    Exp = textscan(fid, ['%s', repmat('%f', 1, n)], 'delimiter', '\t'); % tiny speed-up by not checking for comments
    fclose(fid);
    GeneNames = Exp{1};
    Exp = cat(2, Exp{2:end});
    [NumGenes, NumConditions] = size(Exp);
    fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions);
    Exp = Exp';  % transpose expression matrix from gene-by-sample to sample-by-gene
toc

disp('Reading in motif data!');
tic
    [TF, gene, weight] = textread(motif_file, '%s%s%f');
    TFNames = unique(TF);
    NumTFs = length(TFNames);
    [~,i] = ismember(TF, TFNames);
    [~,j] = ismember(gene, GeneNames);
    RegNet = zeros(NumTFs, NumGenes);
    RegNet(sub2ind([NumTFs, NumGenes], i, j)) = weight;
    fprintf('%d TFs and %d edges!\n', NumTFs, length(weight));
toc

disp('Reading in ppi data!');
tic
    TFCoop = eye(NumTFs);
    if(~isempty(ppi_file))
        [TF1, TF2, weight] = textread(ppi_file, '%s%s%f');
        [~,i] = ismember(TF1, TFNames);
        [~,j] = ismember(TF2, TFNames);
        TFCoop(sub2ind([NumTFs, NumTFs], i, j)) = weight;
        TFCoop(sub2ind([NumTFs, NumTFs], j, i)) = weight;
        fprintf('%d PPIs!\n', length(weight));
    end
toc
 
% Clean up variables to release memory
clear headings n GeneNames TF gene TF1 TF2 weight;

%% ============================================================================
%% Run PUMA
%% ============================================================================
disp('Computing coexpresison network:');
tic; GeneCoReg = Coexpression(Exp); toc;

disp('Normalizing Networks:');
tic
    RegNet = NormalizeNetwork(RegNet);
    GeneCoReg = NormalizeNetwork(GeneCoReg);
    TFCoop = NormalizeNetwork(TFCoop);
toc

if ~isempty(save_temp)
    disp('Saving the transposed expression matrix and normalized networks:');
    tic
        save(fullfile(save_temp, 'expression.transposed.mat'), 'Exp', '-v7.3');  % 2G+
        save(fullfile(save_temp, 'motif.normalized.mat'), 'RegNet', '-v6');  % fast
        save(fullfile(save_temp, 'ppi.normalized.mat'), 'TFCoop', '-v6');  % fast
    toc
end

clear Exp;  % Clean up Exp to release memory (for low-memory machine)

disp('Running PUMA algorithm:');
AgNet = PUMA(RegNet, GeneCoReg, TFCoop, alpha);

%% ============================================================================
%% Saving PUMA network output
%% ============================================================================
if ~isempty(puma_out)
    disp('Saving PUMA network!');
    tic
        [pathstr, name, ext] = fileparts(puma_out);
        switch ext
            case '.txt'
                save(puma_out, 'AgNet', '-ascii');
            case '.tsv'
                save(puma_out, 'AgNet', '-ascii', '-tabs');
            otherwise
                save(puma_out, 'AgNet', '-v6');
        end
    toc
end

disp('All done!');
disp(datestr(now));
