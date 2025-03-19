% Human HSCs

load('../Human-GEM/model/Human-GEM.mat');  % loads model as a structure named "ihuman"

% The second flag indicates if the model should be converted to gene symbols from ENSEMBL. This has to be decided at this point.
% Replace path/to/HumanGEM with your local path to the Human-GEM repo root.
% For use with animal models derived from Human-GEM, such as Mouse-GEM, both the model and paths needs to be replaced. Also, 
% the convert genes flag may be irrelevant depending on the if ENSEMBL genes are used in that model.
  prepData = prepHumanModelForftINIT(ihuman, false, '../Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt', '../Human-GEM/model/reactions.tsv');
  save('prepData.mat', 'prepData')

% replace 'my/path/' with the path on your system, or change to the directory containing the file
gtex_data = readtable('data/humanHSC_YO_TPM.txt'); %load from R
[~, n] = size(gtex_data);
numSamp = n-1; %the first columns are the genes in ENSEMBL format

% take a look at the first few rows and columns of the table
gtex_data(1:5, 1:5)

% setRavenSolver('gurobi') Setting this!
% extract the tissue and gene names
data_struct.genes = gtex_data{:, 1}; % gene names
data_struct.tissues = gtex_data.Properties.VariableNames(2:n); % sample (tissue) names
data_struct.levels = gtex_data{:, 2:n}; % gene TPM values
data_struct.threshold = 1;
models = cell(numSamp, 1);
for i = 1:numSamp
    disp(['Model: ' num2str(i)])
    mres = ftINIT(prepData, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+0'), false,  true,[]);
    mres.id = data_struct.tissues{i};
    models{i,1} = mres
end

save('models.mat', 'models')

