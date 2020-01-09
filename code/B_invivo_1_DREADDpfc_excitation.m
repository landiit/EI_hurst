% B_invivo_1_DREADDpfc_excitation.m

% add nonfractal toolbox to MATLAB path
addpath /Users/mvlombardo/Dropbox/matlab/nonfractal-master/m;

%% load in pheno data
rootpath = '/Users/mvlombardo/Dropbox/Manuscripts/AIMS_Hurst_Sex/reproAnalysis/DREADDexcitation';
datapath = '/Users/mvlombardo/Dropbox/data/ag_EI_DREADD_data/DREADDexcitiation/raw';
pheno_path = fullfile(rootpath,'pheno');
pdata = readtable(fullfile(pheno_path,'tidy_pheno_data.csv'),'delimiter',',');


%% Load in data
sublist = pdata.filename;
ntimepoints = 4040;
tidy_data = zeros(length(sublist),ntimepoints);
for i = 1:length(sublist)
    fname = fullfile(datapath,sublist{i});
    tmp_data = readtable(fname,'ReadVariableNames',false);
    tidy_data(i,:) = tmp_data.Var1';
end

for i = 1:ntimepoints
    volnames{i} = sprintf('vol%03d',i);
end
tidydata2write = cell2table(num2cell(tidy_data),'VariableNames',volnames);
writetable(tidydata2write,fullfile(rootpath,'tidy','tidy_data.csv'),'FileType','text','delimiter',',');


%% compute H on tidy_data_long
Hfilter = 'haar';
waveletType = 'dwt';
lb_param = [-0.5,0];
ub_param = [1.5,10];

window_length = 512;
start_idx = 1:((ntimepoints)-(window_length-1));
end_idx = 512:(ntimepoints);

Hwin = zeros(size(tidy_data,1),length(start_idx));
parfor iwindow = 1:length(start_idx)
    disp(iwindow);
    idx2use = start_idx(iwindow):end_idx(iwindow);
    data2use = tidy_data(:,idx2use);
    Hwin(:,iwindow) = bfn_mfin_ml(data2use', ...
        'filter',Hfilter, ...
        'lb',lb_param,'ub',ub_param, ...
        'wavelet',waveletType, ...
        'verbose',false);
end % for iwindow

for i = 1:size(Hwin,2)
    volnames3{i} = sprintf('window_%04d',i);
end

Hwin_tab = cell2table(num2cell(Hwin),'VariableNames',volnames3);
pdata_final = [pdata Hwin_tab];

fname2write = fullfile(rootpath,'pheno','pheno_data+Hwin_dreaddexcitation.csv');
writetable(pdata_final,fname2write,'FileType','text','delimiter',',');
