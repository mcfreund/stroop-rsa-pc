%% TODO
% add outlier TRs to extraregressors
% crossvalidate using pc50 items only
% test: running glmdenoise on volume, and adding noisepcs as extraregressors


%% setup

% Add path to GLMsingle and fracridge

addpath(genpath('/data/nil-external/ccp/freund/GLMsingle/matlab'));
addpath('/data/nil-external/ccp/freund/fracridge/matlab')

% Start fresh

clear
clc
close all

% constants

opt = struct('wantmemoryoutputs', [0 1 0 1]);
opt.wantglmdenoise = 1;
opt.sessionindicator = [1 1 2 2 3 3];
subjlist = 'wave1_unrel';
glmname = 'glmsingle_wave1';
variables = {'target' 'distractor' 'congruency'};
stimdur = 1;  % stimulus duration (s)
tr = 1.2; % sec
n_vertex = 20484;
idx = {1:216 217:432 433:672};  % bas, pro, rea, both runs
cellnames = {'onoff' 'fithrf' 'glmdenoise', 'rr'};
sessions = {'baseline' 'proactive' 'reactive'};
sess = {'Bas' 'Pro' 'Rea'};
runs = {'run1' 'run2'};
runencs = {'1_AP' '2_PA'};
n_trial = [216 216 240];  % bas pro rea, both runs

proj_path = '/data/nil-external/ccp/freund/stroop-rsa-pc';
path_bluearc = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/';
outputdir = '/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms';
if opt.wantglmdenoise == 1
    path_res_suffix = '_denoise'
else
    path_res_suffix = ''
end

fname_subjs = fullfile(proj_path, 'out', append('subjlist_', subjlist, '.txt'));
subjs = regexp(fileread(fname_subjs), '\r?\n', 'split');



%% get data, design matrices, and 'extra regressors'

subjs=subjs(1:5);

for subj_i = 1:length(subjs)
    
    subj = subjs{subj_i};
    path_design = fullfile(proj_path, 'out', 'glms', subj, 'wave1', 'RESULTS', 'Stroop', glmname);

    for var_i = 1:length(variables)
        
        variable = variables{var_i};
        data = cell(length(sessions), length(runs));
        design = cell(length(sessions), length(runs));
        extreg = cell(length(sessions), length(runs));

        for ses_i = 1:length(sessions)

            session = sessions{ses_i};
            ses = sess{ses_i};

            % read design matrices
            fname_design = fullfile(path_design, append('design_', variable, '.h5'));
            x = h5read(fname_design, append('/', session));
            design(ses_i, 1) = {x(:, :, 1)};
            design(ses_i, 2) = {x(:, :, 2)};
            
            path_extreg = fullfile(proj_path, 'out', 'glms', subj, 'wave1', 'INPUT_DATA', 'Stroop', session);
            
            for run_i = 1:2
                run = runs{run_i};
                runenc = runencs{run_i};
                
                % read extraregressors (TRs to censor)
                fname_fdmask = fullfile(path_extreg, append('movregs_FD_mask_run' + string(run_i) + '.txt'));
                extreg_i = 1 - readmatrix(fname_fdmask{1});
                if sum(extreg_i) == 0  % if no bad TRs, set extrareg to empty cell
                    extreg_i = [];
                else
                    idx_bad_trs = find(extreg_i);
                    extreg_i(idx_bad_trs) = idx_bad_trs;  % set bad TRs to unique value
                    extreg_i = dummyvar(extreg_i + 1);  % buld dummy matrix
                    extreg_i = extreg_i(:, idx_bad_trs + 1);  %% extract cols for bad TR dummy regressors
                    extreg(ses_i, run_i) = {extreg_i};
                end    

                % read giftis and concatenate across hemispheres
                path_data = fullfile(path_bluearc, subj, 'INPUT_DATA', 'Stroop', session);
                fnames_data = fullfile(path_data, append('lpi_scale_tfMRI_Stroop', ses, runenc, '_', {'L' 'R'}, '.func.gii'));
                giiL = gifti(fnames_data{1});
                giiR = gifti(fnames_data{2});
                data(ses_i, run_i) = {[giiL.cdata; giiR.cdata]};
            end

        end
        
        % reshape to 1 x n_fold arrays
        design = reshape(design', 1, []);
        data = reshape(data', 1, []);
        extreg = reshape(extreg', 1, []);
        opt.extraregressors = extreg;

        % fit models
        path_res = fullfile(path_design, append(variable, path_res_suffix));
        [results] = GLMestimatesingletrial(design, data, stimdur, tr, path_res, opt);

        % extract betas, reshape, and save in separate hdf5 file
        idx_res = find(opt.wantmemoryoutputs);
        for b_i = idx_res
            b = results{b_i};
            b = b.modelmd;
            if ~isequal(size(b), [n_vertex 1 1 672])
                error("bad output dims")
            end
            file_name = fullfile(path_res, append('betas_', cellnames{b_i}, '.h5'));
            if isfile(file_name)
                inf = h5info(file_name);
                datasets = {inf.Datasets.Name};
            else
                datasets = {};
            end
            for ses_i = 1:length(sessions)
                b_ses = reshape(b(:, 1, 1, idx{ses_i}), n_vertex, []);
                dset_name = append('/', sessions{ses_i});
                if ~any(strcmp(datasets, sessions{ses_i}))
                    h5create(file_name, dset_name, size(b_ses));
                end
                h5write(file_name, dset_name, b_ses);
            end
        end
               


    end
end

