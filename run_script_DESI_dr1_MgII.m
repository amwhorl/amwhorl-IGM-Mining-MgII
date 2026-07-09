clear
fprintf('Setting paramters ...\n')
num_quasars = 1000;
cataloging = 0;
preloading =0;
sampling   = 0;
null_search = 0; 
plotting   = 1;
processing = 1;
merging    = 0;
EWer       =0;
pltP       =0;
CredInt    =0;
dv_mask    = 350; 
HPCC = 0;
voigtPrep = 0;
maskType = 1;
priorType = 1;
ind_S=1;
saving=1;
cores = 6;
SNR_threshhold = 0;
set_parameters_DESI_dr1_MgII
training_set_name
fprintf('Building catalogs ...\n')
if cataloging == 1
    build_catalog_DESI_dr1_MgII
end
variables_to_load= {'all_QSO_ID_dr1', 'all_RA_dr1', 'all_DEC_dr1', 'all_zqso_dr1'};
load(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_load{:});


%Prepare Voigt function (uses C)
fprintf('preparing voigt.c ...\n')
if voigtPrep == 1 
    cd minFunc_2012
    addpath(genpath(pwd))
    mexAll
    cd ..
    
    mex voigt_iP.c -lcerf
      
end


%Prepare prior catalog, choose dr7
catDR7 = load(sprintf('%s/catalog', processed_directory(releasePrior)));
filter_flagsDR7 = load(sprintf('%s/filter_flags', processed_directory(releasePrior)), ...       %%Error here, filter flags from prior release needed,
'filter_flags');
prior_ind = (filter_flagsDR7.filter_flags==0); 
all_z_MgII1 = catDR7.all_z_MgII1;
% all_REW_1548_DR7 = catDR12.all_EW1;
% all_REW_1550_DR7 = catDR12.all_EW2;


fprintf('Learning model ...\n')

variables_to_load = { 'max_noise_variance', ...
                   'minFunc_options', 'rest_wavelengths', 'mu', ...
                    'initial_M', 'M',  'log_likelihood', ...
                    };

load(sprintf('%s/learned_model-%s', processed_directory(releasePrior),...
                                     training_set_name), variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
% variables_to_load = {'offset_z_samples', 'offset_sigma_samples',...
%                      'log_nciv_samples', 'nciv_samples'};
if sampling==1
    generate_MgII_samples
end

variables_to_load = {'offset_z_samples', 'offset_sigma_samples'
                     'log_nMgII_samples', 'nMgII_samples'};

load(sprintf('%s/MgII_samples_%s.mat', processed_directory(releaseTest), sample_name), variables_to_load{:});

fprintf(sprintf('%d Samples are generated\n', num_MgII_samples));
    
% load preprocessed QSOs
fprintf('Preloading QSOs ...\n')
if preloading == 1
    preload_qsos_DESI_dr1_MgII
end
load(sprintf('%s/preloaded_qsos_%s', processed_directory(releaseTest), testing_set_name));


fprintf('preparing testing and prior indeces ...\n')
test_ind = (filter_flags==0);
if processing==1
    fprintf('processing QSO: %d to %d\n\n', ind_S, ind_S +  num_quasars-1);
    %parpool('local', cores);
    process_qsos_DESI_dr1_MgII
    
end

if merging==1
    mergeProcessesDr16
end

if EWer==1
    EW_dr16_voigt
end

if pltP  ==1
   pltPost;
end

if CredInt==1
    CI_sort
end