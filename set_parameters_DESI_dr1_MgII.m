% set_parameters: sets various parameters for the CIV detection 
% pipeline
% Designed for using DR16 spectra  
%flags for changes
extrapolate_subdla = 0; %0 = off, 1 = on
add_proximity_zone = 0;
integrate          = 1;
optTag = [num2str(integrate), num2str(extrapolate_subdla), num2str(add_proximity_zone)];

% physical constants
mgii_2796_wavelength = 2796.3543;		 % MgII transition wavelength  Å
mgii_2803_wavelength =  2803.5315; 		 % MgII transition wavelength  Å
speed_of_light = 299792458;                   % speed of light                     m s⁻¹

% converts relative velocity in km s^-1 to redshift difference
kms_to_z = @(kms) (kms * 1000) / speed_of_light;

% utility functions for redshifting
emitted_wavelengths = ...
    @(observed_wavelengths, z) (observed_wavelengths / (1 + z));

observed_wavelengths = ...
    @(emitted_wavelengths,  z) ( emitted_wavelengths * (1 + z));

%set figure rendering to avoid black plots
%set(0, 'DefaultFigureRenderer', 'painters')
%set(gcf, 'Color', 'white');


% %-------dr7
% % download Cooksey's dr7 spectra from this page: 
% % http://www.guavanator.uhh.hawaii.edu/~kcooksey/SDSS/CIV/index.html 
% % go to table: "SDSS spectra of the sightlines surveyed for C IV."
releaseTest='dr1';
releasePrior='dr7'; %Need to decide on a prior dataset


% file loading parameters
loading_min_lambda = 1555;          % range of rest wavelengths to load  Å, changed
loading_max_lambda = 2840;                    
% The maximum allowed is set so that even if the peak is redshifted off the end, the
% quasar still has data in the range

% preprocessing parameters
%z_qso_cut      = 2.15;                   % filter out QSOs with z less than this threshold
%z_qso_cut      = 1.7;                    % according to Cooksey z>1.7                      
min_num_pixels = 200;                     % minimum number of non-masked pixels NEED ADJUSTMENT

% NEW NORMALIZTION PARAMETERS
% range of rest wavelengths to use   Å
                                          % Continuum Fitting Redshift
                                          % Range
normalization_min_lambda_1 = 4150;   % Å  z < .6
normalization_max_lambda_1 = 4250;   % Å  z < .6

normalization_min_lambda_2 = 3020;   % Å  .6 < z < 1.0
normalization_max_lambda_2 = 3100;   % Å  .6 < z < 1.0

normalization_min_lambda_3 = 2150;   % Å  1.0 < z < 2.5   was 2150-2250
normalization_max_lambda_3 = 2250;   % Å  1.0 < z < 2.5

normalization_min_lambda_4 = 1420;   % Å  2.5 < z < 4.7
normalization_max_lambda_4 = 1500;   % Å  2.5 < z < 4.7

% null model parameters             XXX DESI has a 0.8A wavelength spacing
min_lambda         = loading_min_lambda+1;                   % range of rest wavelengths to       Å
max_lambda         = loading_max_lambda-1;                    %   model
dlambda            = 0.8;                    % separation of wavelength grid      Å
k                  = 20;                      % rank of non-diagonal contribution
max_noise_variance = 0.5^2;                   % maximum pixel noise allowed during model training
h                  = 2;                     % masking par to remove CIV region 
nAVG               = 0;                     % number of points added between two 
                                            % observed wavelengths to make the Voigt finer
% optimization parameters
minFunc_options =               ...           % optimization options for model fitting
    struct('MaxIter',     10000, ...
           'MaxFunEvals', 10000);

% C4 model parameters: parameter samples (for Quasi-Monte Carlo)
num_MgII_samples         = 20000;                  % number of parameter samples
alpha                    = 0.90;                    % weight of KDE component in mixture
uniform_min_log_nMgII     = 12.5;                   % range of column density samples    [cm⁻²]
uniform_max_log_nMgII     = 16.1;                   % from uniform distribution
fit_min_log_nMgII         = uniform_min_log_nMgII;                   % range of column density samples    [cm⁻²]
fit_max_log_nMgII         = 15.5;                   % from fit to log PDF


min_sigma                = 5e5;                   % cm/s -> b/sqrt(2) -> min Doppler par from Cooksey
max_sigma                = 115e5;                   % cm/s -> b/sqrt(2) -> max Doppler par from Cookse
vCut                     = 3000;                    % maximum cut velocity for CIV system 
RejectionSampling        = 0;
% model prior parameters
prior_z_qso_increase = kms_to_z(30000);       % use QSOs with z < (z_QSO + x) for prior

% instrumental broadening parameters
width = 3;                                    % width of Gaussian broadening (# pixels)
pixel_spacing = 1e-4;                         % wavelength spacing of pixels in dex

% DLA model parameters: absorber range and model
num_lines = 2;                                % number of members of CIV series to use

max_z_cut = kms_to_z(vCut);                   % max z_DLA = z_QSO - max_z_cut

max_z_MgII = @(z_qso, max_z_cut) ...         % determines maximum z_DLA to search
     z_qso - max_z_cut*(1+z_qso);
min_z_cut = kms_to_z(vCut);                   % min z_DLA = z_Ly∞ + min_z_cut
min_z_MgII = @(wavelengths, z_qso) ...         % determines minimum z_DLA to search
    max(min(wavelengths) / mgii_2796_wavelength - 1,                          ...
        observed_wavelengths(1310, z_qso) / mgii_2796_wavelength - 1);
% min_z_c4 = @(wavelengths, z_qso) ...         % determines minimum z_DLA to search
%      min(wavelengths) / civ_1548_wavelength - 1;
train_ratio =0.99;
sample_name = sprintf("N-%d-%d-Sigma-%d-%d-Num-%d",floor(fit_min_log_nMgII*100),floor(100*fit_max_log_nMgII), min_sigma,max_sigma, num_MgII_samples);
training_set_name = 'Seyfert_1555-2840';
testing_set_name = 'chi_2_test';
mkdir(sprintf('output/dr1/plots/%s',testing_set_name))

max_MgII = 10;
dv_mask = 750;%350; % (km/s)% base directory for all data 
%Change to something new on Cosmic
base_directory = 'output';
% utility functions for identifying various directories
distfiles_directory = @(release) ...
   sprintf('%s/%s/distfiles', base_directory, release);

%Need to change 
spectra_directory   = @(release)...
   sprintf('%s/%s/sas/dr1/sdss/spectro/redux/v5_13_0/spectra/lite', base_directory, release);

processed_directory = @(release) ...
   sprintf('%s/%s/processed', base_directory, release);

MgII_catalog_directory = @(name) ...
   sprintf('%s/C4_catalogs/%s/processed', base_directory, name);

   
% replace with @(varargin) (fprintf(varargin{:})) to show debug statements
% fprintf_debug = @(varargin) (fprintf(varargin{:}));
% fprintf_debug = @(varargin) ([]);

%---------DESI dr1-------------

% file_loader_DESI = @(tid)...
%     (read_spec_DESI(sprintf("/home/cosmic/DESI_data/%s.fits",tid)));
file_loader_DESI = @(tid)...
    (read_spec_DESI(sprintf("/home/cosmic/desi-env/sim_abs_mgii/%s.fits",tid)));