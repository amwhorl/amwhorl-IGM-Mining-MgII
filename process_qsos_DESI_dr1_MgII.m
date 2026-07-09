    % process_qsos: run CIV detection algorithm on specified objects

% Averaging     -> yes
% Single model  -> yes
% Multi CIV     -> yes
% Multi Singlet -> no
%
% removing CIV found region after each run 
% also removing the map Z_civ
% load C4 catalog
Full_catalog = ...
    load(sprintf('%s/catalog', processed_directory(releasePrior)));
% train_ind -> those LOSs without any problem (filter_flag==0) that does not have
% Civ and useful for training null model (a model without Civ absorption line)
% prior_ind -> those LOSs with Civ absorption and in the half part of test
% test_ind -> second half without any filter flag for testing null and absorption model
% on the LOSs that we know have Civ or not. So, we asses our algorithm in this way.

if (ischar(prior_ind))
    prior_ind = eval(prior_ind);
end

% My prior_ind here is already those OK sight of lines that have CIV
prior.z_qsos  = Full_catalog.all_zqso(prior_ind);
prior.MgII_ind = prior_ind;
prior.z_MgII = Full_catalog.all_z_MgII3(prior_ind);

% filter out CIVs from prior catalog corresponding to region of spectrum below
% Ly-alpha QSO rest. In Roman's code, for detecting DLAs, instead of
% Ly-alpha they have Lyman-limit
for i = size(prior.z_MgII)
    if (observed_wavelengths(mgii_2796_wavelength , prior.z_MgII(i)) < ...
            observed_wavelengths(min_lambda, prior.z_qsos(i)))
        prior.MgII_ind(i) = false;
    end
end
prior = rmfield(prior, 'z_MgII');

% enable processing specific QSOs via setting to_test_ind
if (ischar(test_ind))
    test_ind = eval(test_ind);   %Test index comes from the inverse of filter_flags
end

all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
all_sigma_pixel    =    all_sigma_pixel(test_ind);
z_qsos             =     all_zqso_dr1(test_ind);
num_quasars        =                numel(z_qsos);
% REW_PM             =          all_EW1(test_ind,:);
% errREW_PM             =       all_errEW1(test_ind,:);



% preprocess model interpolants
% griddedInterpolant does an interpolation and gives a function handle
% based on the grided data. If the data is 2xD, {x,y} are the the same size as
% row  and columns. M is like M=f(x,y) and is like a matrix with each element
% M(i,j) = f(x(i), y(j))
mu_interpolator = ...
    griddedInterpolant(rest_wavelengths,        mu,        'linear');
M_interpolator = ...
    griddedInterpolant({rest_wavelengths, 1:k}, M,         'linear');

% initialize results with nan
min_z_MgIIs                   = nan(num_quasars, 1);
max_z_MgIIs                   = nan(num_quasars, 1);
log_priors_no_MgII            = nan(num_quasars, max_MgII);
log_priors_MgII               = nan(num_quasars, max_MgII);
log_likelihoods_no_MgII       = nan(num_quasars, max_MgII);
sample_log_likelihoods_MgIIL1 = nan(num_quasars, num_MgII_samples);
sample_log_likelihoods_MgIIL2 = nan(num_quasars, num_MgII_samples, max_MgII);
log_likelihoods_MgIIL1        = nan(num_quasars, max_MgII);
log_likelihoods_MgIIL2        = nan(num_quasars, max_MgII);
log_posteriors_no_MgII        = nan(num_quasars, max_MgII);
log_posteriors_MgIIL1         = nan(num_quasars, max_MgII);
log_posteriors_MgIIL2         = nan(num_quasars, max_MgII);
map_N_MgIIL1                  = nan(num_quasars, max_MgII);
map_N_MgIIL2                  = nan(num_quasars, max_MgII);
map_z_MgIIL1                  = nan(num_quasars, max_MgII);
map_z_MgIIL2                  = nan(num_quasars, max_MgII);
map_sigma_MgIIL1              = nan(num_quasars, max_MgII);
map_sigma_MgIIL2              = nan(num_quasars, max_MgII);
p_MgII                        = nan(num_quasars, max_MgII);
p_MgIIL1                      = nan(num_quasars, max_MgII);
p_no_MgII                     = nan(num_quasars, max_MgII);
REW_2796                      = nan(num_quasars, max_MgII);
REW_2803                      = nan(num_quasars, max_MgII);
num_pixel_MgII                = nan(num_quasars, max_MgII, 2);
all_B                         = nan(num_quasars, max_MgII);
fit_chi2                      = nan(num_quasars,1);
MgII_samples                  = (max_sigma-min_sigma)*offset_sigma_samples + min_sigma;
sigma_MgII_samples            = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
ID                            = all_QSO_ID_dr1(test_ind);
% plt_count=0;
% FN_IDs = importdata('FN-list.csv');
% in_Kathy_FN_list = ismember(ID, FN_IDs);
% N_MgII_test = all_N_MgII(test_ind,:);
z_PM_test = all_z_MgII1(test_ind,:);
z_PM_prior = all_z_MgII1(prior_ind,:);
j0=0;

qso_NaN_replacedL2_ind = [];
qso_NaN_replacedL1_ind = [];

all_SNR = nan(numel(all_wavelengths),1);

for quasar_ind = 1:numel(all_wavelengths)

    % if z_PM_test(quasar_ind,1)==-1
    %     continue
    % end

    tic;
    z_qso = z_qsos(quasar_ind);
    fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...', ...
                              quasar_ind, num_quasars, z_qso);
    
    this_wavelengths    =    all_wavelengths{quasar_ind};
    this_flux           =           all_flux{quasar_ind}; 
    this_noise_variance = all_noise_variance{quasar_ind};
    this_pixel_mask     =     all_pixel_mask{quasar_ind};
    this_sigma_pixel    =     all_sigma_pixel{quasar_ind};


    % convert to QSO rest frame
     this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);

    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
       (this_rest_wavelengths <= max_lambda) & (this_sigma_pixel>0);% &...
         %~(abs(this_rest_wavelengths- 1549.48)< 15)& ... % CIV emission line masking --> 30A comes from typical FWHM of CIV emission in SDSS Cite{Monadi & Bird-2022}
         %~(abs(this_rest_wavelengths- 1908.8)< 12) & ... % CIII masking here
         %~(abs(this_rest_wavelengths - 2799.94) < 12) & ... % MgII emission Line Masking
         %~(abs(this_rest_wavelengths- 1393.8)< 6); % SIV emission line masking --> A is FWHM according to

    % keep complete copy of equally spaced wavelengths for absorption computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);
    this_unmasked_sigma_pixel = this_sigma_pixel(unmasked_ind); % avoiding mismathed sizes for padded variavles and NaNs
%     % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
%     % in read_spec_DESI.m
    ind                   =  unmasked_ind & (~this_pixel_mask);
    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_noise_variance   =   this_noise_variance(ind);
    this_sigma_pixel      =      this_sigma_pixel(ind);
    
    % SNR Filter
    SNR = median(this_flux./sqrt(this_noise_variance));
    all_SNR(quasar_ind) = SNR;
    if SNR<SNR_threshhold
        continue
    end

 % MgII existence prior
    less_ind = (prior.z_qsos < (z_qso + prior_z_qso_increase));
    less_systems = z_PM_prior(less_ind,:);

    this_num_quasars = nnz(less_ind);
    this_p_MgII(1) = nnz(less_systems(:,1)>0 )/this_num_quasars;  % at least 1
    for i=2:max_MgII
        this_p_MgII(i) = nnz(less_systems(:,i)>0 )/nnz(less_systems(:,i-1)>0); % at least n given at least n-1
        if (this_p_MgII(i-1)==0)
        this_p_MgII(i) = 0;
        end

    end

    fprintf('\n');
    for i = 1:max_MgII
        fprintf(' ...     p(%i  MgIIs | z_QSO)       : %0.3f\n', i, this_p_MgII(i));
            log_priors_no_MgII(quasar_ind, i) = ...
                log(1 - this_p_MgII(i));
        fprintf(' ...     p(no MgII  | z_QSO)       : %0.3f\n', exp(log_priors_no_MgII(quasar_ind, i)) );
    end

    log_priors_MgII(quasar_ind,:) = log(this_p_MgII(:));

    % interpolate model onto given wavelengths
    this_mu = mu_interpolator( this_rest_wavelengths);
    this_M  =  M_interpolator({this_rest_wavelengths, 1:k});


    % Define Search Range
    min_z_MgIIs(quasar_ind) = min_z_MgII(this_wavelengths, z_qso);
    max_z_MgIIs(quasar_ind) = max_z_MgII(z_qso, max_z_cut);

    % Constrain Samples to Search Range
    sample_z_MgII = ...
        min_z_MgIIs(quasar_ind) +  ...
        (max_z_MgIIs(quasar_ind) - min_z_MgIIs(quasar_ind)) * offset_z_samples;


    % Temperature samples

    % ensure enough pixels are on either side for convolving with
    % instrument profile

   % building a finer wavelength and mask arrays 
   % by adding the mean of ith and ith +1 element

    padded_wavelengths_fine = ...
        [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing/(nAVG+1), ...
        log10(min(this_unmasked_wavelengths)) - pixel_spacing/(nAVG+1),...
        width)';...
        finer(this_unmasked_wavelengths, nAVG)';...
        logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing/(nAVG+1),...
        log10(max(this_unmasked_wavelengths)) + width * pixel_spacing/(nAVG+1),...
        width)'...
        ];

      padded_sigma_pixels_fine = ...
        [this_unmasked_sigma_pixel(1)*ones(width,1);...
        finer(this_unmasked_sigma_pixel, nAVG)';...
        this_unmasked_sigma_pixel(end)*ones(width,1)];

        % when broadening is off
        % padded_wavelengths = this_unmasked_wavelengths;

    % [mask_ind] to retain only unmasked pixels from computed absorption profile
    % this has to be done by using the unmasked_ind which has not yet
    % been applied this_pixel_mask.
    ind = (~this_pixel_mask(unmasked_ind));
    
    % compute probabilities under DLA model for each of the sampled
    % (normalized offset, log(N HI)) pairsl
    lenW_unmasked = length(this_unmasked_wavelengths);
    ind_not_remove = true(size(this_flux));
    absorptionL2_all =1;
    absorptionL1_all =1;

    % Calculate Goodness of Continuum Fit
    bin_size = 20;
    divisable_array_legnth = (length(this_mu) - mod(length(this_mu),20));
    difference_bins = zeros(1,divisable_array_legnth/bin_size + 1); 
    this_difference_weighted = (this_mu - this_flux);
    % divide flux and continuum arrays into bins of [bin_size] length and
    % calculate the mean of each of those bins
    idx_start = 1;
    idx_stop = bin_size;
    for binning_idx = 1:(length(difference_bins) -1)
        difference_bins(binning_idx) = mean(this_difference_weighted(idx_start:idx_stop));
        idx_start = idx_start + bin_size;
        idx_stop = idx_stop + bin_size;
    end
    
    % Custom bins for the tails of both arrays
    difference_bins(end)   = mean(this_difference_weighted(idx_start:end));
    % Calculate Chi square of difference bins
    Chi2 = sum(difference_bins.^2);
    fit_chi2(quasar_ind) = Chi2;

    if null_search == 1 % limits to 1 search if looking for no-absorber spectra to inject with sumilated absorbers
        max_MgII = 1;
    end
    for num_MgII=1:max_MgII
        fprintf('num_MgII:%d\n',num_MgII);
        this_z_2796 = (this_wavelengths / mgii_2796_wavelength) - 1;
        this_z_2803 = (this_wavelengths / mgii_2803_wavelength) - 1;
        if(num_MgII>1)
            if((p_MgII(quasar_ind, num_MgII-1)>p_MgIIL1(quasar_ind, num_MgII-1)) & ...
                (p_MgII(quasar_ind, num_MgII-1)>p_no_MgII(quasar_ind, num_MgII-1)))
                ind_not_remove = ind_not_remove  & ...
                    (abs(this_z_2796 - map_z_MgIIL2(quasar_ind, num_MgII-1))>kms_to_z(dv_mask)*(1+map_z_MgIIL2(quasar_ind, num_MgII-1))) & ...
                    (abs(this_z_2803 - map_z_MgIIL2(quasar_ind, num_MgII-1))>kms_to_z(dv_mask)*(1+map_z_MgIIL2(quasar_ind, num_MgII-1)));
            end

            if((p_MgIIL1(quasar_ind, num_MgII-1)>p_MgII(quasar_ind, num_MgII-1)) & ...
                        (p_MgIIL1(quasar_ind, num_MgII-1)>p_no_MgII(quasar_ind, num_MgII-1)))
                ind_not_remove = ind_not_remove  & ...
                (abs(this_z_2796 - map_z_MgIIL1(quasar_ind, num_MgII-1))>kms_to_z(dv_mask)*(1+map_z_MgIIL1(quasar_ind, num_MgII-1)));
            end

            if((p_no_MgII(quasar_ind, num_MgII-1)>=p_MgII(quasar_ind, num_MgII-1)) & ...
                (p_no_MgII(quasar_ind, num_MgII-1)>=p_MgIIL1(quasar_ind, num_MgII-1)))

                fprintf('No more than %d MgIIs in this spectrum.', num_MgII-1)
                break;
            end

        end
        log_likelihoods_no_MgII(quasar_ind, num_MgII) = ...
        log_mvnpdf_low_rank(this_flux(ind_not_remove), this_mu(ind_not_remove),...
        this_M(ind_not_remove, :), this_noise_variance(ind_not_remove));
        % fprintf('S(this_M(ind_not_remove))=%d-%d\n', size(this_M(ind_not_remove, :)));
        log_posteriors_no_MgII(quasar_ind, num_MgII) = ...
            log_priors_no_MgII(quasar_ind, num_MgII) + log_likelihoods_no_MgII(quasar_ind, num_MgII);

        fprintf(' ... log p(D | z_QSO, no CIV)     : %0.2f\n', ...
        log_likelihoods_no_MgII(quasar_ind, num_MgII));
        fprintf(' ... log p(no CIV | D, z_QSO)     : %0.2f\n', ...
        log_posteriors_no_MgII(quasar_ind, num_MgII));
        parfor i = 1:num_MgII_samples
            % Limitting red-shift in the samples

            num_lines=2;

            absorptionL2_fine = voigt_iP(padded_wavelengths_fine, sample_z_MgII(i), ...
            nMgII_samples(i),num_lines, sigma_MgII_samples(i), padded_sigma_pixels_fine);
            % absorptionL2_fine_NaNs = isnan(absorptionL2_fine); % replaces NaN's created by voigt profile, maybe temp fix
            % absorptionL2_fine(absorptionL2_fine_NaNs) = 1;
            % % get w_r for this sample  
            % if (nnz(absorptionL2_fine_NaNs) > 0) && (i == 1) && (num_MgII == 1) % Keep track of spectra that had NaNs replaced
            %     qso_NaN_replacedL2_ind = [qso_NaN_replacedL2_ind, quasar_ind];
            % end
            % average fine absorption and shrink it to the size of original array
            % as large as the unmasked_wavelengths

            absorptionL2 = Averager(absorptionL2_fine, nAVG, lenW_unmasked);
            absorptionL2 = absorptionL2(ind);
            MgII_muL2     = this_mu     .* absorptionL2;
            MgII_ML2      = this_M      .* absorptionL2;
        
            sample_log_likelihoods_MgIIL2(quasar_ind, i, num_MgII) = ...
            log_mvnpdf_low_rank(this_flux(ind_not_remove),...
                                MgII_muL2(ind_not_remove),...
                                MgII_ML2(ind_not_remove, :), ...
                                this_noise_variance(ind_not_remove));

            num_lines=1;
            % absorptionL1_fine = voigt_iP(finer(padded_wavelengths, nAVG), sample_z_MgII(i), ...
            % nMgII_samples(i),num_lines, sigma_MgII, finer(this_sigma_pixel, nAVG));

            absorptionL1_fine = voigt_iP(padded_wavelengths_fine, sample_z_MgII(i), ...
            nMgII_samples(i),num_lines, sigma_MgII_samples(i), padded_sigma_pixels_fine);
            % absorptionL1_fine_NaNs = isnan(absorptionL1_fine); % replaces NaN's created by voigt profile, maybe temp fix
            % absorptionL1_fine(absorptionL1_fine_NaNs) = 1;
            % if (nnz(absorptionL1_fine_NaNs) > 0) && (i == 1) && (num_MgII == 1) % Keep track of spectra that had NaNs replaced
            %     qso_NaN_replacedL1_ind = [qso_NaN_replacedL1_ind, quasar_ind];
            % end
            % average fine absorption and shrink it to the size of original array
            % as large as the unmasked_wavelengths

            absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
            absorptionL1 = absorptionL1(ind);
            MgII_muL1     = this_mu     .* absorptionL1;
            MgII_ML1      = this_M      .* absorptionL1;
            sample_log_likelihoods_MgIIL1(quasar_ind, i) = ...
            log_mvnpdf_low_rank(this_flux(ind_not_remove),...
            MgII_muL1(ind_not_remove), MgII_ML1(ind_not_remove, :), ...
            this_noise_variance(ind_not_remove));
 

        end

        % compute sample probabilities and log likelihood of DLA model in
        % numerically safe manner for one line
        max_log_likelihoodL1 = max(sample_log_likelihoods_MgIIL1(quasar_ind, :));
        sample_probabilitiesL1 = ...
            exp(sample_log_likelihoods_MgIIL1(quasar_ind, :) - ...
            max_log_likelihoodL1);
        log_likelihoods_MgIIL1(quasar_ind, num_MgII) = ...
            max_log_likelihoodL1 + log(mean(sample_probabilitiesL1));% ...
            % - log(num_MgII_samples)*(num_MgII-1);

        log_posteriors_MgIIL1(quasar_ind, num_MgII) = ...
        log_priors_MgII(quasar_ind, num_MgII) + log_likelihoods_MgIIL1(quasar_ind, num_MgII);

        fprintf(' ... log p(D | z_QSO,    L1)     : %0.2f\n', ...
            log_likelihoods_MgIIL1(quasar_ind, num_MgII));
        fprintf(' ... log p(L1 | D, z_QSO)        : %0.2f\n', ...
            log_posteriors_MgIIL1(quasar_ind, num_MgII));

        % compute sample probabilities and log likelihood of DLA model in
        % numerically safe manner for  doublet 
        max_log_likelihoodL2 = max(sample_log_likelihoods_MgIIL2(quasar_ind, :, num_MgII));
        sample_probabilitiesL2 = ...
            exp(sample_log_likelihoods_MgIIL2(quasar_ind, :, num_MgII)  ... 
            - max_log_likelihoodL2);
        log_likelihoods_MgIIL2(quasar_ind, num_MgII) = ...
            max_log_likelihoodL2 + log(mean(sample_probabilitiesL2));%...
            % - log(num_MgII_samples)*(num_MgII-1);


        log_posteriors_MgIIL2(quasar_ind, num_MgII) = ...
            log_priors_MgII(quasar_ind, num_MgII) + log_likelihoods_MgIIL2(quasar_ind, num_MgII);

        fprintf(' ... log p(D | z_QSO,    CIV)     : %0.2f\n', ...
            log_likelihoods_MgIIL2(quasar_ind, num_MgII));
        fprintf(' ... log p(CIV | D, z_QSO)        : %0.2f\n', ...
            log_posteriors_MgIIL2(quasar_ind, num_MgII));
        [~, maxindL1] = max(sample_log_likelihoods_MgIIL1(quasar_ind, :));
        map_z_MgIIL1(quasar_ind, num_MgII )    = sample_z_MgII(maxindL1);        
        map_N_MgIIL1(quasar_ind, num_MgII)  = log_nMgII_samples(maxindL1);
        map_sigma_MgIIL1(quasar_ind, num_MgII)  = sigma_MgII_samples(maxindL1);
        % fprintf('L1\nmap(N): %.2f, map(z_MgII): %.2f, map(b/1e5): %.2f\n',map_N_MgIIL1(quasar_ind, num_MgII),...
            % map_z_MgIIL1(quasar_ind, num_MgII), map_sigma_MgIIL1(quasar_ind, num_MgII)/1e5);



        [~, maxindL2] = max(sample_log_likelihoods_MgIIL2(quasar_ind, :, num_MgII));
        map_z_MgIIL2(quasar_ind, num_MgII)    = sample_z_MgII(maxindL2);        
        map_N_MgIIL2(quasar_ind, num_MgII)  = log_nMgII_samples(maxindL2);
        map_sigma_MgIIL2(quasar_ind, num_MgII)  = sigma_MgII_samples(maxindL2);
        % fprintf('L2\nmap(N): %.2f, map(z_MgII): %.2f, map(b/1e5): %.2f\n',...
        % map_N_MgIIL2(quasar_ind, num_MgII), map_z_MgIIL2(quasar_ind, num_MgII),...
        % map_sigma_MgIIL2(quasar_ind, num_MgII)/1e5);

        max_log_posteriors = max([log_posteriors_no_MgII(quasar_ind, num_MgII), log_posteriors_MgIIL1(quasar_ind, num_MgII), log_posteriors_MgIIL2(quasar_ind,num_MgII)], [], 2);

        model_posteriors = ...
                exp(bsxfun(@minus, ...           
                [log_posteriors_no_MgII(quasar_ind, num_MgII), log_posteriors_MgIIL1(quasar_ind, num_MgII), log_posteriors_MgIIL2(quasar_ind, num_MgII)], ...
                max_log_posteriors));
        model_posteriors = ...
        bsxfun(@times, model_posteriors, 1 ./ sum(model_posteriors, 2));

        p_no_MgII(quasar_ind, num_MgII) = model_posteriors(1);
        p_MgIIL1(quasar_ind, num_MgII)  = model_posteriors(2);
        p_MgII(quasar_ind, num_MgII)    = 1 - p_no_MgII(quasar_ind, num_MgII) -...
                                    p_MgIIL1(quasar_ind, num_MgII);

        c4_pixel_ind1 = abs(this_wavelengths - (1+map_z_MgIIL2(quasar_ind, num_MgII))*mgii_2796_wavelength)<3;
        c4_pixel_ind2 = abs(this_wavelengths - (1+map_z_MgIIL2(quasar_ind, num_MgII))*mgii_2803_wavelength)<3;
        num_pixel_MgII(quasar_ind, num_MgII, 1) = nnz(c4_pixel_ind1);
        num_pixel_MgII(quasar_ind, num_MgII, 2) = nnz(c4_pixel_ind2);
        fprintf('CIV pixels:[%d, %d]\n', num_pixel_MgII(quasar_ind, num_MgII, :)); 

        % fprintf('s(fine_L1)-%d-%d\n', size(absorptionL1_fine)) 
        % fprintf('s(unmasked)-%d-%d\n', size(this_unmasked_wavelengths))                                    
        % fprintf('s(emitted_finer_unmasked)-%d-%d\n', size(emitted_wavelengths(finer(this_unmasked_wavelengths, nAVG), z_qso)))                                    


            % REW Calculations for both lines


            aL1_fine = voigt_iP(padded_wavelengths_fine,... % singlet absorbtion profile -> whichever line is most probable
                                         map_z_MgIIL2(quasar_ind, num_MgII), ...
                                         10^map_N_MgIIL2(quasar_ind, num_MgII), 1,...
                                         map_sigma_MgIIL2(quasar_ind, num_MgII), ...
                                         padded_sigma_pixels_fine);
            % aL1_fine_NaNs = isnan(aL1_fine); % replaces NaN's created by voigt profile, maybe temp fix
            % aL1_fine(aL1_fine_NaNs) = 1;
            aL1 = Averager(aL1_fine, nAVG, lenW_unmasked);
            aL1 = aL1(ind); 

            aL2_fine = voigt_iP(padded_wavelengths_fine,... % doublet absorbtion profile
                                         map_z_MgIIL2(quasar_ind, num_MgII), ...
                                         10^map_N_MgIIL2(quasar_ind, num_MgII), 2,...
                                         map_sigma_MgIIL2(quasar_ind, num_MgII), ...
                                         padded_sigma_pixels_fine);
            % aL2_fine_NaNs = isnan(aL2_fine); % replaces NaN's created by voigt profile, maybe temp fix
            % aL2_fine(aL2_fine_NaNs) = 1;
            aL2 = Averager(aL2_fine, nAVG, lenW_unmasked);
            aL2 = aL2(ind);
            % checks to ensure that we are calculating the correct REW for
            % the correct line.
            [~,abs_minL1] = min(aL1);
            [~,abs_minL2] = min(aL2./aL1);
            if abs_minL1 > abs_minL2 % If the singlet line is a higher ind, its at a higher wavelegnth, and so its the 2803 line
                REW_2803(quasar_ind, num_MgII) = trapz(this_unmasked_wavelengths(ind), 1-aL1)/(1+map_z_MgIIL2(quasar_ind, num_MgII));
                REW_2796(quasar_ind, num_MgII) = trapz(this_unmasked_wavelengths(ind), 1-(aL2./aL1))/(1+map_z_MgIIL2(quasar_ind, num_MgII));
            else % Vice versa -> the singlet line was the 2796 line
                REW_2796(quasar_ind, num_MgII) = trapz(this_unmasked_wavelengths(ind), 1-aL1)/(1+map_z_MgIIL2(quasar_ind, num_MgII));
                REW_2803(quasar_ind, num_MgII) = trapz(this_unmasked_wavelengths(ind), 1-(aL2./aL1))/(1+map_z_MgIIL2(quasar_ind, num_MgII));
            end
            fprintf('REW(%d,%d)=%e\n', quasar_ind, num_MgII, REW_2796(quasar_ind, num_MgII));


         if(plotting==1) 
            % plotting

            this_ID = ID{quasar_ind};
            max_log_posteriors = max([log_posteriors_no_MgII(quasar_ind, num_MgII),...
                log_posteriors_MgIIL1(quasar_ind, num_MgII), ...
                log_posteriors_MgIIL2(quasar_ind,num_MgII)], [], 2);

            mu_post = mean_posterior_continuum(this_flux, this_mu,...
        this_M, this_noise_variance);

            num_lines=1;
            absorptionL1_fine= voigt_iP(padded_wavelengths_fine,...
                            map_z_MgIIL1(quasar_ind, num_MgII),1.2*(10^map_N_MgIIL1(quasar_ind, num_MgII)),...
                            num_lines, map_sigma_MgIIL2(quasar_ind, num_MgII), padded_sigma_pixels_fine);
            absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
            absorptionL1 = absorptionL1(ind);
            MgII_muL1    = mu_post     .* absorptionL1;

            MgII_muL2    = mu_post     .* aL2;
            
            if (null_search == 0) || ((null_search == 1) && (p_no_MgII(quasar_ind, num_MgII) > 0.5))
                % Different ways of quanitifying how good the continuum/doublet fit is
                lambda_cut_range = (this_wavelengths > ((map_z_MgIIL2(quasar_ind,num_MgII) + 1)*2796.4) -100) & ...
                 (this_wavelengths < ((map_z_MgIIL2(quasar_ind,num_MgII) + 1)*2796.4) + 100);
                normalization_cut_range = (this_wavelengths > (2150*(z_qso + 1))) & ...
                 (this_wavelengths < (2250*(z_qso + 1)));
    
  
                % ttl = sprintf('ID:%s, zQSO:%.2f, P(CIV)=%.2f, P(S)=%.2f, z_{CIV}=%.6f\nz_{PM}=[%.4f,%.4f,%.4f,%.4f]\nREW_{PM}=[%.3f,%.3f,%.3f,%.3f]\n errREW_{PM}=[%.3f,%.3f,%.3f,%.3f], REW(GP)=%.3f, err(GP)=%.3f',  ...
                %     this_ID, z_qso, p_MgII(quasar_ind, num_MgII), p_MgIIL1(quasar_ind, num_MgII),  map_z_MgIIL2(quasar_ind, num_MgII), ...
                %     z_PM_test(quasar_ind,1:4),...
                %     REW_PM(quasar_ind,1:4), errREW_PM(quasar_ind,1:4),...
                %     REW_2796_dr7_flux(quasar_ind, num_MgII), ErrREW_1548_flux(quasar_ind, num_MgII));
                % DZ = abs(z_PM_test(quasar_ind, 1:4) - map_z_MgIIL2(quasar_ind,num_MgII));
                % 
                % dv = DZ./(1+z_PM_test(quasar_ind, 1:4))*speed_of_light/1e3;
                
                ttl = sprintf(['ID:%s, zQSO:%.2f\n P(MgII)=%.2f, P(S)=%.2f, z_{MgII}=%.6f, S/N=%.2f, sigma=%.3f, ' ...
                    'N=%.2f\n W_{0}^{2796}=%.3f, W_{0}^{2803}=%.3f, \\chi^{2}=%.4f, log(P(D|M_{N}))=%.4f'],  ...
                    this_ID, z_qso, p_MgII(quasar_ind, num_MgII), p_MgIIL1(quasar_ind, num_MgII),  ...
                    map_z_MgIIL2(quasar_ind, num_MgII), SNR, map_sigma_MgIIL2(quasar_ind, num_MgII), ...
                    map_N_MgIIL2(quasar_ind, num_MgII), REW_2796(quasar_ind,num_MgII),...
                    REW_2803(quasar_ind,num_MgII),fit_chi2(quasar_ind),log_likelihoods_no_MgII(quasar_ind,1));
                ttl; 
    
                % dz_Doppler = kms_to_z(sqrt(2)*4*map_sigma_MgIIL2(quasar_ind, num_MgII)/1e5); % in km to z
                % dz_mid = (mgii_2803_wavelength - mgii_2796_wavelength)*0.5/mgii_2796_wavelength; % cm/s
                % z_EWhigh = min(map_z_MgIIL2(quasar_ind, num_MgII) + dz_Doppler*(1+z_qso), map_z_MgIIL2(quasar_ind, num_MgII) + dz_mid*(1+z_qso));
                % z_EWlow = map_z_MgIIL2(quasar_ind, num_MgII) - dz_Doppler*(1+z_qso); 
    
                % z_PM_test_plot = z_PM_test(quasar_ind,:);
                % z_PM_test_plot = z_PM_test_plot(z_PM_test_plot>0);
                fid = sprintf('output/dr1/plots/%s/ind-%d-MgII-%s-%d.png',testing_set_name, quasar_ind,this_ID, num_MgII);
                % ind_zoomL2 = (abs(this_z_2796-map_z_MgIIL2(quasar_ind, num_MgII))<20*kms_to_z(map_sigma_MgIIL2(quasar_ind, num_MgII)/1e5)*(1+z_qso));
                % ind_zoomL1 = (abs(this_z_2796-map_z_MgIIL1(quasar_ind, num_MgII))<20*kms_to_z(map_sigma_MgIIL2(quasar_ind, num_MgII)/1e5)*(1+z_qso));
                
                pltQSO(this_flux, this_wavelengths, mu_post, MgII_muL2, MgII_muL1,...
                    this_noise_variance, ttl, fid,ind_not_remove, lambda_cut_range,normalization_cut_range)
                % figure(num_MgII)
                % clf("reset")
                % p = plot((this_wavelengths(ind_not_remove)/2.7964e+03) -1,MgII_muL2(ind_not_remove));
                % end
    
                % fid = sprintf('muPlot/mu-id-%s.png', this_ID);
                % ttl = sprintf('ID:%s, zQSO:%.2f',  this_ID, z_qso);
                % plt_mu(this_flux, this_wavelengths, this_mu, z_qso, this_M, ttl, fid)
                ttl;
            end
        end



    end

    fprintf(' took %0.3fs.\n', toc);

end
% compute model posteriors in numerically safe manner



% save results
variables_to_save = {'releaseTest', 'training_set_name', ...
    'prior_ind', 'releasePrior', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_MgIIs', 'max_z_MgIIs', ...
    'log_priors_no_MgII', 'log_priors_MgII', ...
    'log_likelihoods_no_MgII',  ...
    'sample_log_likelihoods_MgIIL2', 'log_likelihoods_MgIIL2'...
    'log_posteriors_no_MgII', 'log_posteriors_MgIIL1', 'log_posteriors_MgIIL2',...
    'model_posteriors', 'p_no_MgII', 'p_MgIIL1' ...
    'map_z_MgIIL2', 'map_N_MgIIL2', 'map_sigma_MgIIL2' ,'p_MgII', 'REW_2796',...
    'REW_2803','all_SNR','all_B','fit_chi2'};

filename = sprintf('%s/processed_sigma_125_Spline_2200-2300_%s.mat', ...
    processed_directory(releaseTest), ...
    testing_set_name);

save(filename, variables_to_save{:}, '-v7.3');

toc