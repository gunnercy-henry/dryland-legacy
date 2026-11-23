%% 01_quantify_legacy_effects.m
% Reproducible analysis for quantifying productivity legacy effects after
% hydroclimatic extremes.
%
% REQUIRED INPUT FILE (edit input_file below):
%   NDVI : rows x cols x years annual productivity (NDVI)
%   PPT      : rows x cols x years annual water availability
%   AI             : rows x cols aridity index
%
% OPTIONAL INPUT VARIABLES:
%   burned_fraction : rows x cols x years OR rows x cols, fraction [0,1]
%                     used only for fire-exclusion sensitivity
%
% OUTPUTS:
%   legacy_results.mat containing:
%     event_class_all          : rows x cols x years x threshold
%                                +1 sustained, 0 non-detectable,
%                                -1 reversed, NaN no event
%     legacy_event_all         : observed aligned event-level lambda_d
%     legacy_mean_map          : mean observed lambda_d across all real events
%     pixel_class_maps         : rows x cols x threshold
%     freq_sustained           : pixel-level event frequency
%     freq_reversed
%     freq_nondetectable
%     detectable_freq_sustained
%     detectable_freq_reversed
%     event_count_map
%     pixel_null_lower_maps / pixel_null_upper_maps
%     sensitivity structures
%
% MAIN ANALYSIS:
%   Real isolated event: |W_z(t)| > 1.5 and both adjacent years < 1.5
%   Null pool: normal -> normal, both years |W_z| < 1
%   Detectability: Q10/Q90
%   Pixel null bootstrap: 1000 iterations
%
% NOTE:
%   Statistical classification is performed at the original grid resolution.
%   No plotting or display-only spatial dilation is included here.

clear; clc;
rng(1,'twister');

%% ------------------------------------------------------------------------
% 1. USER SETTINGS
% -------------------------------------------------------------------------
input_file  = 'input_legacy_data.mat';
output_file = 'legacy_results.mat';

extreme_threshold      = 1.5;
event_normal_threshold = 1.5;  % fixed for real isolated events
main_null_threshold    = 1.0;
min_valid_years        = 10;
min_control            = 5;
min_events_pixel       = 1;
nboot                  = 1000;

% Main threshold is Q10/Q90. Others are retained for sensitivity analysis.
lower_percentiles = [10 25 5 2.5];
upper_percentiles = [90 75 95 97.5];
threshold_names   = {'10-90','25-75','5-95','2.5-97.5'};
main_q = 1;

% Null-definition sensitivity. Main analysis uses |W_z| < 1.
null_thresholds = [1.0 0.5];
null_names      = {'normal_lt_1','normal_lt_0p5'};
main_null_k = 1;

%% ------------------------------------------------------------------------
% 2. LOAD AND CHECK INPUTS
% -------------------------------------------------------------------------
S = load(input_file);
required = {'NDVI','PPT','AI'};
for k = 1:numel(required)
    if ~isfield(S,required{k})
        error('Missing required variable "%s" in %s.',required{k},input_file);
    end
end

P_raw = double(S.NDVI);
W_raw = double(S.PPT);
AI    = double(S.AI);

if ~isequal(size(P_raw),size(W_raw))
    error('NDVI and PPT must have identical dimensions.');
end
[rows,cols,nYears] = size(P_raw);
if ~isequal(size(AI),[rows,cols])
    error('AI must be rows x cols and match the first two data dimensions.');
end

has_fire = isfield(S,'burned_fraction');
if has_fire
    burned_fraction = double(S.burned_fraction);
end

%% ------------------------------------------------------------------------
% 3. PREALLOCATE CORE TIME-SERIES PRODUCTS
% -------------------------------------------------------------------------
P_z          = nan(rows,cols,nYears,'single');
W_z          = nan(rows,cols,nYears,'single');
P_pred       = nan(rows,cols,nYears,'single');
legacy_raw   = nan(rows,cols,nYears,'single'); % lambda(t)=Pobs-Ppred
water_pvalue = nan(rows,cols,'single');
water_beta   = nan(rows,cols,'single');

% Study-domain mask: drylands + valid productivity.
ndvi_mask = ~isnan(mean(P_raw,3,'omitnan'));
dryland_mask = AI < 0.65 & ndvi_mask;

%% ------------------------------------------------------------------------
% 4. DETREND, STANDARDIZE, AND FIT CONTEMPORANEOUS P ~ W
% -------------------------------------------------------------------------
fprintf('Preprocessing and fitting contemporaneous P ~ W ...\n');
for i = 1:rows
    if mod(i,20)==0 || i==1 || i==rows
        fprintf('  row %d / %d\n',i,rows);
    end
    for j = 1:cols
        if ~dryland_mask(i,j), continue; end

        p = squeeze(P_raw(i,j,:));
        w = squeeze(W_raw(i,j,:));
        valid = isfinite(p) & isfinite(w);
        if sum(valid) < min_valid_years, continue; end

        p_dt = local_detrend_nan(p);
        w_dt = local_detrend_nan(w);
        sp = std(p_dt,'omitnan');
        sw = std(w_dt,'omitnan');
        if ~isfinite(sp) || ~isfinite(sw) || sp==0 || sw==0, continue; end

        pz = (p_dt - mean(p_dt,'omitnan')) ./ sp;
        wz = (w_dt - mean(w_dt,'omitnan')) ./ sw;
        valid2 = isfinite(pz) & isfinite(wz);
        if sum(valid2) < min_valid_years, continue; end

        mdl = fitlm(wz(valid2),pz(valid2));
        pred = nan(nYears,1);
        pred(valid2) = predict(mdl,wz(valid2));
        res = pz - pred;

        P_z(i,j,:)        = single(pz);
        W_z(i,j,:)        = single(wz);
        P_pred(i,j,:)     = single(pred);
        legacy_raw(i,j,:) = single(res);
        water_beta(i,j)   = single(mdl.Coefficients.Estimate(2));
        water_pvalue(i,j) = single(mdl.Coefficients.pValue(2));
    end
end

% Main manuscript domain additionally requires significant water sensitivity.
analysis_mask = dryland_mask & isfinite(water_pvalue) & water_pvalue < 0.05;

%% ------------------------------------------------------------------------
% 5. IDENTIFY REAL ISOLATED EXTREME EVENTS AND OBSERVED lambda_d
% -------------------------------------------------------------------------
legacy_event_all = nan(rows,cols,nYears,'single');
event_sign_all   = nan(rows,cols,nYears,'single');
event_count_map  = zeros(rows,cols,'uint16');
legacy_mean_map  = nan(rows,cols,'single');

fprintf('Identifying isolated extreme events ...\n');
for i = 1:rows
    for j = 1:cols
        if ~analysis_mask(i,j), continue; end
        wz  = squeeze(double(W_z(i,j,:)));
        lam = squeeze(double(legacy_raw(i,j,:)));
        vals = [];

        for t = 2:(nYears-1)
            if any(~isfinite(wz([t-1 t t+1]))) || ~isfinite(lam(t+1)), continue; end
            is_extreme = abs(wz(t)) > extreme_threshold;
            isolated_before = abs(wz(t-1)) < event_normal_threshold;
            isolated_after  = abs(wz(t+1)) < event_normal_threshold;
            if ~(is_extreme && isolated_before && isolated_after), continue; end

            sgn = sign(wz(t));
            lambda_d = lam(t+1) * sgn;
            legacy_event_all(i,j,t) = single(lambda_d);
            event_sign_all(i,j,t)   = single(sgn);
            vals(end+1,1) = lambda_d; 
        end

        if numel(vals) >= min_events_pixel
            event_count_map(i,j) = uint16(numel(vals));
            legacy_mean_map(i,j) = single(mean(vals,'omitnan'));
        end
    end
end

%% ------------------------------------------------------------------------
% 6. EVENT-LEVEL EMPIRICAL NULL CLASSIFICATION
% -------------------------------------------------------------------------
nThreshold = numel(lower_percentiles);
event_class_all = nan(rows,cols,nYears,nThreshold,'single');
% +1 sustained, 0 non-detectable, -1 reversed
control_count_map = zeros(rows,cols,'uint16');

fprintf('Classifying individual events using matched empirical nulls ...\n');
for i = 1:rows
    for j = 1:cols
        if ~analysis_mask(i,j), continue; end
        wz  = squeeze(double(W_z(i,j,:)));
        lam = squeeze(double(legacy_raw(i,j,:)));

        control_res = local_control_pool(wz,lam,main_null_threshold);
        control_count_map(i,j) = uint16(numel(control_res));
        if numel(control_res) < min_control, continue; end

        event_years = find(isfinite(squeeze(legacy_event_all(i,j,:))));
        for ii = 1:numel(event_years)
            t = event_years(ii);
            sgn = double(event_sign_all(i,j,t));
            obs = double(legacy_event_all(i,j,t));
            signed_null = control_res .* sgn;

            for q = 1:nThreshold
                lo = prctile(signed_null,lower_percentiles(q));
                hi = prctile(signed_null,upper_percentiles(q));
                if obs < lo
                    cls = -1;
                elseif obs > hi
                    cls = 1;
                else
                    cls = 0;
                end
                event_class_all(i,j,t,q) = single(cls);
            end
        end
    end
end

%% ------------------------------------------------------------------------
% 7. EVENT FREQUENCIES PER PIXEL
% -------------------------------------------------------------------------
freq_sustained    = nan(rows,cols,nThreshold,'single');
freq_reversed     = nan(rows,cols,nThreshold,'single');
freq_nondetectable= nan(rows,cols,nThreshold,'single');
detectable_freq_sustained = nan(rows,cols,nThreshold,'single');
detectable_freq_reversed  = nan(rows,cols,nThreshold,'single');

for q = 1:nThreshold
    for i = 1:rows
        for j = 1:cols
            c = squeeze(double(event_class_all(i,j,:,q)));
            c = c(isfinite(c));
            if isempty(c), continue; end
            n = numel(c);
            freq_sustained(i,j,q)     = single(sum(c==1)/n);
            freq_reversed(i,j,q)      = single(sum(c==-1)/n);
            freq_nondetectable(i,j,q) = single(sum(c==0)/n);
            det = c(c~=0);
            if ~isempty(det)
                detectable_freq_sustained(i,j,q) = single(sum(det==1)/numel(det));
                detectable_freq_reversed(i,j,q)  = single(sum(det==-1)/numel(det));
            end
        end
    end
end

%% ------------------------------------------------------------------------
% 8. PIXEL-LEVEL BOOTSTRAP NULL FOR MEAN lambda_d
% -------------------------------------------------------------------------
nNull = numel(null_thresholds);
pixel_null_lower_maps = nan(rows,cols,nThreshold,nNull,'single');
pixel_null_upper_maps = nan(rows,cols,nThreshold,nNull,'single');
pixel_class_maps      = nan(rows,cols,nThreshold,nNull,'single');
% +1 sustained, 0 non-detectable, -1 reversed

fprintf('Building pixel-level bootstrap null distributions ...\n');
for k = 1:nNull
    null_thr = null_thresholds(k);
    fprintf('  null definition |W_z| < %.2f\n',null_thr);
    for i = 1:rows
        for j = 1:cols
            if ~analysis_mask(i,j) || ~isfinite(legacy_mean_map(i,j)), continue; end

            wz  = squeeze(double(W_z(i,j,:)));
            lam = squeeze(double(legacy_raw(i,j,:)));
            control_res = local_control_pool(wz,lam,null_thr);
            if numel(control_res) < min_control, continue; end

            signs = squeeze(double(event_sign_all(i,j,:)));
            signs = signs(isfinite(signs));
            nEvent = numel(signs);
            if nEvent < min_events_pixel, continue; end

            idx = randi(numel(control_res),nEvent,nboot);
            sampled = control_res(idx) .* signs;
            null_means = mean(sampled,1,'omitnan');

            obs = double(legacy_mean_map(i,j));
            for q = 1:nThreshold
                lo = prctile(null_means,lower_percentiles(q));
                hi = prctile(null_means,upper_percentiles(q));
                pixel_null_lower_maps(i,j,q,k) = single(lo);
                pixel_null_upper_maps(i,j,q,k) = single(hi);
                if obs < lo
                    cls = -1;
                elseif obs > hi
                    cls = 1;
                else
                    cls = 0;
                end
                pixel_class_maps(i,j,q,k) = single(cls);
            end
        end
    end
end

%% ------------------------------------------------------------------------
% 9. OPTIONAL FIRE-EXCLUSION SENSITIVITY
% -------------------------------------------------------------------------
fire_sensitivity = struct();
if has_fire
    if ndims(burned_fraction)==3
        fire_affected = max(burned_fraction,[],3,'omitnan') > 0.50;
    elseif isequal(size(burned_fraction),[rows,cols])
        fire_affected = burned_fraction > 0.50;
    else
        warning('burned_fraction has incompatible dimensions; fire sensitivity skipped.');
        fire_affected = false(rows,cols);
    end
    fire_sensitivity.fire_affected = fire_affected;
    fire_sensitivity.legacy_mean_with_fire = legacy_mean_map;
    tmp = legacy_mean_map;
    tmp(fire_affected) = NaN;
    fire_sensitivity.legacy_mean_without_fire = tmp;
end

%% ------------------------------------------------------------------------
% 10. ALTERNATIVE BASELINE SENSITIVITY: MEAN P DURING NON-EXTREME YEARS
% -------------------------------------------------------------------------
legacy_mean_altbaseline = nan(rows,cols,'single');
for i = 1:rows
    for j = 1:cols
        if ~analysis_mask(i,j), continue; end
        pz = squeeze(double(P_z(i,j,:)));
        wz = squeeze(double(W_z(i,j,:)));
        if all(~isfinite(pz)) || all(~isfinite(wz)), continue; end
        nonextreme = isfinite(pz) & isfinite(wz) & abs(wz) < extreme_threshold;
        if ~any(nonextreme), continue; end
        baseline = mean(pz(nonextreme),'omitnan');
        vals = [];
        for t = 2:(nYears-1)
            if any(~isfinite(wz([t-1 t t+1]))) || ~isfinite(pz(t+1)), continue; end
            if abs(wz(t)) > extreme_threshold && ...
               abs(wz(t-1)) < event_normal_threshold && ...
               abs(wz(t+1)) < event_normal_threshold
                vals(end+1,1) = (pz(t+1)-baseline) * sign(wz(t)); 
            end
        end
        if ~isempty(vals)
            legacy_mean_altbaseline(i,j) = single(mean(vals,'omitnan'));
        end
    end
end

%% ------------------------------------------------------------------------
% 11. SUMMARY FOR MAIN ANALYSIS
% -------------------------------------------------------------------------
main_event_class = event_class_all(:,:,:,main_q);
vec = double(main_event_class(:));
vec = vec(isfinite(vec));
fprintf('\nMAIN EVENT-LEVEL CLASSIFICATION (Q10/Q90)\n');
fprintf('  N events        = %d\n',numel(vec));
fprintf('  Sustained       = %.2f%%\n',100*sum(vec==1)/numel(vec));
fprintf('  Reversed        = %.2f%%\n',100*sum(vec==-1)/numel(vec));
fprintf('  Non-detectable  = %.2f%%\n',100*sum(vec==0)/numel(vec));

det = vec(vec~=0);
if ~isempty(det)
    fprintf('  Among detectable: sustained %.2f%%, reversed %.2f%%\n', ...
        100*sum(det==1)/numel(det),100*sum(det==-1)/numel(det));
end

main_pixel_class = pixel_class_maps(:,:,main_q,main_null_k);
pv = double(main_pixel_class(:));
pv = pv(isfinite(pv));
fprintf('\nMAIN PIXEL-LEVEL CLASSIFICATION (null |W_z|<1; Q10/Q90)\n');
fprintf('  N pixels        = %d\n',numel(pv));
fprintf('  Sustained       = %.2f%%\n',100*sum(pv==1)/numel(pv));
fprintf('  Reversed        = %.2f%%\n',100*sum(pv==-1)/numel(pv));
fprintf('  Non-detectable  = %.2f%%\n',100*sum(pv==0)/numel(pv));

%% ------------------------------------------------------------------------
% 12. SAVE
% -------------------------------------------------------------------------
save(output_file, ...
    'P_z','W_z','P_pred','legacy_raw','water_beta','water_pvalue', ...
    'dryland_mask','analysis_mask', ...
    'legacy_event_all','event_sign_all','event_count_map','legacy_mean_map', ...
    'event_class_all','control_count_map', ...
    'freq_sustained','freq_reversed','freq_nondetectable', ...
    'detectable_freq_sustained','detectable_freq_reversed', ...
    'pixel_null_lower_maps','pixel_null_upper_maps','pixel_class_maps', ...
    'legacy_mean_altbaseline','fire_sensitivity', ...
    'extreme_threshold','event_normal_threshold','main_null_threshold', ...
    'lower_percentiles','upper_percentiles','threshold_names', ...
    'null_thresholds','null_names','main_q','main_null_k','nboot', ...
    '-v7.3');

fprintf('\nSaved: %s\n',output_file);

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function x_dt = local_detrend_nan(x)
    x = double(x(:));
    tt = (1:numel(x))';
    good = isfinite(x);
    x_dt = nan(size(x));
    if sum(good) < 2, return; end
    b = polyfit(tt(good),x(good),1);
    trend = polyval(b,tt);
    x_dt(good) = x(good) - trend(good);
end

function control_res = local_control_pool(wz,lambda_raw,null_thr)
    wz = double(wz(:));
    lambda_raw = double(lambda_raw(:));
    nYears = numel(wz);
    control_res = [];
    for tc = 1:(nYears-1)
        if ~isfinite(wz(tc)) || ~isfinite(wz(tc+1)) || ~isfinite(lambda_raw(tc+1))
            continue;
        end
        if abs(wz(tc)) < null_thr && abs(wz(tc+1)) < null_thr
            control_res(end+1,1) = lambda_raw(tc+1); 
        end
    end
end
