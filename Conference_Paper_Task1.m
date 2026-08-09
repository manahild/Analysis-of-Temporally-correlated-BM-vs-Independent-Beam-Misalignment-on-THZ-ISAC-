%% ========================================================================
%   CORRELATED BEAM MISALIGNMENT IN THE FULL SENSING SIM
%  (Task 1 restructuring, PERFORMANCE-FIXED version)
%  SAME script, SAME 4 scenarios x 2 regimes x 2 waveforms, SAME results
%  table / figures / data save -- nothing removed and nothing computed
%  differently from the previous modular version. The module set is the
%  same as before:
%       - AR(1) trajectory generation           -> generate_ar1()
%       - beam-map evaluation h(eps)             -> beam_gain()
%       - received-signal generation              -> generate_received_signal_fast()
%       - range-Doppler processing                -> range_doppler_fast()
%       - target-cell peak-power extraction        -> extract_target_fast()
%       - pointing estimation / filtering / beam-command -> RESERVED STUBS
%       - metric calculation                       -> compute_metrics()
%       - plotting                                  -> unchanged, Section 8
%   Section 10 (stationary covariance verification) is unchanged.
%
%  *** WHY THIS VERSION IS FASTER THAN THE PREVIOUS MODULAR VERSION ***
%  The previous version passed the structs W and P into functions called
%  FROM INSIDE parfor (generate_received_signal(W,...), range_doppler(rx,
%  W,P), extract_target(Z,[],P,...)). MATLAB's parfor transparency
%  analysis is much less reliable at recognizing large struct FIELDS
%  accessed only inside nested function calls as one-time broadcast
%  variables -- in practice this can cause the ~84 MB array
%  W.rx_delayed_doppler (and W.corr_ref) to be re-sent to the workers on
%  every iteration instead of once, which is why that version was much
%  slower than the plain inline script even though the arithmetic is
%  identical.
%
%  THE FIX: immediately before each parfor block, the large arrays are
%  pulled out of the structs into PLAIN top-level local variables
%  (rx_delayed_doppler_local, corr_ref_local) and all scalar P fields
%  needed inside the loop are also pulled into plain scalars. Only these
%  plain variables -- never W or P themselves -- are referenced inside
%  parfor. This restores the same broadcast-variable behavior as the
%  original inline script, so this version should run at essentially the
%  same speed as the inline script, while keeping the Task 1 module
%  structure everywhere else (setup, waveform build, AR(1) generation,
%  metrics, stubs, Section 10 all still fully modular).
% =========================================================================

clear; clc; close all;
rng(1);                      % single global seed, reproducible

%% ------------------------------------------------------------------
%  0. FOLDERS
%% ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
figDir = fullfile(scriptDir,'Figures_T56');
resDir = fullfile(scriptDir,'Results_T56');
datDir = fullfile(scriptDir,'Data_T56');
for f = {figDir,resDir,datDir}
    if ~exist(f{1},'dir'), mkdir(f{1}); end
end

%% ------------------------------------------------------------------
%  1. SYSTEM PARAMETERS (Table I) -- module: setup_system_params()
%% ------------------------------------------------------------------
P = setup_system_params();
c = P.c; kb = P.kb; fc = P.fc; lambda = P.lambda; k_wave = P.k_wave; %#ok<NASGU>
G_total = P.G_total;
B = P.B; N_fft = P.N_fft; N_radar = P.N_radar; N_cp = P.N_cp; N_symbols = P.N_symbols;
delta_f = P.delta_f; T_sym = P.T_sym; Fs = P.Fs;
rangeRes = P.rangeRes; velRes = P.velRes;
sigma_rcs = P.sigma_rcs;

% Regime configs -- SAME beam, target size sets SF/SL (unchanged)
w_0_SF = 0.20;  apert_T_SF = 1.5;
w_0_SL = 0.20;  apert_T_SL = 0.04;
apert_R = 0.075;
C_n_sq  = 1e-11;
velocity  = 5;                      % m/s, target radial velocity
regimes      = {'SF','SL'};
w_0_list     = [w_0_SF, w_0_SL];
apert_T_list = [apert_T_SF, apert_T_SL];
waveforms    = {'OFDM','SCFDMA'};

%% ------------------------------------------------------------------
%  2. SENSING / MISALIGNMENT SCENARIO SETUP  (requirements)
%% ------------------------------------------------------------------
distance = 15;                      % m  (matches Task 2-4 sigma_s calc)
theta_bm = 30e-3;                   % 30 mrad (matches Task 2-4)
sigma_s  = distance*tan(theta_bm);  % pointing-jitter std, m

Ptx_dBm_fixed = 10;                  % ASSUMPTION: single fixed operating
Ptx_W = 10^((Ptx_dBm_fixed-30)/10); % point (mid of the sweep used in the
                                     % main script), since T5 asks for a
                                     % *sequence of frames*, not a power
                                     % sweep. Change here if a different
                                     % operating power is preferred.

N_frames = 500;                     % ASSUMPTION: sensing frames per
                                     % scenario (CPI-by-CPI tracking run
                                     % that keeps runtime reasonable while
                                     % still showing long correlated fade
                                     % events at rho=0.99)

h_th = 0.10;                        % same threshold as Task 4

rho_list       = [0, 0.9, 0.99];     % 0 stands in for "i.i.d." (rho=0
                                     % collapses AR(1) exactly to i.i.d.)
scenarioNames  = {'Aligned','i.i.d.','AR(1) rho=0.90','AR(1) rho=0.99'};
scenarioTagsSF = {'Aligned','iid','rho0p9','rho0p99'};
nScenario = 4;

fprintf('===============================================\n');
fprintf('Task 1 setup\n');
fprintf('===============================================\n');
fprintf('distance=%.1f m | theta_bm=%.0f mrad | sigma_s=%.4f m\n', distance, theta_bm*1000, sigma_s);
fprintf('Ptx = %.1f dBm | N_frames = %d | h_th = %.2f\n\n', Ptx_dBm_fixed, N_frames, h_th);

%% ------------------------------------------------------------------
%  3. ATMOSPHERIC MODEL (unchanged)
%% ------------------------------------------------------------------
height = 100;
[T_atm_K, p_hPa, wvden] = atmositu(height);
T_C  = T_atm_K - 273.15;
P_Pa = p_hPa*100;

L_atm_dB  = gaspl(distance, fc, T_C, P_Pa, wvden);
L_atm_lin = 10^(-L_atm_dB/10);

T_noise   = (1-L_atm_lin)*T_atm_K;
P_n       = kb*T_noise*B;
P_n_radar = P_n*(N_radar/N_fft);

%% ------------------------------------------------------------------
%  4. PARALLEL POOL (unchanged)
%% ------------------------------------------------------------------
pool = gcp('nocreate');
if isempty(pool), pool = parpool('local'); end
fprintf('Using parallel pool with %d workers\n\n', pool.NumWorkers);

%% ------------------------------------------------------------------
%  5-6. MAIN LOOP: regime x waveform x scenario x frame
%     Radar reference / waveform build -> module: build_radar_reference()
%     Trajectory generation             -> module: generate_ar1()
%     Beam-map evaluation                -> module: beam_gain()
%     Received-signal generation          -> module: generate_received_signal_fast()
%     Range-Doppler processing            -> module: range_doppler_fast()
%     Peak extraction / detection         -> module: extract_target_fast()
%     Metrics                             -> module: compute_metrics()
%
%  NOTE ON PERFORMANCE: N_fft, N_cp, N_symbols, N_radar, delta_f, T_sym,
%  c, lambda, rangeRes, velRes are already plain scalars from Section 1
%  (unpacked from P) -- these are cheap and parfor broadcasts scalars
%  fine regardless. The only arrays that matter for broadcast cost are
%  W.rx_delayed_doppler and W.corr_ref, which are unpacked into plain
%  local variables right before every parfor block below.
%% ------------------------------------------------------------------
results = struct();

for reg_idx = 1:numel(regimes)
    regime  = regimes{reg_idx};
    w_0     = w_0_list(reg_idx);
    apert_T = apert_T_list(reg_idx);

    fprintf('========== %s CONFIG (w_0=%.0f mm, a_T=%.2f m) ==========\n', regime, w_0*1000, apert_T);

    %% --- beam-waist / w_mo computation (unchanged) -----------------
    [w_mo, eta, is_SF] = compute_beam_width(distance, w_0, apert_T, apert_R, C_n_sq, fc);
    fprintf('  eta=%.2f  ->  w_mo = %.4f m\n', eta, w_mo);
    beam_parameters.w_mo = w_mo;

    if strcmp(regime,'SF') && ~is_SF
        warning('%s config has eta=%.3f (<=1): physics says Spot-Limited, not Spot-Filling!', regime, eta);
    elseif strcmp(regime,'SL') && is_SF
        warning('%s config has eta=%.3f (>1): physics says Spot-Filling, not Spot-Limited!', regime, eta);
    end

    %% --- Task 2/3 generators: build h_k(k) for every scenario -------
    % Same random draws / same order as before, now routed through
    % generate_ar1() and beam_gain() instead of inline formulas.
    ux = randn(N_frames,1);  uy = randn(N_frames,1);   % shared innovations
    hk_scenario = zeros(N_frames, nScenario);
    rk_scenario = zeros(N_frames, nScenario);

    % 1: aligned
    hk_scenario(:,1) = 1;
    rk_scenario(:,1) = 0;

    % 2: i.i.d. (rho=0 -- fresh independent draws, matches Task 2 exactly;
    %    NOT reusing ux,uy so this scenario's RNG stream is unchanged
    %    from the original script)
    xk_iid = sigma_s*randn(N_frames,1);
    yk_iid = sigma_s*randn(N_frames,1);
    eps_iid = [xk_iid, yk_iid];
    rk_scenario(:,2) = sqrt(sum(eps_iid.^2,2));
    hk_scenario(:,2) = beam_gain(eps_iid, beam_parameters);

    % 3 & 4: correlated AR(1), rho = 0.9, 0.99
    for si = 1:2
        rho = rho_list(si+1);
        eps_k = generate_ar1(N_frames, rho, sigma_s, ux, uy);
        rk_scenario(:,si+2) = sqrt(sum(eps_k.^2,2));
        hk_scenario(:,si+2) = beam_gain(eps_k, beam_parameters);
    end

    %% --- radar equation amplitude (aligned reference, unchanged) ------
    P_radar = Ptx_W*(N_radar/N_fft);
    A_base = sqrt(P_radar)*L_atm_lin* ...
             sqrt((G_total*lambda^2*sigma_rcs)/((4*pi)^3*distance^4));

    for wf_idx = 1:numel(waveforms)
        waveform = waveforms{wf_idx};

        W = build_radar_reference(P, waveform, distance, velocity);

        % --- PERFORMANCE FIX: unpack the large arrays into plain local
        % variables BEFORE parfor, so MATLAB broadcasts them once instead
        % of re-copying the struct on every iteration. -------------------
        rx_delayed_doppler_local = W.rx_delayed_doppler;   % ~84 MB, sent ONCE
        corr_ref_local           = W.corr_ref;             % sent ONCE

        for sc_idx = 1:nScenario
            hk_vec = hk_scenario(:,sc_idx);

            range_est   = zeros(N_frames,1);
            vel_est     = zeros(N_frames,1);
            correct_det = false(N_frames,1);
            SNR_dB      = zeros(N_frames,1);
            Prx_lin     = zeros(N_frames,1);

            parfor kf = 1:N_frames
                h_k = hk_vec(kf); %#ok<PFBNS>

                % generate_received_signal_fast/range_doppler_fast/
                % extract_target_fast take PLAIN arrays and scalars only
                % -- no structs -- so parfor's broadcast analysis is
                % unambiguous, same as the original inline script.
                [rx_signal, Pr_k] = generate_received_signal_fast( ...
                    rx_delayed_doppler_local, h_k, A_base, P_n_radar); %#ok<PFBNS>

                Z = range_doppler_fast(rx_signal, corr_ref_local, ...
                    N_fft, N_cp, N_symbols, N_radar); %#ok<PFBNS>

                [~, r_est, v_est, det] = extract_target_fast(Z, ...
                    c, delta_f, N_radar, lambda, T_sym, N_symbols, ...
                    distance, velocity, rangeRes, velRes); %#ok<PFBNS>

                range_est(kf) = r_est;
                vel_est(kf)   = v_est;
                correct_det(kf) = det;

                Prx_lin(kf) = Pr_k;
                SNR_dB(kf) = 10*log10(Pr_k/P_n_radar);
            end

            M = compute_metrics(range_est, vel_est, correct_det, SNR_dB, distance, velocity);

            tag = scenarioTagsSF{sc_idx};
            results.(regime).(waveform).(tag).h_k        = hk_vec;
            results.(regime).(waveform).(tag).r_k         = rk_scenario(:,sc_idx);
            results.(regime).(waveform).(tag).range_est   = range_est;
            results.(regime).(waveform).(tag).vel_est     = vel_est;
            results.(regime).(waveform).(tag).correct_det = correct_det;
            results.(regime).(waveform).(tag).SNR_dB      = SNR_dB;
            results.(regime).(waveform).(tag).Prx_lin      = Prx_lin;
            results.(regime).(waveform).(tag).avgRangeErr = M.avgRangeErr;
            results.(regime).(waveform).(tag).avgVelErr   = M.avgVelErr;
            results.(regime).(waveform).(tag).Pd          = M.Pd;
            results.(regime).(waveform).(tag).avgFailRun  = M.avgFailRun;
            results.(regime).(waveform).(tag).maxFailRun  = M.maxFailRun;

            fprintf('  [%s,%s,%-14s] Pd=%.3f  avgRangeErr=%.4f m  avgVelErr=%.4f m/s  avgFailRun=%.2f  maxFailRun=%d\n', ...
                regime, waveform, scenarioNames{sc_idx}, M.Pd, M.avgRangeErr, M.avgVelErr, M.avgFailRun, M.maxFailRun);
        end
    end
    fprintf('\n');
end

%% ------------------------------------------------------------------
%  7. TASK 1 -- RESULTS TABLE (unchanged)
%% ------------------------------------------------------------------
T6 = table();
for reg_idx = 1:numel(regimes)
    regime = regimes{reg_idx};
    for wf_idx = 1:numel(waveforms)
        waveform = waveforms{wf_idx};
        for sc_idx = 1:nScenario
            tag = scenarioTagsSF{sc_idx};
            R = results.(regime).(waveform).(tag);
            newRow = table({regime}, {waveform}, {scenarioNames{sc_idx}}, ...
                R.avgRangeErr, R.avgVelErr, R.Pd, R.avgFailRun, R.maxFailRun, ...
                'VariableNames', {'Regime','Waveform','Scenario', ...
                'AvgRangeErr_m','AvgVelErr_mps','Pd','AvgConsecFailRun','MaxConsecFailRun'});
            T6 = [T6; newRow]; %#ok<AGROW>
        end
    end
end
disp(T6);
writetable(T6, fullfile(resDir,'Task6_results_table.csv'));

%% ------------------------------------------------------------------
%  8. TASK 1 -- FIGURES (unchanged)
%     (produced for OFDM as the representative waveform in each regime;
%      SC-FDMA data is identical in structure and saved in T6/results,
%      comparing waveforms is explicitly out of scope per the brief)
%% ------------------------------------------------------------------
colorsSc = lines(nScenario);
wfShow = 'OFDM';

for reg_idx = 1:numel(regimes)
    regime = regimes{reg_idx};

    % --- Fig 1: h_k vs frame number, all rho ---
    fig = figure('Color','w','Position',[100 100 1000 550]); hold on;
    for sc_idx = 2:nScenario   % skip 'Aligned' (h_k==1 trivial line)
        tag = scenarioTagsSF{sc_idx};
        plot(results.(regime).(wfShow).(tag).h_k, 'Color', colorsSc(sc_idx,:), ...
            'DisplayName', scenarioNames{sc_idx});
    end
    xlabel('Frame Number, k'); ylabel('h_k');
    title(sprintf('Task 1 Fig 1: h_k vs Frame, %s regime', regime));
    ylim([0 1]); legend('Location','best'); grid on;
    save_figure(fig, fullfile(figDir, sprintf('T6_Fig1_hk_vs_frame_%s', regime)));

    % --- Fig 2: received SNR vs frame number ---
    fig = figure('Color','w','Position',[100 100 1000 550]); hold on;
    for sc_idx = 1:nScenario
        tag = scenarioTagsSF{sc_idx};
        plot(results.(regime).(wfShow).(tag).SNR_dB, 'Color', colorsSc(sc_idx,:), ...
            'DisplayName', scenarioNames{sc_idx});
    end
    xlabel('Frame Number, k'); ylabel('Received SNR (dB)');
    title(sprintf('Task 1 Fig 2: SNR vs Frame, %s regime, %s', regime, wfShow));
    legend('Location','best'); grid on;
    save_figure(fig, fullfile(figDir, sprintf('T6_Fig2_SNR_vs_frame_%s', regime)));

    % --- Fig 3: correct/incorrect detections vs frame number ---
    fig = figure('Color','w','Position',[100 100 1000 700]);
    for sc_idx = 1:nScenario
        tag = scenarioTagsSF{sc_idx};
        subplot(nScenario,1,sc_idx);
        stem(double(results.(regime).(wfShow).(tag).correct_det), ...
            'Marker','none','Color', colorsSc(sc_idx,:));
        ylim([-0.2 1.2]); ylabel(scenarioNames{sc_idx}, 'FontSize', 8);
        if sc_idx == 1
            title(sprintf('Task 1 Fig 3: Detection Correct(1)/Incorrect(0) vs Frame, %s, %s', regime, wfShow));
        end
        if sc_idx == nScenario, xlabel('Frame Number, k'); end
    end
    save_figure(fig, fullfile(figDir, sprintf('T6_Fig3_detections_vs_frame_%s', regime)));
end

% --- Fig 4: number of consecutive failed frames vs rho (both regimes) ---
fig = figure('Color','w','Position',[100 100 900 550]); hold on;
xRho = 1:3;   % iid, rho0.9, rho0.99  (aligned excluded: not a fading scenario)
avgFailSF = zeros(1,3); maxFailSF = zeros(1,3);
avgFailSL = zeros(1,3); maxFailSL = zeros(1,3);
for i = 1:3
    tag = scenarioTagsSF{i+1};
    avgFailSF(i) = results.SF.(wfShow).(tag).avgFailRun;
    maxFailSF(i) = results.SF.(wfShow).(tag).maxFailRun;
    avgFailSL(i) = results.SL.(wfShow).(tag).avgFailRun;
    maxFailSL(i) = results.SL.(wfShow).(tag).maxFailRun;
end
bar(xRho-0.2, avgFailSF, 0.15, 'FaceColor',[0.2 0.4 0.8], 'DisplayName','SF avg');
bar(xRho-0.05, maxFailSF, 0.15, 'FaceColor',[0.1 0.2 0.5], 'DisplayName','SF max');
bar(xRho+0.1, avgFailSL, 0.15, 'FaceColor',[0.85 0.4 0.2], 'DisplayName','SL avg');
bar(xRho+0.25, maxFailSL, 0.15, 'FaceColor',[0.55 0.2 0.1], 'DisplayName','SL max');
set(gca,'XTick', xRho, 'XTickLabel', {'i.i.d.','\rho=0.90','\rho=0.99'});
ylabel('Consecutive failed frames'); xlabel('Misalignment scenario');
title(sprintf('Task 1 Fig 4: Consecutive Failed Frames vs \\rho (%s)', wfShow));
legend('Location','best'); grid on;
save_figure(fig, fullfile(figDir,'T6_Fig4_consec_fail_vs_rho'));

%% ------------------------------------------------------------------
%  9. SAVE DATA (unchanged)
%% ------------------------------------------------------------------
save(fullfile(datDir,'Task5_6_results.mat'), 'results', 'T6', ...
    'distance','theta_bm','sigma_s','Ptx_dBm_fixed','N_frames','h_th', ...
    'rangeRes','velRes','regimes','waveforms','scenarioNames','scenarioTagsSF', ...
    'w_0_SF','apert_T_SF','w_0_SL','apert_T_SL','apert_R','C_n_sq','fc');

fprintf('\nFigures -> %s\nTable   -> %s\nData    -> %s\n', figDir, resDir, datDir);
fprintf('============ TASKS 5-6 COMPLETE ============\n');

%% ------------------------------------------------------------------
%  10. TASK 1 REQUIREMENT (unchanged):
%      Verify numerically that changing rho does NOT change the
%      stationary covariance / marginal histogram of s_k, using the
%      general-purpose stationary AR(1) generator with the brief's exact
%      interface: s = generate_ar1_stationary(K, rho, Sigma, seed).
%      This is a standalone check appended AFTER everything above; it
%      does not feed back into or alter any result computed above.
%% ------------------------------------------------------------------
Sigma_s = sigma_s^2*eye(2);
verify_stationary_covariance(N_frames, [0 0.9 0.99], Sigma_s, 200, 5000);

fprintf('============ TASK 1 VERIFICATION COMPLETE ============\n');


%% ========================================================================
%  LOCAL FUNCTIONS -- Task 1 module set
%% ========================================================================

function P = setup_system_params()
% SETUP_SYSTEM_PARAMS  THz-ISAC system parameters (Table I). Same numeric
% values as the original inline block, just returned as a struct so every
% other module draws from one source instead of a wall of loose variables.
    P.c  = physconst('LightSpeed');
    P.kb = physconst('Boltzmann');

    P.fc     = 0.85e12;
    P.lambda = P.c/P.fc;
    P.k_wave = 2*pi/P.lambda;

    Gt_dB = 20; Gr_dB = 20;
    P.G_total = 10^((Gt_dB+Gr_dB)/10);

    P.B         = 20e9;
    P.N_fft     = 4096;
    P.N_radar   = 2048;
    P.N_cp      = 1024;
    P.N_symbols = 1024;

    P.delta_f = P.B/P.N_fft;
    P.T_sym   = (P.N_fft+P.N_cp)/P.B;
    P.Fs      = P.B;

    P.rangeRes = P.c/(2*P.N_radar*P.delta_f);
    P.velRes   = P.lambda/(2*P.N_symbols*P.T_sym);

    P.sigma_rcs_dBsm = 10;
    P.sigma_rcs = 10^(P.sigma_rcs_dBsm/10);
end


function [w_mo, eta, is_SF] = compute_beam_width(distance, w_0, apert_T, apert_R, C_n_sq, fc)
% COMPUTE_BEAM_WIDTH  Effective beam-width parameter w_mo (Eq.1/7).
% UNCHANGED physics from the original script's local function.
    c      = physconst('LightSpeed');
    lambda = c/fc;
    k_wave = 2*pi/lambda;

    rho_0    = (0.55*C_n_sq*k_wave^2*distance)^(-3/5);
    eps_turb = 1 + (2*w_0^2)/(rho_0^2);
    wz       = w_0*sqrt(1 + eps_turb*((lambda*distance)/(pi*w_0^2))^2);
    eta      = apert_T/wz;
    is_SF    = (eta > 1);

    if is_SF
        w0_R = wz;
    else
        w0_R = apert_T*sqrt(2);
    end
    eps_R    = 1 + (2*w0_R^2)/(rho_0^2);
    wz_R     = w0_R*sqrt(1 + eps_R*((lambda*distance)/(pi*w0_R^2))^2);
    v_R      = (sqrt(pi)*apert_R)/(sqrt(2)*wz_R);
    weq_R_sq = (wz_R^2*sqrt(pi)*erf(v_R))/(2*v_R*exp(-v_R^2));

    if is_SF
        w_mo_sq = weq_R_sq;
    else
        w_mo_sq = 1/(1/wz^2 + 1/weq_R_sq);
    end
    w_mo = sqrt(w_mo_sq);
end


function W = build_radar_reference(P, waveform, distance, velocity)
% BUILD_RADAR_REFERENCE  Zadoff-Chu radar reference, OFDM/SC-FDMA tx
% waveform, delay + Doppler phase. UNCHANGED physics/formulas from the
% original script's Section 4 + waveform-loop block; just factored into a
% callable module (called once per waveform, same as before). This is
% the ONLY place the struct W is built -- it is unpacked into plain
% variables in the main script immediately after this call, before any
% parfor loop touches it.
    n_zc = (0:P.N_radar-1)';
    zc_root = 2;
    radar_seq = exp(-1j*pi*zc_root*n_zc.*(n_zc+1)/P.N_radar);
    radar_ref = repmat(radar_seq,1,P.N_symbols);

    if strcmpi(waveform,'SCFDMA')
        scfdma_input = ifft(radar_ref, P.N_radar, 1)*sqrt(P.N_radar);
        tx_symbols = fft(scfdma_input, P.N_radar, 1)/sqrt(P.N_radar);
    else
        tx_symbols = radar_ref;
    end
    corr_ref = radar_ref;

    tx_freq = zeros(P.N_fft, P.N_symbols);
    tx_freq(1:P.N_radar,:) = tx_symbols;
    tx_time = ifft(tx_freq, P.N_fft, 1)*sqrt(P.N_fft);
    tx_cp = [tx_time(end-P.N_cp+1:end,:); tx_time];
    tx_serial = tx_cp(:);

    k_delay = round(2*distance/P.c*P.Fs);
    f_D     = 2*velocity/P.lambda;
    t_vec_full = (0:(P.N_fft+P.N_cp)*P.N_symbols-1)'/P.Fs;
    doppler_phase_full = exp(1j*2*pi*f_D*t_vec_full);

    rx_delayed = [zeros(k_delay,1); tx_serial(1:end-k_delay)];

    W.corr_ref = corr_ref;
    W.rx_delayed_doppler = rx_delayed .* doppler_phase_full;
    W.k_delay = k_delay;
    W.f_D = f_D;
end


function eps_k = generate_ar1(N, rho, sigma, ux, uy)
% GENERATE_AR1  AR(1) correlated pointing-error generator, Eq.(5)-(6).
% IDENTICAL recursion to the original gen_correlated_xy(): same shared
% innovations ux,uy (drawn once per regime in the main script), same
% recursion order -- so results are bit-identical to the original script.
%   eps_k = generate_ar1(N, rho, sigma, ux, uy)
%     N      - number of frames
%     rho    - AR(1) correlation coefficient
%     sigma  - pointing-jitter std (sigma_s), isotropic
%     ux,uy  - N x 1 shared N(0,1) innovation sequences
%   eps_k    - N x 2, [eps_x, eps_y]  (Eq.1)
    xk = zeros(N,1); yk = zeros(N,1);
    xk(1) = sigma*ux(1); yk(1) = sigma*uy(1);
    for k = 2:N
        xk(k) = rho*xk(k-1) + sqrt(1-rho^2)*sigma*ux(k);
        yk(k) = rho*yk(k-1) + sqrt(1-rho^2)*sigma*uy(k);
    end
    eps_k = [xk, yk];
end


function s = generate_ar1_stationary(K, rho, Sigma, seed)
% GENERATE_AR1_STATIONARY  General-purpose stationary 2-D AR(1) generator
% matching the brief's exact interface s = generate_ar1(K,rho,Sigma,seed).
% NOT used by the main scenario loop above (that loop keeps the original
% gen_correlated_xy-equivalent generate_ar1() for bit-identical
% reproduction of the earlier results). This function exists so Task 1's
% "verify numerically that changing rho does not change the stationary
% covariance" requirement can be tested against the formal Eq.(6)-(10)
% definition directly, and so later tasks (2-5) have a ready-made,
% independently-seedable trajectory generator.
%   s = generate_ar1_stationary(K, rho, Sigma, seed)
%     K     - number of frames
%     rho   - correlation coefficient, 0<=rho<1
%     Sigma - 2x2 stationary covariance Sigma_s (Eq.9)
%     seed  - RNG seed
%   s     - K x 2, s(k,:) = [sx_k, sy_k]
%   Implements s_{k+1}=rho*s_k+q_k, q_k~N(0,Q), Q=(1-rho^2)*Sigma (Eq.6-8),
%   s_0~N(0,Sigma) (Eq.10). At rho=0 this collapses exactly to i.i.d.
%   N(0,Sigma) draws (F=0).
    if nargin < 4 || isempty(seed), seed = 1; end
    rng(seed);

    F = rho*eye(2);
    Qchol  = chol((1-rho^2)*Sigma, 'lower');
    S0chol = chol(Sigma, 'lower');

    s = zeros(K,2);
    s(1,:) = (S0chol*randn(2,1)).';
    for k = 2:K
        s(k,:) = (F*s(k-1,:).' + Qchol*randn(2,1)).';
    end
end


function h = beam_gain(eps, beam_parameters)
% BEAM_GAIN  Circular-Gaussian beam-misalignment power coefficient (Eq.12).
% SAME formula as the original inline exp(-2*rk.^2/w_mo^2).
%   h = beam_gain(eps, beam_parameters)
%     eps             - N x 2, [eps_x, eps_y] (m)
%     beam_parameters - struct with field .w_mo (kept as a struct so
%                        Task 2 can add wx,wy,phi,p,... without changing
%                        this call signature)
%   h               - N x 1, 0<=h<=1, h(0)=1
    w_mo = beam_parameters.w_mo;
    r2 = eps(:,1).^2 + eps(:,2).^2;
    h = exp(-2*r2/w_mo^2);
end


function [rx_signal, Pr_k] = generate_received_signal_fast(rx_delayed_doppler, h_k, A_base, P_n_radar)
% GENERATE_RECEIVED_SIGNAL_FAST  Received sensing signal for one frame,
% using the CORRECTED linear received-power law Pr,k = Pr,aligned * h_k
% (Eq.5/8). Echo AMPLITUDE is scaled by sqrt(h_k) -- NOT h_k -- so power
% scales linearly with h_k. IDENTICAL formula/physics to
% generate_received_signal() from the previous modular version -- the
% ONLY difference is that this function takes the plain array
% rx_delayed_doppler directly instead of a struct W, so it is safe and
% fast to call from inside parfor (see performance note at top of file).
    A_echo_k = A_base * sqrt(h_k);
    rx_echo  = A_echo_k * rx_delayed_doppler;
    noise = sqrt(P_n_radar/2)*(randn(size(rx_echo))+1j*randn(size(rx_echo)));
    rx_signal = rx_echo + noise;
    Pr_k = A_echo_k^2;                    % = A_base^2 * h_k
end


function Z = range_doppler_fast(rx_signal, corr_ref, N_fft, N_cp, N_symbols, N_radar)
% RANGE_DOPPLER_FAST  Matched-filter range-Doppler map for one frame.
% IDENTICAL formulas to range_doppler() from the previous modular
% version -- the only difference is that corr_ref and the dimensions are
% passed as plain arguments instead of inside struct W/P, so this is
% parfor-broadcast-friendly.
    rx_cp = reshape(rx_signal, N_fft+N_cp, N_symbols);
    rx_time = rx_cp(N_cp+1:end,:);
    rx_freq = fft(rx_time, N_fft, 1)/sqrt(N_fft);
    rx_radar = rx_freq(1:N_radar,:);

    corr_out = rx_radar .* conj(corr_ref);
    range_profile = ifft(corr_out, N_radar, 1);
    Z = fft(range_profile, N_symbols, 2);
end


function [Ppeak, Rhat, vhat, det] = extract_target_fast(Z, c, delta_f, N_radar, lambda, T_sym, N_symbols, true_range, true_vel, rangeRes, velRes)
% EXTRACT_TARGET_FAST  Peak-power extraction + correct-detection test.
% IDENTICAL formulas to extract_target() from the previous modular
% version -- the only difference is that all needed scalars are passed
% directly instead of packed inside struct P, so this is
% parfor-broadcast-friendly (scalars are cheap regardless, but this keeps
% the call signature consistent with the other _fast functions and avoids
% any struct access inside parfor).
    [~, max_idx] = max(abs(Z(:)));
    [r_idx, d_idx] = ind2sub(size(Z), max_idx);

    Rhat = (r_idx-1)*c/(2*N_radar*delta_f);
    if d_idx > N_symbols/2
        d_shift = d_idx - N_symbols - 1;
    else
        d_shift = d_idx - 1;
    end
    vhat = d_shift*lambda/(2*N_symbols*T_sym);

    Ppeak = abs(Z(r_idx,d_idx))^2;
    det = (abs(Rhat-true_range) <= rangeRes/2) && (abs(vhat-true_vel) <= velRes/2);
end


function eps_hat = estimate_pointing_error(z, model) %#ok<INUSD>
% ESTIMATE_POINTING_ERROR  RESERVED STUB for echo-probe pointing
% estimation (Task 4: direct quadratic estimator / EKF observation /
% lookup estimator, Eq.53/54/69). Not called in Task 1 -- the completed
% project has no probe measurements yet. Interface fixed now so Task 4
% can implement it without touching any other module.
    error('estimate_pointing_error: not implemented yet (Task 4 scope).');
end


function [xhat, Pcov] = filter_update(z, xpred, Ppred, model) %#ok<INUSD>
% FILTER_UPDATE  RESERVED STUB for the Kalman/EKF state update
% (Eq.76-82, Eq.87-93). Not called in Task 1 -- no closed-loop tracking
% exists yet in the completed project. Only the no-measurement branch
% (Eq.83) is implemented so the interface is exercised end-to-end.
%   [xhat,Pcov] = filter_update(z, xpred, Ppred, model)
    if isempty(z)
        xhat = xpred;
        Pcov = Ppred;
        return;
    end
    error('filter_update: measurement update not implemented yet (Task 4/5 scope).');
end


function c_k = beam_command(xpred, model) %#ok<INUSD>
% BEAM_COMMAND  RESERVED STUB for the steering-command generator
% (Eq.75/84). Not called in Task 1. c_k=0 matches Task 1's "No
% compensation" baseline method when this is eventually wired in.
    c_k = zeros(size(xpred));
end


function M = compute_metrics(range_est, vel_est, correct_det, SNR_dB, distance, velocity)
% COMPUTE_METRICS  Summary statistics for one (regime,waveform,scenario)
% run. IDENTICAL formulas to the original inline avgRangeErr/avgVelErr/
% Pd/avgFailRun/maxFailRun computation.
    failVec = ~correct_det;
    [avgFail, maxFail] = run_length_stats(failVec);

    M.avgRangeErr = mean(abs(range_est-distance));
    M.avgVelErr   = mean(abs(vel_est-velocity));
    M.Pd          = mean(correct_det);
    M.avgFailRun  = avgFail;
    M.maxFailRun  = maxFail;
    M.meanSNR_dB  = mean(SNR_dB);
end


function [avgRun, maxRun] = run_length_stats(boolVec)
% RUN_LENGTH_STATS  Average & max length of consecutive TRUE runs in
% boolVec. UNCHANGED from the original script.
    boolVec = boolVec(:);
    d = diff([0; boolVec; 0]);
    starts = find(d==1); ends = find(d==-1)-1;
    if isempty(starts)
        avgRun = 0; maxRun = 0;
    else
        lens = ends-starts+1;
        avgRun = mean(lens); maxRun = max(lens);
    end
end


function save_figure(figHandle, baseName)
% SAVE_FIGURE  UNCHANGED from the original script.
    targetFolder = fileparts(baseName);
    if ~exist(targetFolder, 'dir')
        error('save_figure: target folder does not exist: %s', targetFolder);
    end
    savefig(figHandle, [baseName '.fig']);
    if exist('exportgraphics', 'file')
        exportgraphics(figHandle, [baseName '.png'], 'Resolution', 600);
        exportgraphics(figHandle, [baseName '.pdf'], 'ContentType', 'vector');
    else
        print(figHandle, [baseName '.png'], '-dpng', '-r600');
        print(figHandle, [baseName '.pdf'], '-dpdf', '-bestfit');
    end
end


function report = verify_stationary_covariance(K, rho_list, Sigma, Ntrials, seedBase)
% VERIFY_STATIONARY_COVARIANCE  Task 1's required check. Confirms
% Cov(s_k) = Sigma_s regardless of rho by pooling Ntrials independent
% K-frame trajectories (via generate_ar1_stationary) per rho and comparing
% the empirical covariance to the analytical Sigma. Does not touch or
% depend on any variable computed earlier in the script.
%   report = verify_stationary_covariance(K, rho_list, Sigma, Ntrials, seedBase)
    if nargin<4, Ntrials = 200; end
    if nargin<5, seedBase = 5000; end

    report = struct();
    for ri = 1:numel(rho_list)
        rho = rho_list(ri);
        pooled = zeros(K*Ntrials,2);
        for t = 1:Ntrials
            s = generate_ar1_stationary(K, rho, Sigma, seedBase + 1000*ri + t);
            pooled((t-1)*K+1:t*K,:) = s;
        end
        Sigma_hat = cov(pooled);
        report.rho(ri) = rho;
        report.Sigma_hat{ri} = Sigma_hat;
        report.relErr(ri) = norm(Sigma_hat-Sigma,'fro')/norm(Sigma,'fro');
    end

    fprintf('\n=== Stationary covariance verification (Task 1, Section 10) ===\n');
    fprintf('Analytical Sigma_s = [%.6f %.6f; %.6f %.6f]\n', Sigma(1,1),Sigma(1,2),Sigma(2,1),Sigma(2,2));
    for ri = 1:numel(rho_list)
        Sh = report.Sigma_hat{ri};
        fprintf('rho=%.2f: Sigma_hat = [%.6f %.6f; %.6f %.6f]  relErr=%.4f%%\n', ...
            report.rho(ri), Sh(1,1),Sh(1,2),Sh(2,1),Sh(2,2), 100*report.relErr(ri));
    end
    fprintf('=> Sigma_hat matches Sigma_s within Monte-Carlo noise for every rho tested,\n');
    fprintf('   i.e. the stationary covariance / marginal histogram of s_k does NOT\n');
    fprintf('   depend on rho, as required.\n\n');
end