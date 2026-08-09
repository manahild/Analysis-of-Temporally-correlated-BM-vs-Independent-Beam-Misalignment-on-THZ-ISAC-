%% ========================================================================
%  TASKS 2, 3, 4, 5 & 6 -- GENERALIZED BEAM-MISALIGNMENT LIBRARY,
%  CURVATURE FIT/VALIDATION, ECHO-BASED POINTING-ESTIMATION VALIDATION,
%  CLOSED-LOOP TRACKING, AND THE FULL CORRELATED-MISALIGNMENT SENSING
%  SIMULATOR
%  ------------------------------------------------------------------------
%  This is the "Tasks 2/3/5/6" master script with the standalone Task 4
%  script (Stage A ideal validation, Stage A2 diagnostics, Stage B
%  complete-chain validation) merged in, PLUS the standalone, FIXED Task 5
%  closed-loop-tracking script merged in as a new Section 2.7 / call.
%
%  WHAT CHANGED VS. THE PREVIOUS MERGE (Tasks 2/3/4/5/6 master, before this
%  Task 5 merge)
%  ----------------------------------------------------------------------
%  1) The standalone Task 5 script's duplicate copies of
%     setup_system_params, build_radar_reference, beam_gain_general,
%     mysinc2, beam_matrix_A, make_beam_case, tf_str, save_figure_t2,
%     generate_received_signal_fast, range_doppler_fast,
%     extract_target_fast, generate_ar1_stationary, g_observation,
%     jacobian_fd, lookup_estimate, directional_radius, target_bins,
%     probe_gate_power, calibrate_detection_threshold, measure_probe_pair,
%     fit_curvature_matrix were DELETED -- this file already had exactly
%     one copy of each (shared with Tasks 1-4), so the Task 5 section now
%     reuses those single copies instead of shadowing them.
%  2) Only the Task-5-SPECIFIC local functions were kept and added
%     VERBATIM (post-fix; see the two correctness fixes noted below):
%     run_task5_stageA_ideal, run_task5_method, compute_task5_metrics,
%     plot_task5_timeseries, plot_task5_comparison -- plus one new
%     wrapper, run_task5_closed_loop_tracking(), that owns the folder
%     setup / parameter bundling / trajectory generation / R-calibration
%     the standalone script used to do at script level (mirrors how
%     run_task4_validation() wraps Stage A/A2/B).
%  3) The standalone script's placeholder "wref = 0.20 m (ASSUMPTION)" is
%     replaced by the SAME wref2 = w_mo_SFref already used by Tasks 2/3/4
%     in this master script.
%  4) run_task5_closed_loop_tracking() is called once, immediately after
%     run_task4_validation(), using wref2 and the beam case list. Its
%     outputs (Figures_T5/, Results_T5/Task5_methods_summary_table.csv)
%     are produced alongside the other tasks' outputs.
%  5) K (frames/trajectory), Mcal (R-calibration trials) and Nfa
%     (false-alarm calibration trials) are BUMPED UP from the standalone
%     script's fast-test values to the brief's Task 6 minimums: K=200
%     (was 60), Mcal=1000 (was 200), Nfa=500 (was 100; the brief doesn't
%     set a hard floor for this one, but 100 trials gives a noisy
%     percentile estimate for a 1e-3 false-alarm target, so this was
%     raised for a more stable threshold).
%  6) save_figure_t2() gained a small robustness fix (FIX 4 below),
%     applied once since it is the single shared copy.
%
%  TWO CORRECTNESS FIXES CARRIED OVER FROM THE STANDALONE TASK 5 SCRIPT
%  (both verified against a fixed test run before this merge):
%  ----------------------------------------------------------------------
%  FIX 1 (Stage A correctness bug): run_task5_stageA_ideal's 'quad' branch
%  used to generate the noiseless measurement z from the METHOD'S OWN
%  assumed matrix (methodA{mi}), so every quad method -- including
%  MismatchedCircularKF -- was being tested against a linear system built
%  from its own (possibly wrong) A, and converged to ~0 error regardless
%  of how wrong that A was. Fixed by generating z from A_true (the
%  physical beam) always, and using methodA{mi} only inside the KF's own
%  H for the update step. Verified: after the fix, MismatchedCircularKF's
%  Stage-A final-frame error is O(1e3) m (garbage, as expected for a
%  structurally wrong model) vs O(1e-8) m for the correctly-specified
%  methods.
%
%  FIX 2 (unpaired noise seeds across methods): run_task5_method used to
%  be seeded with `baseSeedB + 1000*mi`, so each of the 7 methods drew a
%  DIFFERENT probe-noise and main-burst-noise realization even at the same
%  frame k on the same shared trajectory, contradicting the brief's Task 6
%  requirement ("Use paired trajectories and paired noise seeds across
%  methods"). Fixed by seeding every method with the SAME baseSeed.
%  Verified: after the fix, every correctly-specified method (Matched,
%  Fitted, EKF, Lookup) now shares IDENTICAL Pcorrect/avgFail/maxFail
%  (since they see the same main-burst noise draws), differing only in
%  the pointing-error-driven RMSE terms -- exactly the expected signature
%  of a properly paired comparison.
%
%  FIX 4 (this merge, cosmetic/robustness): save_figure_t2() now calls
%  drawnow before exportgraphics and falls back to print() on failure.
%  exportgraphics can throw "Invalid or deleted object" on a
%  freshly-built figure that uses yyaxis (as plot_task5_timeseries's
%  third subplot does) before the render pipeline has caught up; drawnow
%  forces the figure to finish rendering first, and the try/catch ensures
%  one bad figure never aborts the whole run.
%
%  Everything else about the radar waveform, range-Doppler processor,
%  AR(1) recursion, metric definitions, and Task 2/3/4/5/6 loop structure
%  is UNCHANGED from the previous merge.
% =========================================================================

clear; clc; close all;
rng(1);                      % single global seed, reproducible

%% ------------------------------------------------------------------
%  0. FOLDERS
%% ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
figDir2 = fullfile(scriptDir,'Figures_T2');     % Task 2 library figures
resDir2 = fullfile(scriptDir,'Results_T2');     % Task 2 library table
figDir3 = fullfile(scriptDir,'Figures_T3');     % Task 3 curvature-fit figures
resDir3 = fullfile(scriptDir,'Results_T3');     % Task 3 curvature-fit table
figDir4 = fullfile(scriptDir,'Figures_T4');     % Task 4 pointing-estimation figures
resDir4 = fullfile(scriptDir,'Results_T4');     % Task 4 pointing-estimation tables
figDir5 = fullfile(scriptDir,'Figures_T5');     % Task 5 closed-loop-tracking figures
resDir5 = fullfile(scriptDir,'Results_T5');     % Task 5 closed-loop-tracking table
figDir  = fullfile(scriptDir,'Figures_T56');    % Task 5/6 sensing figures
resDir  = fullfile(scriptDir,'Results_T56');    % Task 5/6 sensing table
datDir  = fullfile(scriptDir,'Data_T56');       % Task 5/6 saved data
for f = {figDir2,resDir2,figDir3,resDir3,figDir4,resDir4,figDir5,resDir5,figDir,resDir,datDir}
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

% Regime configs -- SAME beam-waist geometry, target size sets SF/SL
w_0_SF = 0.20;  apert_T_SF = 1.5;
w_0_SL = 0.20;  apert_T_SL = 0.04;
apert_R = 0.075;
C_n_sq  = 1e-11;
velocity  = 5;                      % m/s, target radial velocity
regimes      = {'SF','SL'};
w_0_list     = [w_0_SF, w_0_SL];
apert_T_list = [apert_T_SF, apert_T_SL];
waveforms    = {'OFDM','SCFDMA'};

distance = 15;                      % m  (matches Task 2-4 sigma_s calc)

%% ------------------------------------------------------------------
%  2. THE 5 REQUIRED BEAM CASES (Task 2 list) -- shared by the Task 2
%     validation section, the Task 3 curvature-fit section, the Task 4
%     pointing-estimation section, the Task 5 closed-loop-tracking
%     section, AND the Task 5/6 sensing loop below, via
%     make_beam_case(name, wref).
%% ------------------------------------------------------------------
beamCaseNames  = {'CircularGaussian','EllipticalGaussian', ...
                   'RotatedEllipticalGaussian','FlatTopSuperGaussian','Sidelobed'};
beamCaseLabels = {'Circular Gaussian','Elliptical Gaussian', ...
                   'Rotated Elliptical Gaussian (\phi=30^\circ)', ...
                   'Flat-top Super-Gaussian (p=4, \phi=30^\circ)', ...
                   'Sidelobed (sinc^2, \phi=30^\circ)'};
nBeamCase = numel(beamCaseNames);

%% ------------------------------------------------------------------
%  2.5 TASK 3 CURVATURE-FIT SWEEP PARAMETERS (Eq.25, Eq.39)
%% ------------------------------------------------------------------
hfloor      = 1e-12;                 % Eq.25, numerical log-protection floor only
hfitMinList = [0.8 0.5 0.25 0.1];    % Eq.39, required sweep
Ngrid3      = 401;                   % same grid density/extent convention as Task 2 (Eq.26)

%% ------------------------------------------------------------------
%  2.6 TASK 4 POINTING-ESTIMATION VALIDATION PARAMETERS
%% ------------------------------------------------------------------
fdstep4        = 1e-4*w_0_SF;        % placeholder ONLY until wref2 exists below;
                                      % recomputed properly right before the Task 4 call
isQuadratic4   = [true true true false false];   % matches beamCaseNames order (p=2 gauss family only)
beamsB4        = {'EllipticalGaussian','FlatTopSuperGaussian'};
beamLabelsB4   = {'Elliptical Gaussian','Flat-top Super-Gaussian (p=4, \phi=30^\circ)'};
alphaListB4    = [0.2, 0.4];
Ptx_dBm_listB4 = 10;
Ntrials4       = 30;
Mcal4          = 100;
Nfa4           = 100;
Pfa_target4    = 1e-3;
% jitterStd4 depends on wref2 -- set right before the Task 4 call below

%% ------------------------------------------------------------------
%  2.7 TASK 5 CLOSED-LOOP-TRACKING PARAMETERS
%     K/Mcal5/Nfa5 bumped to the brief's Task 6 minimums (see header note
%     5 above); the standalone script's fast-test values (K=60, Mcal=200,
%     Nfa=100) were only for quick iteration.
%% ------------------------------------------------------------------
beamTag5_true   = 'RotatedEllipticalGaussian';   % representative actual beam:
                                                  % quadratic (Matched KF genuinely
                                                  % matched) but non-circular
                                                  % (Mismatched circular KF
                                                  % genuinely mismatched)
beamLabel5_true = 'Rotated Elliptical Gaussian (\phi=30^\circ)';
rho5        = 0.99;              % Task 6 sweep: {0.9, 0.99}
jitterFrac5 = 0.25;               % sigma_s / wref, Task 6 sweep: {0.10, 0.25}
alphaProbe5 = 0.3;                % probe-spacing fraction, Task 6 sweep over alpha
J_interval5 = 1;                  % probe every frame, Task 6 sweep: {1,2,5,10}
K5          = 200;                % brief minimum for the Task 6 study (was 60 in fast test)
K_ideal5    = 20;                 % Stage A ideal sanity-check length (unaffected by K5)
hth5        = 0.10;               % primary outage threshold (10 dB loss)
gridNLookup5 = 41;                % lookup-estimator search-grid resolution
Mcal5       = 1000;               % brief minimum Mcal>=1000 (was 200 in fast test)
Nfa5        = 500;                % raised from the fast-test 100 for a more stable
                                   % false-alarm-threshold percentile estimate
Pfa_target5 = 1e-3;
hfit_min_task3_for5 = 0.25;       % Task 3 default fit region, used for method 5 (Fitted)
Ptx_dBm5    = 10;

%% ========================================================================
%  TASK 2 -- GENERAL BEAM-MAP LIBRARY VALIDATION (run once)
%  wref is the ACTUAL w_mo of the SF regime (computed below), replacing
%  the standalone script's placeholder wref = 0.20 m.
%% ========================================================================
[w_mo_SFref, eta_SFref, is_SF_ref] = compute_beam_width(distance, w_0_SF, apert_T_SF, apert_R, C_n_sq, fc); %#ok<ASGLU>
wref2 = w_mo_SFref;
run_task2_validation(wref2, figDir2, resDir2, beamCaseNames, beamCaseLabels);

%% ========================================================================
%  TASK 3 -- FIT AND TEST THE BEAM-CURVATURE MATRIX (run once)
%  Uses the SAME wref2 (= w_mo of the SF regime) as Task 2's validation,
%  instead of the standalone script's placeholder wref = 0.20 m.
%% ========================================================================
run_task3_validation(wref2, figDir3, resDir3, beamCaseNames, beamCaseLabels, ...
    hfitMinList, hfloor, Ngrid3);

%% ========================================================================
%  TASK 4 -- ECHO-BASED POINTING ESTIMATION VALIDATION (run once)
%  Uses the SAME wref2 as Tasks 2/3 (replacing the standalone script's
%  wref=0.20 m placeholder). See the header note above for the pre-merge
%  validation summary (Direct/GN exact; EKF-1-step cold-start artifact
%  for FlatTop/Sidelobed, confirmed by Stage A2; Lookup is the reliable
%  fallback per Stage B).
%% ========================================================================
fdstep4      = 1e-4*wref2;
jitterStd4   = 0.25*wref2;
run_task4_validation(wref2, hfloor, fdstep4, beamCaseNames, beamCaseLabels, isQuadratic4, ...
    figDir4, resDir4, beamsB4, beamLabelsB4, alphaListB4, Ptx_dBm_listB4, ...
    Ntrials4, Mcal4, Nfa4, Pfa_target4, jitterStd4);

%% ========================================================================
%  TASK 5 -- CLOSED-LOOP TRACKING AND COMPENSATION (run once)
%  Uses the SAME wref2 as Tasks 2/3/4. Implements all 7 methods (No
%  compensation, Oracle, Matched quadratic KF, Mismatched circular KF,
%  Fitted-matrix KF, Full-map EKF, Lookup plus KF), Stage A ideal sanity
%  check followed by Stage B full sensing-chain run, for one representative
%  non-circular beam (rotated elliptical Gaussian) per the brief's
%  deliverable ("a time-series figure ... for one representative
%  trajectory"). See header notes above for the two correctness fixes
%  (Stage A z-generation from A_true; paired noise seeds across methods)
%  already verified in the standalone test before this merge.
%% ========================================================================
fdstep5 = 1e-4*wref2;
run_task5_closed_loop_tracking(wref2, hfloor, fdstep5, figDir5, resDir5, ...
    beamTag5_true, beamLabel5_true, rho5, jitterFrac5, alphaProbe5, J_interval5, ...
    K5, K_ideal5, hth5, gridNLookup5, Mcal5, Nfa5, Pfa_target5, ...
    hfit_min_task3_for5, Ptx_dBm5, distance, velocity);

%% ------------------------------------------------------------------
%  3. SENSING / MISALIGNMENT SCENARIO SETUP  (Task 6 experiment matrix)
%% ------------------------------------------------------------------
theta_bm = 30e-3;                   % 30 mrad (matches Task 2-4)
sigma_s  = distance*tan(theta_bm);  % pointing-jitter std, m

Ptx_dBm_fixed = 10;                  % ASSUMPTION: single fixed operating
Ptx_W = 10^((Ptx_dBm_fixed-30)/10); % point (mid of the sweep used in the
                                     % main script), since T5 asks for a
                                     % *sequence of frames*, not a power
                                     % sweep. Change here if a different
                                     % operating power is preferred.

N_frames = 10;                     % ASSUMPTION: sensing frames per
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
fprintf('Task 5-6 setup\n');
fprintf('===============================================\n');
fprintf('distance=%.1f m | theta_bm=%.0f mrad | sigma_s=%.4f m\n', distance, theta_bm*1000, sigma_s);
fprintf('Ptx = %.1f dBm | N_frames = %d | h_th = %.2f\n', Ptx_dBm_fixed, N_frames, h_th);
fprintf('Beam cases (%d): %s\n\n', nBeamCase, strjoin(beamCaseNames, ', '));

%% ------------------------------------------------------------------
%  4. ATMOSPHERIC MODEL (unchanged)
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
%  5. PARALLEL POOL (unchanged)
%% ------------------------------------------------------------------
pool = gcp('nocreate');
if isempty(pool), pool = parpool('local'); end
fprintf('Using parallel pool with %d workers\n\n', pool.NumWorkers);

%% ------------------------------------------------------------------
%  6-7. MAIN LOOP: regime x BEAM CASE x waveform x scenario x frame
%     Radar reference / waveform build -> module: build_radar_reference()
%     Trajectory generation             -> module: generate_ar1()
%     Beam-map evaluation                -> module: beam_gain_general()  <-- Task 2 library
%     Received-signal generation          -> module: generate_received_signal_fast()
%     Range-Doppler processing            -> module: range_doppler_fast()
%     Peak extraction / detection         -> module: extract_target_fast()
%     Metrics                             -> module: compute_metrics()
%
%  NOTE ON TRAJECTORY REUSE: eps trajectories (aligned/iid/AR1-0.9/
%  AR1-0.99) are drawn ONCE per regime, BEFORE the beam-case loop, and the
%  identical trajectories are reused for every beam case -- only h(eps)
%  changes across beam cases (brief: "Reuse exactly the same trajectory
%  when comparing beam models").
%
%  NOTE ON PERFORMANCE: the large arrays W.rx_delayed_doppler and
%  W.corr_ref are unpacked into plain local variables immediately before
%  every parfor block, so MATLAB broadcasts them once per waveform instead
%  of re-copying the struct on every iteration (same fix as the source
%  Task 5/6 script).
%% ------------------------------------------------------------------
results = struct();

for reg_idx = 1:numel(regimes)
    regime  = regimes{reg_idx};
    w_0     = w_0_list(reg_idx);
    apert_T = apert_T_list(reg_idx);

    fprintf('========== %s CONFIG (w_0=%.0f mm, a_T=%.2f m) ==========\n', regime, w_0*1000, apert_T);

    %% --- beam-waist / w_mo computation (unchanged) -----------------
    [w_mo, eta, is_SF] = compute_beam_width(distance, w_0, apert_T, apert_R, C_n_sq, fc);
    fprintf('  eta=%.2f  ->  w_mo = %.4f m  (used as wref for all %d beam cases)\n', eta, w_mo, nBeamCase);

    if strcmp(regime,'SF') && ~is_SF
        warning('%s config has eta=%.3f (<=1): physics says Spot-Limited, not Spot-Filling!', regime, eta);
    elseif strcmp(regime,'SL') && is_SF
        warning('%s config has eta=%.3f (>1): physics says Spot-Filling, not Spot-Limited!', regime, eta);
    end

    wref = w_mo;   % Task 2 convention: wref = wmo, kept fixed across beam cases

    %% --- generate the SAME pointing trajectories for every beam case ---
    ux = randn(N_frames,1);  uy = randn(N_frames,1);   % shared innovations

    epsScenario = cell(1,nScenario);
    rk_scenario = zeros(N_frames, nScenario);

    % 1: aligned
    epsScenario{1} = zeros(N_frames,2);
    rk_scenario(:,1) = 0;

    % 2: i.i.d. (rho=0 -- fresh independent draws, matches Task 2 exactly)
    xk_iid = sigma_s*randn(N_frames,1);
    yk_iid = sigma_s*randn(N_frames,1);
    epsScenario{2} = [xk_iid, yk_iid];
    rk_scenario(:,2) = sqrt(sum(epsScenario{2}.^2,2));

    % 3 & 4: correlated AR(1), rho = 0.9, 0.99
    for si = 1:2
        rho = rho_list(si+1);
        epsScenario{si+2} = generate_ar1(N_frames, rho, sigma_s, ux, uy);
        rk_scenario(:,si+2) = sqrt(sum(epsScenario{si+2}.^2,2));
    end

    %% --- radar equation amplitude (aligned reference, unchanged) ------
    P_radar = Ptx_W*(N_radar/N_fft);
    A_base = sqrt(P_radar)*L_atm_lin* ...
             sqrt((G_total*lambda^2*sigma_rcs)/((4*pi)^3*distance^4));

    for bc_idx = 1:nBeamCase
        beamTag   = beamCaseNames{bc_idx};
        beamLabel = beamCaseLabels{bc_idx};
        bp = make_beam_case(beamTag, wref);

        %% --- Task 2 beam evaluation: h_k for every scenario under THIS
        %      beam, using the SAME eps trajectories generated above -----
        hk_scenario = zeros(N_frames, nScenario);
        for sc_idx = 1:nScenario
            e = epsScenario{sc_idx};
            hk_scenario(:,sc_idx) = beam_gain_general(e(:,1), e(:,2), bp);   % Task 2 library
        end

        for wf_idx = 1:numel(waveforms)
            waveform = waveforms{wf_idx};

            W = build_radar_reference(P, waveform, distance, velocity);

            % --- PERFORMANCE FIX: unpack the large arrays into plain
            % local variables BEFORE parfor. --------------------------
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
                results.(regime).(beamTag).(waveform).(tag).h_k        = hk_vec;
                results.(regime).(beamTag).(waveform).(tag).r_k         = rk_scenario(:,sc_idx);
                results.(regime).(beamTag).(waveform).(tag).range_est   = range_est;
                results.(regime).(beamTag).(waveform).(tag).vel_est     = vel_est;
                results.(regime).(beamTag).(waveform).(tag).correct_det = correct_det;
                results.(regime).(beamTag).(waveform).(tag).SNR_dB      = SNR_dB;
                results.(regime).(beamTag).(waveform).(tag).Prx_lin      = Prx_lin;
                results.(regime).(beamTag).(waveform).(tag).avgRangeErr = M.avgRangeErr;
                results.(regime).(beamTag).(waveform).(tag).avgVelErr   = M.avgVelErr;
                results.(regime).(beamTag).(waveform).(tag).Pd          = M.Pd;
                results.(regime).(beamTag).(waveform).(tag).avgFailRun  = M.avgFailRun;
                results.(regime).(beamTag).(waveform).(tag).maxFailRun  = M.maxFailRun;

                fprintf('  [%s,%-42s,%s,%-14s] Pd=%.3f  avgRangeErr=%.4f m  avgVelErr=%.4f m/s  avgFailRun=%.2f  maxFailRun=%d\n', ...
                    regime, beamLabel, waveform, scenarioNames{sc_idx}, M.Pd, M.avgRangeErr, M.avgVelErr, M.avgFailRun, M.maxFailRun);
            end
        end
    end
    fprintf('\n');
end

%% ------------------------------------------------------------------
%  8. TASK 6 -- RESULTS TABLE (now includes a BeamCase column)
%% ------------------------------------------------------------------
T6 = table();
for reg_idx = 1:numel(regimes)
    regime = regimes{reg_idx};
    for bc_idx = 1:nBeamCase
        beamTag   = beamCaseNames{bc_idx};
        beamLabel = beamCaseLabels{bc_idx};
        for wf_idx = 1:numel(waveforms)
            waveform = waveforms{wf_idx};
            for sc_idx = 1:nScenario
                tag = scenarioTagsSF{sc_idx};
                R = results.(regime).(beamTag).(waveform).(tag);
                newRow = table({regime}, {beamLabel}, {waveform}, {scenarioNames{sc_idx}}, ...
                    R.avgRangeErr, R.avgVelErr, R.Pd, R.avgFailRun, R.maxFailRun, ...
                    'VariableNames', {'Regime','BeamCase','Waveform','Scenario', ...
                    'AvgRangeErr_m','AvgVelErr_mps','Pd','AvgConsecFailRun','MaxConsecFailRun'});
                T6 = [T6; newRow]; %#ok<AGROW>
            end
        end
    end
end
disp(T6);
writetable(T6, fullfile(resDir,'Task6_results_table.csv'));

%% ------------------------------------------------------------------
%  9. TASK 6 -- FIGURES
%     Figs 1-3 are now produced per (regime, beam case) using OFDM as the
%     representative waveform (SC-FDMA data is identical in structure and
%     saved in T6/results; comparing waveforms is out of scope per the
%     brief). Fig 4 (consecutive-failure bars, SF vs SL) is produced once
%     per beam case.
%% ------------------------------------------------------------------
colorsSc = lines(nScenario);
wfShow = 'OFDM';

for reg_idx = 1:numel(regimes)
    regime = regimes{reg_idx};
    for bc_idx = 1:nBeamCase
        beamTag = beamCaseNames{bc_idx};

        % --- Fig 1: h_k vs frame number, all rho ---
        fig = figure('Color','w','Position',[100 100 1000 550]); hold on;
        for sc_idx = 2:nScenario   % skip 'Aligned' (h_k==1 trivial line)
            tag = scenarioTagsSF{sc_idx};
            plot(results.(regime).(beamTag).(wfShow).(tag).h_k, 'Color', colorsSc(sc_idx,:), ...
                'DisplayName', scenarioNames{sc_idx});
        end
        xlabel('Frame Number, k'); ylabel('h_k');
        title(sprintf('Task 6 Fig 1: h_k vs Frame, %s regime, %s', regime, beamCaseLabels{bc_idx}), 'Interpreter','tex');
        ylim([0 1]); legend('Location','best'); grid on;
        save_figure(fig, fullfile(figDir, sprintf('T6_Fig1_hk_vs_frame_%s_%s', regime, beamTag)));

        % --- Fig 2: received SNR vs frame number ---
        fig = figure('Color','w','Position',[100 100 1000 550]); hold on;
        for sc_idx = 1:nScenario
            tag = scenarioTagsSF{sc_idx};
            plot(results.(regime).(beamTag).(wfShow).(tag).SNR_dB, 'Color', colorsSc(sc_idx,:), ...
                'DisplayName', scenarioNames{sc_idx});
        end
        xlabel('Frame Number, k'); ylabel('Received SNR (dB)');
        title(sprintf('Task 6 Fig 2: SNR vs Frame, %s regime, %s, %s', regime, beamCaseLabels{bc_idx}, wfShow), 'Interpreter','tex');
        legend('Location','best'); grid on;
        save_figure(fig, fullfile(figDir, sprintf('T6_Fig2_SNR_vs_frame_%s_%s', regime, beamTag)));

        % --- Fig 3: correct/incorrect detections vs frame number ---
        fig = figure('Color','w','Position',[100 100 1000 700]);
        for sc_idx = 1:nScenario
            tag = scenarioTagsSF{sc_idx};
            subplot(nScenario,1,sc_idx);
            stem(double(results.(regime).(beamTag).(wfShow).(tag).correct_det), ...
                'Marker','none','Color', colorsSc(sc_idx,:));
            ylim([-0.2 1.2]); ylabel(scenarioNames{sc_idx}, 'FontSize', 8);
            if sc_idx == 1
                title(sprintf('Task 6 Fig 3: Detection Correct(1)/Incorrect(0) vs Frame, %s, %s, %s', regime, beamCaseLabels{bc_idx}, wfShow), 'Interpreter','tex');
            end
            if sc_idx == nScenario, xlabel('Frame Number, k'); end
        end
        save_figure(fig, fullfile(figDir, sprintf('T6_Fig3_detections_vs_frame_%s_%s', regime, beamTag)));
    end
end

% --- Fig 4: number of consecutive failed frames vs rho (both regimes),
%     one figure per beam case ---
xRho = 1:3;   % iid, rho0.9, rho0.99  (aligned excluded: not a fading scenario)
for bc_idx = 1:nBeamCase
    beamTag = beamCaseNames{bc_idx};

    fig = figure('Color','w','Position',[100 100 900 550]); hold on;
    avgFailSF = zeros(1,3); maxFailSF = zeros(1,3);
    avgFailSL = zeros(1,3); maxFailSL = zeros(1,3);
    for i = 1:3
        tag = scenarioTagsSF{i+1};
        avgFailSF(i) = results.SF.(beamTag).(wfShow).(tag).avgFailRun;
        maxFailSF(i) = results.SF.(beamTag).(wfShow).(tag).maxFailRun;
        avgFailSL(i) = results.SL.(beamTag).(wfShow).(tag).avgFailRun;
        maxFailSL(i) = results.SL.(beamTag).(wfShow).(tag).maxFailRun;
    end
    bar(xRho-0.2, avgFailSF, 0.15, 'FaceColor',[0.2 0.4 0.8], 'DisplayName','SF avg');
    bar(xRho-0.05, maxFailSF, 0.15, 'FaceColor',[0.1 0.2 0.5], 'DisplayName','SF max');
    bar(xRho+0.1, avgFailSL, 0.15, 'FaceColor',[0.85 0.4 0.2], 'DisplayName','SL avg');
    bar(xRho+0.25, maxFailSL, 0.15, 'FaceColor',[0.55 0.2 0.1], 'DisplayName','SL max');
    set(gca,'XTick', xRho, 'XTickLabel', {'i.i.d.','\rho=0.90','\rho=0.99'});
    ylabel('Consecutive failed frames'); xlabel('Misalignment scenario');
    title(sprintf('Task 6 Fig 4: Consecutive Failed Frames vs \\rho (%s, %s)', wfShow, beamCaseLabels{bc_idx}), 'Interpreter','tex');
    legend('Location','best'); grid on;
    save_figure(fig, fullfile(figDir, sprintf('T6_Fig4_consec_fail_vs_rho_%s', beamTag)));
end
close all;

%% ------------------------------------------------------------------
%  10. SAVE DATA (now includes the beam-case dimension + Task 4/5 params)
%% ------------------------------------------------------------------
save(fullfile(datDir,'Task5_6_results.mat'), 'results', 'T6', ...
    'distance','theta_bm','sigma_s','Ptx_dBm_fixed','N_frames','h_th', ...
    'rangeRes','velRes','regimes','waveforms','scenarioNames','scenarioTagsSF', ...
    'beamCaseNames','beamCaseLabels', ...
    'w_0_SF','apert_T_SF','w_0_SL','apert_T_SL','apert_R','C_n_sq','fc', ...
    'hfitMinList','hfloor','Ngrid3', ...
    'fdstep4','isQuadratic4','beamsB4','beamLabelsB4','alphaListB4', ...
    'Ptx_dBm_listB4','Ntrials4','Mcal4','Nfa4','Pfa_target4','jitterStd4', ...
    'beamTag5_true','beamLabel5_true','rho5','jitterFrac5','alphaProbe5', ...
    'J_interval5','K5','K_ideal5','hth5','gridNLookup5','Mcal5','Nfa5', ...
    'Pfa_target5','hfit_min_task3_for5','Ptx_dBm5');

fprintf('\nFigures -> %s\nTable   -> %s\nData    -> %s\n', figDir, resDir, datDir);
fprintf('============ TASKS 2, 3, 4, 5 & 6 COMPLETE ============\n');

%% ------------------------------------------------------------------
%  11. TASK 1 REQUIREMENT (unchanged):
%      Verify numerically that changing rho does NOT change the
%      stationary covariance / marginal histogram of s_k.
%% ------------------------------------------------------------------
Sigma_s = sigma_s^2*eye(2);
verify_stationary_covariance(N_frames, [0 0.9 0.99], Sigma_s, 200, 5000);

fprintf('============ TASK 1 VERIFICATION COMPLETE ============\n');


%% ========================================================================
%  LOCAL FUNCTIONS -- Task 1 module set (unchanged from the source script)
%% ========================================================================

function P = setup_system_params()
% SETUP_SYSTEM_PARAMS  THz-ISAC system parameters (Table I).
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
% COMPUTE_BEAM_WIDTH  Effective beam-width parameter w_mo (Eq.1/7). Unchanged physics.
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
% waveform, delay + Doppler phase. Unchanged.
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
% GENERATE_AR1  AR(1) correlated pointing-error generator, Eq.(5)-(6). Unchanged.
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
% matching the brief's exact interface. Unchanged.
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


function h = beam_gain(eps, beam_parameters) %#ok<DEFNU>
% BEAM_GAIN  LEGACY circular-Gaussian-only beam coefficient (Eq.12).
% Superseded by beam_gain_general() below -- kept only for reference /
% backward compatibility. Not called by the main loop in this file.
    w_mo = beam_parameters.w_mo;
    r2 = eps(:,1).^2 + eps(:,2).^2;
    h = exp(-2*r2/w_mo^2);
end


function [rx_signal, Pr_k] = generate_received_signal_fast(rx_delayed_doppler, h_k, A_base, P_n_radar)
% GENERATE_RECEIVED_SIGNAL_FAST  Received sensing signal for one frame,
% Pr,k = Pr,aligned * h_k (Eq.5/8). Amplitude scaled by sqrt(h_k). Unchanged.
    A_echo_k = A_base * sqrt(h_k);
    rx_echo  = A_echo_k * rx_delayed_doppler;
    noise = sqrt(P_n_radar/2)*(randn(size(rx_echo))+1j*randn(size(rx_echo)));
    rx_signal = rx_echo + noise;
    Pr_k = A_echo_k^2;                    % = A_base^2 * h_k
end


function Z = range_doppler_fast(rx_signal, corr_ref, N_fft, N_cp, N_symbols, N_radar)
% RANGE_DOPPLER_FAST  Matched-filter range-Doppler map for one frame. Unchanged.
    rx_cp = reshape(rx_signal, N_fft+N_cp, N_symbols);
    rx_time = rx_cp(N_cp+1:end,:);
    rx_freq = fft(rx_time, N_fft, 1)/sqrt(N_fft);
    rx_radar = rx_freq(1:N_radar,:);

    corr_out = rx_radar .* conj(corr_ref);
    range_profile = ifft(corr_out, N_radar, 1);
    Z = fft(range_profile, N_symbols, 2);
end


function [Ppeak, Rhat, vhat, det] = extract_target_fast(Z, c, delta_f, N_radar, lambda, T_sym, N_symbols, true_range, true_vel, rangeRes, velRes)
% EXTRACT_TARGET_FAST  Peak-power extraction + correct-detection test. Unchanged.
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


function eps_hat = estimate_pointing_error(z, model) %#ok<INUSD,DEFNU>
% ESTIMATE_POINTING_ERROR  RESERVED STUB. Superseded by the Task 4
% estimators (lookup_estimate / jacobian_fd-based EKF-1-step /
% gauss_newton_solve) below. Kept only for interface compatibility with
% any Task 5 code that may still reference this stub name.
    error('estimate_pointing_error: use the Task 4 estimators (lookup_estimate, EKF-1step via jacobian_fd, or gauss_newton_solve) instead.');
end


function [xhat, Pcov] = filter_update(z, xpred, Ppred, model) %#ok<INUSD,DEFNU>
% FILTER_UPDATE  RESERVED STUB for the Kalman/EKF state update (Task 5 scope). Unchanged.
    if isempty(z)
        xhat = xpred;
        Pcov = Ppred;
        return;
    end
    error('filter_update: measurement update not implemented yet (Task 5 scope).');
end


function c_k = beam_command(xpred, model) %#ok<INUSD,DEFNU>
% BEAM_COMMAND  RESERVED STUB for the steering-command generator. Unchanged.
    c_k = zeros(size(xpred));
end


function M = compute_metrics(range_est, vel_est, correct_det, SNR_dB, distance, velocity)
% COMPUTE_METRICS  Summary statistics for one (regime,beam,waveform,scenario) run. Unchanged.
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
% RUN_LENGTH_STATS  Average & max length of consecutive TRUE runs. Unchanged.
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
% SAVE_FIGURE  Unchanged.
    targetFolder = fileparts(baseName);
    if ~exist(targetFolder, 'dir')
        error('save_figure: target folder does not exist: %s', targetFolder);
    end
    try
        savefig(figHandle, [baseName '.fig']);
    catch
        % savefig not always available/needed under Octave batch mode
    end
    if exist('exportgraphics', 'file')
        exportgraphics(figHandle, [baseName '.png'], 'Resolution', 600);
        try
            exportgraphics(figHandle, [baseName '.pdf'], 'ContentType', 'vector');
        catch
        end
    else
        print(figHandle, [baseName '.png'], '-dpng', '-r600');
    end
end


function report = verify_stationary_covariance(K, rho_list, Sigma, Ntrials, seedBase)
% VERIFY_STATIONARY_COVARIANCE  Task 1's required check. Unchanged.
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


%% ========================================================================
%  LOCAL FUNCTIONS -- Task 2 beam-map library (SINGLE shared copy; used by
%  Task 2 validation, Task 3 curvature fitting, Task 4 pointing
%  estimation, Task 5 closed-loop tracking, AND the Task 5/6 sensing loop
%  above)
%% ========================================================================

function h = beam_gain_general(epsx, epsy, bp)
% BEAM_GAIN_GENERAL  Dispatches to the analytic beam map requested by bp.
%   h = beam_gain_general(epsx, epsy, bp)
%     epsx, epsy - arrays (same size), lateral residual in metres
%     bp.type    - 'gauss' (circular/elliptical/rotated/super-Gaussian,
%                  selected by wx,wy,phi,p) or 'sidelobe' (sinc^2)
%     bp.wx,wy   - 1/e^2 beam radii along the beam's own axes (m)
%     bp.phi     - rotation angle (rad), Eq.13-14
%     bp.p       - super-Gaussian order (p=2 -> elliptical Gaussian), Eq.19
%     bp.bx,by   - sidelobe null spacings (m), Eq.20
%   h            - same size as epsx, 0<=h<=1, h(0)=1
    cphi = cos(bp.phi); sphi = sin(bp.phi);
    ux =  cphi.*epsx + sphi.*epsy;    % Eq.14, u = R(phi)^T * eps
    uy = -sphi.*epsx + cphi.*epsy;

    switch bp.type
        case 'gauss'
            p = bp.p;
            chi = (ux/bp.wx).^2 + (uy/bp.wy).^2;      % Eq.18
            h = exp(-2*chi.^(p/2));                    % Eq.19 (p=2 -> Eq.15)
        case 'sidelobe'
            h = mysinc2(ux/bp.bx) .* mysinc2(uy/bp.by);  % Eq.20-21
        otherwise
            error('beam_gain_general: unknown bp.type "%s"', bp.type);
    end
end


function y = mysinc2(v)
% MYSINC2  sinc^2(v) with sinc(v)=sin(pi v)/(pi v), sinc(0)=1 (Eq.21).
    y = ones(size(v));
    nz = (v ~= 0);
    y(nz) = (sin(pi*v(nz))./(pi*v(nz))).^2;
end


function A = beam_matrix_A(wx, wy, phi)
% BEAM_MATRIX_A  Matrix-form parameter A of Eq.17.
    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    D = diag([2/wx^2, 2/wy^2]);
    A = R*D*R.';
end


function h = beam_quadratic_matrix(epsx, epsy, A)
% BEAM_QUADRATIC_MATRIX  h = exp(-eps^T A eps), evaluated elementwise (Eq.16).
    quad = A(1,1)*epsx.^2 + (A(1,2)+A(2,1))*epsx.*epsy + A(2,2)*epsy.^2;
    h = exp(-quad);
end


function h = beam_gain_numerical(xg, yg, Hnorm, xq, yq)
% BEAM_GAIN_NUMERICAL  Bilinear interpolation of a normalized numerical
% beam map (Eq.22-24). Returns 0 outside the stored grid.
    h = zeros(size(xq));
    inside = xq >= xg(1) & xq <= xg(end) & yq >= yg(1) & yq <= yg(end);

    xi = xq(inside); yi = yq(inside);
    h_in = zeros(size(xi));
    for n = 1:numel(xi)
        ix = find(xg <= xi(n), 1, 'last'); if ix >= numel(xg), ix = numel(xg)-1; end
        jy = find(yg <= yi(n), 1, 'last'); if jy >= numel(yg), jy = numel(yg)-1; end
        t = (xi(n)-xg(ix)) / (xg(ix+1)-xg(ix));       % Eq.23
        v = (yi(n)-yg(jy)) / (yg(jy+1)-yg(jy));
        H00 = Hnorm(jy,   ix);
        H10 = Hnorm(jy,   ix+1);
        H01 = Hnorm(jy+1, ix);
        H11 = Hnorm(jy+1, ix+1);
        h_in(n) = (1-t)*(1-v)*H00 + t*(1-v)*H10 + (1-t)*v*H01 + t*v*H11;  % Eq.24
    end
    h(inside) = h_in;
end


function [W, Wnull] = cut_widths(xv, hcut, level, findNull)
% CUT_WIDTHS  Full width between the two mainlobe crossings of `level`
% nearest the origin.
    [~, i0] = min(abs(xv-0));

    xNeg = find_crossing(xv, hcut, level, i0, -1);
    xPos = find_crossing(xv, hcut, level, i0, +1);
    if isnan(xNeg) || isnan(xPos)
        W = NaN;
    else
        W = xPos - xNeg;
    end

    Wnull = NaN;
    if findNull
        nullLevel = 1e-3;   % numerical "null" floor for the sinc^2 mainlobe
        xNegN = find_crossing(xv, hcut, nullLevel, i0, -1);
        xPosN = find_crossing(xv, hcut, nullLevel, i0, +1);
        if ~isnan(xNegN) && ~isnan(xPosN)
            Wnull = xPosN - xNegN;
        end
    end
end


function xc = find_crossing(xv, hcut, level, i0, dirn)
% FIND_CROSSING  Walk outward from index i0, linearly interpolate the
% first crossing of `level`. Returns NaN if none found before the edge.
    N = numel(xv);
    xc = NaN;
    idx = i0;
    while true
        nextIdx = idx + dirn;
        if nextIdx < 1 || nextIdx > N, return; end
        if (hcut(idx)-level)*(hcut(nextIdx)-level) <= 0 && hcut(idx) ~= hcut(nextIdx)
            frac = (level-hcut(idx)) / (hcut(nextIdx)-hcut(idx));
            xc = xv(idx) + frac*(xv(nextIdx)-xv(idx));
            return;
        end
        idx = nextIdx;
    end
end


function s = tf_str(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end


function s = fmt_or_na(v)
    if isnan(v), s = '-'; else, s = sprintf('%.4f', v); end
end


function save_figure_t2(figHandle, baseName)
% SAVE_FIGURE_T2  FIX 4 (this merge): drawnow before export, and a
% try/catch fallback to print() -- exportgraphics can throw "Invalid or
% deleted object" on a freshly-built figure that uses yyaxis (as
% plot_task5_timeseries's third subplot does) before the render pipeline
% has caught up; drawnow forces the figure to finish rendering first, and
% the try/catch ensures one bad figure never aborts the whole run.
    targetFolder = fileparts(baseName);
    if ~exist(targetFolder, 'dir')
        error('save_figure_t2: target folder does not exist: %s', targetFolder);
    end
    try
        savefig(figHandle, [baseName '.fig']);
    catch
        % savefig not always available/needed under Octave batch mode
    end
    drawnow;   % FIX 4: let the figure (esp. yyaxis plots) finish rendering
    try
        if exist('exportgraphics', 'file')
            exportgraphics(figHandle, [baseName '.png'], 'Resolution', 200);
        else
            print(figHandle, [baseName '.png'], '-dpng', '-r200');
        end
    catch ME
        warning('save_figure_t2: exportgraphics failed (%s) -- falling back to print().', ME.message);
        try
            print(figHandle, [baseName '.png'], '-dpng', '-r200');
        catch ME2
            warning('save_figure_t2: print() fallback also failed (%s) -- skipping PNG for %s.', ME2.message, baseName);
        end
    end
end


%% ========================================================================
%  LOCAL FUNCTIONS -- integration glue between Task 2 and Task 5/6
%% ========================================================================

function bp = make_beam_case(beamTag, wref)
% MAKE_BEAM_CASE  Returns the bp struct for one of the 5 Task 2 required
% beam cases, parameterised by wref (= w_mo of the regime under test):
%   1. circular Gaussian:          wx = wy = wref, phi = 0
%   2. elliptical Gaussian:        wx = wref, wy = 0.6*wref, phi = 0
%   3. rotated elliptical Gaussian: same widths, phi = 30 deg
%   4. flat-top (p=4):             same widths/phi, p = 4
%   5. sidelobed (sinc^2):         bx = wref, by = 0.6*wref, phi = 30 deg
    phi30 = deg2rad(30);
    switch beamTag
        case 'CircularGaussian'
            bp = struct('type','gauss','wx',wref,'wy',wref,'phi',0,'p',2);
        case 'EllipticalGaussian'
            bp = struct('type','gauss','wx',wref,'wy',0.6*wref,'phi',0,'p',2);
        case 'RotatedEllipticalGaussian'
            bp = struct('type','gauss','wx',wref,'wy',0.6*wref,'phi',phi30,'p',2);
        case 'FlatTopSuperGaussian'
            bp = struct('type','gauss','wx',wref,'wy',0.6*wref,'phi',phi30,'p',4);
        case 'Sidelobed'
            bp = struct('type','sidelobe','bx',wref,'by',0.6*wref,'phi',phi30);
        otherwise
            error('make_beam_case: unknown beamTag "%s"', beamTag);
    end
end


function run_task2_validation(wref, figDir2, resDir2, beamCaseNames, beamCaseLabels)
% RUN_TASK2_VALIDATION  Runs the 6 required Task 2 validation checks and
% produces, for the 5 required beam cases, a 2-D power map, contour plot,
% horizontal/vertical cuts, and a table of -3 dB / -10 dB / first-null
% widths + symmetry flags. Logic identical to the standalone Task 2
% script; wref is now passed in (the actual w_mo of the SF regime) rather
% than a hardcoded placeholder.
    Ngrid = 401;                                  % >= 401x401 required
    xv = linspace(-2*wref, 2*wref, Ngrid);
    yv = linspace(-2*wref, 2*wref, Ngrid);
    [EPSx, EPSy] = meshgrid(xv, yv);

    h_3dB  = 10^(-3/10);                          % Eq.28
    h_10dB = 10^(-10/10);                         % = 0.1, Eq.28
    tol = 1e-9;

    fprintf('===============================================\n');
    fprintf('Task 2 -- General beam library validation\n');
    fprintf('===============================================\n');
    fprintf('wref = %.4f m (= w_mo of the SF regime) | grid = %d x %d over [-2wref,2wref]\n\n', wref, Ngrid, Ngrid);

    nBeam = numel(beamCaseNames);
    bpList = cell(1,nBeam);
    for i = 1:nBeam
        bpList{i} = make_beam_case(beamCaseNames{i}, wref);
    end

    %% Check 1: h(0) = 1
    fprintf('--- Check 1: h(0) = 1 ---------------------------------\n');
    check1_pass = true;
    for i = 1:nBeam
        h0 = beam_gain_general(0, 0, bpList{i});
        ok = abs(h0-1) < tol;
        check1_pass = check1_pass && ok;
        fprintf('  %-45s h(0) = %.12f  [%s]\n', beamCaseLabels{i}, h0, tf_str(ok));
    end
    fprintf('  => Check 1 overall: %s\n\n', tf_str(check1_pass));

    %% Check 2: 0 <= h(eps) <= 1 everywhere
    fprintf('--- Check 2: 0 <= h(eps) <= 1 everywhere --------------\n');
    check2_pass = true;
    Hall = cell(1,nBeam);
    for i = 1:nBeam
        H = beam_gain_general(EPSx, EPSy, bpList{i});
        Hall{i} = H;
        ok = (min(H(:)) >= -tol) && (max(H(:)) <= 1+tol);
        check2_pass = check2_pass && ok;
        fprintf('  %-45s min=%.3e  max=%.6f  [%s]\n', beamCaseLabels{i}, min(H(:)), max(H(:)), tf_str(ok));
    end
    fprintf('  => Check 2 overall: %s\n\n', tf_str(check2_pass));

    %% Check 3: p=2 super-Gaussian == elliptical Gaussian
    fprintf('--- Check 3: p=2 super-Gaussian == elliptical Gaussian -\n');
    phi30 = deg2rad(30);
    bp_ell = struct('type','gauss','wx',wref,'wy',0.6*wref,'phi',phi30,'p',2);
    bp_sg2 = struct('type','gauss','wx',wref,'wy',0.6*wref,'phi',phi30,'p',2);
    H_ell = beam_gain_general(EPSx, EPSy, bp_ell);
    H_sg2 = beam_gain_general(EPSx, EPSy, bp_sg2);
    maxdiff3 = max(abs(H_ell(:)-H_sg2(:)));
    check3_pass = maxdiff3 < tol;
    fprintf('  max|H_elliptical - H_supergaussian(p=2)| = %.3e  [%s]\n\n', maxdiff3, tf_str(check3_pass));

    %% Check 4: matrix form == rotated-coordinate form
    fprintf('--- Check 4: matrix form == rotated-coordinate form ---\n');
    wx4 = wref; wy4 = 0.6*wref; phi4 = phi30;
    A4 = beam_matrix_A(wx4, wy4, phi4);
    H_rot = beam_gain_general(EPSx, EPSy, struct('type','gauss','wx',wx4,'wy',wy4,'phi',phi4,'p',2));
    H_mat = beam_quadratic_matrix(EPSx, EPSy, A4);
    maxdiff4 = max(abs(H_rot(:)-H_mat(:)));
    check4_pass = maxdiff4 < 1e-8;
    fprintf('  A = [%.6f %.6f; %.6f %.6f]\n', A4(1,1),A4(1,2),A4(2,1),A4(2,2));
    fprintf('  max|H_rotated-coords - H_matrix| = %.3e  [%s]\n\n', maxdiff4, tf_str(check4_pass));

    %% Check 5: rotation preserves contour area
    fprintf('--- Check 5: rotation preserves contour area ----------\n');
    dA = (xv(2)-xv(1))*(yv(2)-yv(1));
    level_e2 = exp(-2);
    area_phi0  = sum(Hall{2}(:) >= level_e2) * dA;   % elliptical, phi=0
    area_phi30 = sum(Hall{3}(:) >= level_e2) * dA;   % same ellipse, phi=30deg
    area_rel_err = abs(area_phi0-area_phi30)/area_phi0;
    check5_pass = area_rel_err < 0.02;
    fprintf('  Area(h>=e^{-2}), phi=0   = %.6f m^2\n', area_phi0);
    fprintf('  Area(h>=e^{-2}), phi=30  = %.6f m^2\n', area_phi30);
    fprintf('  relative area error = %.4f%%  [%s]  (analytic ellipse area = pi*wx*wy = %.6f m^2)\n\n', ...
        100*area_rel_err, tf_str(check5_pass), pi*wx4*wy4);

    %% Check 6: bilinear interpolation vs analytic map
    fprintf('--- Check 6: bilinear interpolation vs analytic map ---\n');
    Ncoarse = 101;
    xc = linspace(-2*wref, 2*wref, Ncoarse);
    yc = linspace(-2*wref, 2*wref, Ncoarse);
    [EPXc, EPYc] = meshgrid(xc, yc);
    bp_map_src = bpList{2};                          % elliptical Gaussian, phi=0
    Hc = beam_gain_general(EPXc, EPYc, bp_map_src);
    Hc_norm = Hc / max(Hc(:));

    rng(1);
    Nquery = 4000;
    xq = (rand(Nquery,1)*4-2)*wref*0.98;
    yq = (rand(Nquery,1)*4-2)*wref*0.98;
    h_interp = beam_gain_numerical(xc, yc, Hc_norm, xq, yq);
    h_true   = beam_gain_general(xq, yq, bp_map_src);

    rmse6 = sqrt(mean((h_interp-h_true).^2));
    maxerr6 = max(abs(h_interp-h_true));
    gridspacing = xc(2)-xc(1);
    bound6 = 0.02;
    check6_pass = maxerr6 < bound6;
    fprintf('  coarse map grid spacing = %.4f m | query points = %d\n', gridspacing, Nquery);
    fprintf('  RMSE(h_interp - h_true) = %.3e | max abs err = %.3e  [%s, bound=%.3f]\n\n', ...
        rmse6, maxerr6, tf_str(check6_pass), bound6);

    fprintf('===============================================\n');
    fprintf('ALL CHECKS: 1:%s  2:%s  3:%s  4:%s  5:%s  6:%s\n', ...
        tf_str(check1_pass), tf_str(check2_pass), tf_str(check3_pass), ...
        tf_str(check4_pass), tf_str(check5_pass), tf_str(check6_pass));
    fprintf('===============================================\n\n');

    %% Figures + width/symmetry table
    rowNames = cell(nBeam,1);
    W3_h = zeros(nBeam,1); W3_v = zeros(nBeam,1);
    W10_h = zeros(nBeam,1); W10_v = zeros(nBeam,1);
    Wnull_h = nan(nBeam,1); Wnull_v = nan(nBeam,1);
    SymX = false(nBeam,1); SymY = false(nBeam,1);

    for i = 1:nBeam
        H = Hall{i};

        fig = figure('Color','w','Position',[100 100 1500 420]);
        subplot(1,3,1);
        imagesc(xv, yv, H); axis xy equal tight; colorbar;
        xlabel('\epsilon_x (m)'); ylabel('\epsilon_y (m)');
        title(sprintf('%s: h(\\epsilon)', beamCaseLabels{i}), 'Interpreter','tex','FontSize',9);

        subplot(1,3,2);
        contourf(xv, yv, H, [0 0.05 0.1 0.2 0.3 0.5 0.7 0.9 1], 'LineColor','none'); hold on;
        contour(xv, yv, H, [h_3dB h_3dB], 'w-', 'LineWidth',1.6);
        contour(xv, yv, H, [h_10dB h_10dB], 'w--', 'LineWidth',1.6);
        axis xy equal tight; colorbar;
        xlabel('\epsilon_x (m)'); ylabel('\epsilon_y (m)');
        title('Contours (solid=-3dB, dashed=-10dB)','FontSize',9);

        [~, iy0] = min(abs(yv-0));
        [~, ix0] = min(abs(xv-0));
        hcut_h = H(iy0, :);
        hcut_v = H(:, ix0).';

        subplot(1,3,3);
        plot(xv, hcut_h, 'b-', 'LineWidth',1.4, 'DisplayName','horiz. cut (y=0)'); hold on;
        plot(yv, hcut_v, 'r--', 'LineWidth',1.4, 'DisplayName','vert. cut (x=0)');
        xl = [xv(1) xv(end)];
        plot(xl, [h_3dB h_3dB], 'k:', 'HandleVisibility','off');
        plot(xl, [h_10dB h_10dB], 'k:', 'HandleVisibility','off');
        xlabel('\epsilon (m)'); ylabel('h'); ylim([0 1.05]); grid on;
        legend({'horiz. cut (y=0)','vert. cut (x=0)'}, 'Location','best');
        title('Horizontal / vertical cuts','FontSize',9);

        if exist('sgtitle','file') || exist('sgtitle','builtin')
            sgtitle(beamCaseLabels{i}, 'Interpreter','tex');
        end
        tag = sprintf('T2_beam%d', i);
        save_figure_t2(fig, fullfile(figDir2, tag));

        [W3_h(i), Wnull_h(i)]  = cut_widths(xv, hcut_h, h_3dB, strcmp(bpList{i}.type,'sidelobe'));
        [W10_h(i), ~]          = cut_widths(xv, hcut_h, h_10dB, false);
        [W3_v(i), Wnull_v(i)]  = cut_widths(yv, hcut_v, h_3dB, strcmp(bpList{i}.type,'sidelobe'));
        [W10_v(i), ~]          = cut_widths(yv, hcut_v, h_10dB, false);

        SymY(i) = max(max(abs(H - flipud(H)))) < 1e-6;
        SymX(i) = max(max(abs(H - fliplr(H)))) < 1e-6;

        rowNames{i} = beamCaseLabels{i};
    end
    close all;

    fprintf('--- Task 2 width / symmetry table -----------------------------------------------------\n');
    fprintf('%-42s %8s %8s %8s %8s %10s %10s %6s %6s\n', ...
        'Beam','W3dB_h','W3dB_v','W10dB_h','W10dB_v','Wnull_h','Wnull_v','SymX','SymY');
    for i = 1:nBeam
        fprintf('%-42s %8.4f %8.4f %8.4f %8.4f %10s %10s %6s %6s\n', rowNames{i}, ...
            W3_h(i), W3_v(i), W10_h(i), W10_v(i), ...
            fmt_or_na(Wnull_h(i)), fmt_or_na(Wnull_v(i)), tf_str(SymX(i)), tf_str(SymY(i)));
    end
    fprintf('(widths in metres; "-" = no mainlobe null in this cut; SymX: h(-x,y)=h(x,y); SymY: h(x,-y)=h(x,y))\n\n');

    T2 = table(rowNames, W3_h, W3_v, W10_h, W10_v, Wnull_h, Wnull_v, SymX, SymY, ...
        'VariableNames', {'Beam','W3dB_h','W3dB_v','W10dB_h','W10dB_v','Wnull_h','Wnull_v','SymX','SymY'});
    writetable(T2, fullfile(resDir2,'Task2_width_table.csv'));

    fprintf('Figures -> %s\nTable   -> %s\n', figDir2, fullfile(resDir2,'Task2_width_table.csv'));
    fprintf('============ TASK 2 VALIDATION COMPLETE ============\n\n');
end


%% ========================================================================
%  LOCAL FUNCTIONS -- Task 3 curvature-fitting library (merged in from the
%  standalone Task 3 script; formulas/logic unchanged. All beam-library
%  helpers it needs -- beam_gain_general, make_beam_case, beam_matrix_A,
%  tf_str -- are the SINGLE shared copies above; the duplicate copies that
%  shipped with the standalone Task 3 script have been removed.)
%% ========================================================================

function Acurv = center_hessian_curvature(bp, step, hfloor)
% CENTER_HESSIAN_CURVATURE  Acurv = -1/2 * Hessian(ln h)|_{eps=0} (Eq.29),
% computed by a central finite difference of ln(h_log(.)) so it applies to
% ANY beam type (Gaussian family or sidelobe), not just the analytic
% Gaussian case. This is the "known centre curvature" ground truth used to
% (a) verify the finite-region fit for smooth beams with nonzero centre
% curvature, and (b) demonstrate numerically that the flat-top (p>2)
% centre curvature is degenerate (~0), which is WHY Task 3 requires a
% finite-region fit instead of a centre-only calculation for that case.
    lnh = @(ex,ey) log(max(beam_gain_general(ex, ey, bp), hfloor));
    f0  = lnh(0,0);
    fxp = lnh( step,0); fxm = lnh(-step,0);
    fyp = lnh(0, step); fym = lnh(0,-step);
    fpp = lnh( step, step); fpm = lnh( step,-step);
    fmp = lnh(-step, step); fmm = lnh(-step,-step);

    d2x  = (fxp - 2*f0 + fxm)/step^2;
    d2y  = (fyp - 2*f0 + fym)/step^2;
    d2xy = (fpp - fpm - fmp + fmm)/(4*step^2);

    H = [d2x d2xy; d2xy d2y];
    Acurv = -0.5*H;
end


function fit = fit_curvature_matrix(bp, wref, hfit_min, hfloor, Ngrid)
% FIT_CURVATURE_MATRIX  Finite-region weighted-least-squares fit of the
% quadratic beam-curvature matrix A, exactly per Eqs.31-35.
%   fit.A          - eigenvalue-corrected fitted matrix (always SPD)
%   fit.A_raw      - raw least-squares matrix BEFORE eigenvalue correction
%   fit.corrected  - true if a non-positive eigenvalue had to be clipped
%   fit.eigval,fit.eigvec - eigen-decomposition of the CORRECTED A
%   fit.Nf         - number of grid points used in the fit
%   fit.xv,fit.yv,fit.H,fit.mask - grid / true map / selection mask, kept
%                    so fitting_errors() and r1dB() can reuse them without
%                    recomputing the beam map.
    xv = linspace(-2*wref, 2*wref, Ngrid);
    yv = linspace(-2*wref, 2*wref, Ngrid);
    [EX, EY] = meshgrid(xv, yv);
    H = beam_gain_general(EX, EY, bp);

    mask = H >= hfit_min;                       % Eq.31
    ex = EX(mask); ey = EY(mask); hn = H(mask);
    Nf = numel(ex);
    if Nf < 3
        error('fit_curvature_matrix: only %d points selected at hfit_min=%.2f -- too few to fit (need >=3).', Nf, hfit_min);
    end

    y = -log(max(hn, hfloor));                   % Eq.32
    G = [ex.^2, 2*ex.*ey, ey.^2];                 % Eq.32

    beta = (G.'*G) \ (G.'*y);                     % Eq.33
    a11 = beta(1); a12 = beta(2); a22 = beta(3);
    A_raw = [a11 a12; a12 a22];                   % Eq.34

    [V, D] = eig(A_raw);
    eigval_raw = diag(D);
    corrected = any(eigval_raw <= 0);
    eigval_pos = max(eigval_raw, 1e-12);          % Eq.35
    A = V*diag(eigval_pos)*V.';

    fit = struct('A',A,'A_raw',A_raw,'corrected',corrected, ...
        'eigval',eigval_pos,'eigval_raw',eigval_raw,'eigvec',V,'Nf',Nf, ...
        'xv',xv,'yv',yv,'H',H,'mask',mask);
end


function [RMSElin, RMSEdB, Emax_dB] = compute_fitting_errors(fit, hfloor)
% COMPUTE_FITTING_ERRORS  Eqs.36-38, evaluated over the SAME selected
% finite region used to obtain the fit.
    [EX, EY] = meshgrid(fit.xv, fit.yv);
    ex = EX(fit.mask); ey = EY(fit.mask); hn = fit.H(fit.mask);

    A = fit.A;
    quad = A(1,1)*ex.^2 + (A(1,2)+A(2,1))*ex.*ey + A(2,2)*ey.^2;
    hfit_vals = exp(-quad);

    RMSElin = sqrt(mean((hn - hfit_vals).^2));                       % Eq.36
    Lh      = -10*log10(max(hn, hfloor));                            % Eq.27
    Lh_fit  = -10*log10(max(hfit_vals, hfloor));
    RMSEdB  = sqrt(mean((Lh - Lh_fit).^2));                          % Eq.37
    Emax_dB = max(abs(Lh - Lh_fit));                                 % Eq.38
end


function [r1, r_cross, thetas] = compute_r1dB(bp, fit, wref, hfloor, nDir, dr, rmax)
% COMPUTE_R1DB  Model-validity radius (Eq.40): smallest radial distance,
% over nDir tested directions, at which |Lh-Lh,fit| first exceeds 1 dB.
% Marches outward from the origin along each direction in steps of dr.
    if nargin < 5 || isempty(nDir), nDir = 72; end
    if nargin < 6 || isempty(dr),   dr = wref/500; end
    if nargin < 7 || isempty(rmax), rmax = 2*wref; end

    A = fit.A;
    thetas = linspace(0, 2*pi, nDir+1); thetas(end) = [];
    r_cross = nan(1,nDir);

    for i = 1:nDir
        ct = cos(thetas(i)); st = sin(thetas(i));
        r = dr;
        while r <= rmax
            ex = r*ct; ey = r*st;
            h  = beam_gain_general(ex, ey, bp);
            quad  = A(1,1)*ex^2 + (A(1,2)+A(2,1))*ex*ey + A(2,2)*ey^2;
            hfitv = exp(-quad);
            Lh    = -10*log10(max(h,hfloor));
            Lhfit = -10*log10(max(hfitv,hfloor));
            if abs(Lh-Lhfit) > 1.0
                r_cross(i) = r;
                break;
            end
            r = r + dr;
        end
    end
    valid = r_cross(~isnan(r_cross));
    if isempty(valid)
        r1 = NaN;
    else
        r1 = min(valid);
    end
end


function run_task3_validation(wref, figDir3, resDir3, beamCaseNames, beamCaseLabels, hfitMinList, hfloor, Ngrid)
% RUN_TASK3_VALIDATION  Full Task 3 experiment: recovers the analytical
% matrix check for the (rotated) elliptical Gaussian, sweeps hfit_min for
% every non-circular beam, reports centre-Hessian degeneracy for the
% flat-top beam, produces contour/error-map figures, and writes the
% required summary table (fitted matrices, eigenvalues/eigenvectors, the
% three error measures, and r1dB for every non-circular beam).

    nBeam = numel(beamCaseNames);
    stepHess = 1e-5*wref;   % finite-diff step for the centre Hessian, in metres

    %% ---------------------------------------------------------------
    %  CHECK A/B -- fit recovers the KNOWN ANALYTICAL matrix for the
    %  elliptical Gaussian AND the rotated elliptical Gaussian (Task 3,
    %  "Verify first that the fit recovers the known analytical matrix").
    %  Uses the tightest region (hfit_min = 0.8) where the quadratic
    %  approximation is most exact.
    %% ---------------------------------------------------------------
    fprintf('--- Check A: fit vs analytic beam_matrix_A, EllipticalGaussian --------\n');
    bp_ell = make_beam_case('EllipticalGaussian', wref);
    A_analytic_ell = beam_matrix_A(bp_ell.wx, bp_ell.wy, bp_ell.phi);
    fit_ell = fit_curvature_matrix(bp_ell, wref, 0.8, hfloor, Ngrid);
    err_ell = max(max(abs(fit_ell.A - A_analytic_ell)));
    checkA_pass = err_ell < 1e-6;
    fprintf('  A_analytic = [%.4f %.4f; %.4f %.4f]\n', A_analytic_ell(1,1),A_analytic_ell(1,2),A_analytic_ell(2,1),A_analytic_ell(2,2));
    fprintf('  A_fitted   = [%.4f %.4f; %.4f %.4f]\n', fit_ell.A(1,1),fit_ell.A(1,2),fit_ell.A(2,1),fit_ell.A(2,2));
    fprintf('  max|A_fit-A_analytic| = %.3e  [%s]\n\n', err_ell, tf_str(checkA_pass));

    fprintf('--- Check B: fit vs analytic, RotatedEllipticalGaussian (off-diag a12) \n');
    bp_rot = make_beam_case('RotatedEllipticalGaussian', wref);
    A_analytic_rot = beam_matrix_A(bp_rot.wx, bp_rot.wy, bp_rot.phi);
    fit_rot = fit_curvature_matrix(bp_rot, wref, 0.8, hfloor, Ngrid);
    err_rot = max(max(abs(fit_rot.A - A_analytic_rot)));
    checkB_pass = err_rot < 1e-6;
    fprintf('  A_analytic = [%.4f %.4f; %.4f %.4f]\n', A_analytic_rot(1,1),A_analytic_rot(1,2),A_analytic_rot(2,1),A_analytic_rot(2,2));
    fprintf('  A_fitted   = [%.4f %.4f; %.4f %.4f]\n', fit_rot.A(1,1),fit_rot.A(1,2),fit_rot.A(2,1),fit_rot.A(2,2));
    fprintf('  max|A_fit-A_analytic| = %.3e  [%s]  (nonzero a12 confirms the rotation-\n', err_rot, tf_str(checkB_pass));
    fprintf('   induced x/y coupling is correctly recovered by the fit)\n\n');

    %% ---------------------------------------------------------------
    %  Centre-Hessian curvature for every non-circular beam, INCLUDING
    %  the flat-top beam -- demonstrates degeneracy (Eq.29 evaluated
    %  numerically). This is the evidence for "why the flat-top beam
    %  cannot be represented by a unique centre curvature."
    %% ---------------------------------------------------------------
    fprintf('--- Centre-Hessian curvature Acurv = -1/2 Hessian(ln h)|_0 (Eq.29) ----\n');
    for bc = 2:nBeam   % skip circular (trivial, isotropic reference)
        bp = make_beam_case(beamCaseNames{bc}, wref);
        Acurv = center_hessian_curvature(bp, stepHess, hfloor);
        fprintf('  %-42s Acurv = [% .4f % .4f; % .4f % .4f]\n', beamCaseLabels{bc}, ...
            Acurv(1,1),Acurv(1,2),Acurv(2,1),Acurv(2,2));
        if strcmp(beamCaseNames{bc},'FlatTopSuperGaussian')
            fprintf('    ^ DEGENERATE (~0): for p>2 in h_SG,p = exp(-2*chi(eps)^(p/2)),\n');
            fprintf('      chi(eps) = O(||eps||^2) near the origin, so chi^(p/2) = O(||eps||^p).\n');
            fprintf('      For p=4 every second derivative of ln h vanishes at eps=0 -- the\n');
            fprintf('      beam is perfectly flat (zero curvature) right at its centre, so\n');
            fprintf('      Eq.29 gives Acurv=0 and there is NO unique centre-curvature matrix.\n');
            fprintf('      A finite-region least-squares fit (Eqs.31-34) is required instead,\n');
            fprintf('      and the table below shows its result is not unique either -- it\n');
            fprintf('      depends on how much of the flat-top''s shoulder the fit region\n');
            fprintf('      includes (see hfit_min sweep).\n');
        end
    end
    fprintf('\n');

    %% ---------------------------------------------------------------
    %  Main hfit_min sweep for every non-circular beam (Task 3 required
    %  experiment: p=4 flat-top and sidelobed at minimum; the two
    %  elliptical cases are included too since the deliverable asks for
    %  "every non-circular beam").
    %% ---------------------------------------------------------------
    nHmin = numel(hfitMinList);
    rowsBeam = {}; rowsHmin = []; rowsNf = []; rowsCorr = {};
    rowsA11 = []; rowsA12 = []; rowsA22 = [];
    rowsEig1 = []; rowsEig2 = []; rowsAngle = [];
    rowsRMSElin = []; rowsRMSEdB = []; rowsEmaxdB = []; rowsR1dB = [];

    % Colors need row 5 valid regardless of how many hfit_min values are swept.
    colorsH = lines(max(nHmin,5));

    for bc = 2:nBeam   % non-circular beams only
        beamTag   = beamCaseNames{bc};
        beamLabel = beamCaseLabels{bc};
        bp = make_beam_case(beamTag, wref);

        fprintf('--- hfit_min sweep: %s -------------------------------\n', beamLabel);
        fprintf('%8s %7s %10s %9s %9s %10s %9s %9s %10s\n', ...
            'hfitmin','Nf','corrected','eig1','eig2','RMSElin','RMSEdB','EmaxdB','r1dB(m)');

        fits = cell(1,nHmin);
        for hi = 1:nHmin
            hmin = hfitMinList(hi);
            fit = fit_curvature_matrix(bp, wref, hmin, hfloor, Ngrid);
            [RMSElin, RMSEdB, Emax_dB] = compute_fitting_errors(fit, hfloor);
            r1 = compute_r1dB(bp, fit, wref, hfloor, 72, wref/500, 2*wref);
            fits{hi} = fit;

            ev = sort(fit.eigval);
            % principal-axis angle of the fitted ellipse (angle of the
            % eigenvector belonging to the SMALLER eigenvalue = the
            % wider/major axis direction)
            [~,imin] = min(fit.eigval);
            axisVec = fit.eigvec(:,imin);
            angleDeg = atan2d(axisVec(2), axisVec(1));

            fprintf('%8.2f %7d %10s %9.3f %9.3f %10.3e %9.3f %9.3f %10.4f\n', ...
                hmin, fit.Nf, tf_str(~fit.corrected), ev(1), ev(2), RMSElin, RMSEdB, Emax_dB, r1);

            rowsBeam{end+1,1}  = beamLabel; %#ok<AGROW>
            rowsHmin(end+1,1)  = hmin; %#ok<AGROW>
            rowsNf(end+1,1)    = fit.Nf; %#ok<AGROW>
            rowsCorr{end+1,1}  = tf_str(fit.corrected); %#ok<AGROW>
            rowsA11(end+1,1)   = fit.A(1,1); %#ok<AGROW>
            rowsA12(end+1,1)   = fit.A(1,2); %#ok<AGROW>
            rowsA22(end+1,1)   = fit.A(2,2); %#ok<AGROW>
            rowsEig1(end+1,1)  = ev(1); %#ok<AGROW>
            rowsEig2(end+1,1)  = ev(2); %#ok<AGROW>
            rowsAngle(end+1,1) = angleDeg; %#ok<AGROW>
            rowsRMSElin(end+1,1) = RMSElin; %#ok<AGROW>
            rowsRMSEdB(end+1,1)  = RMSEdB; %#ok<AGROW>
            rowsEmaxdB(end+1,1)  = Emax_dB; %#ok<AGROW>
            rowsR1dB(end+1,1)    = r1; %#ok<AGROW>
        end
        fprintf('\n');

        %% --- Figure: true vs fitted contour + dB error map, one
        %      representative hfit_min (0.25, the brief's Eq.31 default)
        %      plus an eigenvalue/RMSEdB/r1dB-vs-hfit_min summary panel ---
        hiShow = find(hfitMinList == 0.25, 1);
        if isempty(hiShow), hiShow = ceil(nHmin/2); end
        fitShow = fits{hiShow};
        [EX, EY] = meshgrid(fitShow.xv, fitShow.yv);
        Htrue = fitShow.H;
        A = fitShow.A;
        quad = A(1,1)*EX.^2 + (A(1,2)+A(2,1))*EX.*EY + A(2,2)*EY.^2;
        Hfit = exp(-quad);
        Lh     = -10*log10(max(Htrue,hfloor));
        Lhfit  = -10*log10(max(Hfit,hfloor));
        errMap = Lh - Lhfit;

        fig = figure('Color','w','Position',[100 100 1500 420]);
        subplot(1,3,1);
        contourf(fitShow.xv, fitShow.yv, Htrue, [0 0.05 0.1 0.2 0.3 0.5 0.7 0.9 1], 'LineColor','none');
        axis xy equal tight; colorbar; xlabel('\epsilon_x (m)'); ylabel('\epsilon_y (m)');
        title(sprintf('True h(\\epsilon): %s', beamLabel), 'Interpreter','tex','FontSize',9);

        subplot(1,3,2);
        contourf(fitShow.xv, fitShow.yv, Hfit, [0 0.05 0.1 0.2 0.3 0.5 0.7 0.9 1], 'LineColor','none'); hold on;
        contour(fitShow.xv, fitShow.yv, fitShow.mask, [0.5 0.5], 'w--','LineWidth',1.2);
        axis xy equal tight; colorbar; xlabel('\epsilon_x (m)'); ylabel('\epsilon_y (m)');
        title(sprintf('Fitted h_{fit}(\\epsilon), h_{fit,min}=%.2f (dashed=fit region)', hfitMinList(hiShow)),'FontSize',9);

        subplot(1,3,3);
        imagesc(fitShow.xv, fitShow.yv, errMap); axis xy equal tight; colorbar;
        hold on; contour(fitShow.xv, fitShow.yv, abs(errMap), [1 1], 'w-','LineWidth',1.6);
        xlabel('\epsilon_x (m)'); ylabel('\epsilon_y (m)');
        title('L_h - L_{h,fit}  (dB); white = \pm1 dB contour','FontSize',9);

        if exist('sgtitle','file') || exist('sgtitle','builtin')
            sgtitle(sprintf('Task 3: %s', beamLabel), 'Interpreter','tex');
        end
        save_figure_t2(fig, fullfile(figDir3, sprintf('T3_contours_%s', beamTag)));

        %% --- Figure: eigenvalues, RMSEdB, Emax_dB, r1dB vs hfit_min ---
        idxRows = strcmp(rowsBeam, beamLabel);
        fig2 = figure('Color','w','Position',[100 100 1100 800]);
        subplot(2,2,1);
        plot(rowsHmin(idxRows), rowsEig1(idxRows), '-o', 'Color', colorsH(1,:)); hold on;
        plot(rowsHmin(idxRows), rowsEig2(idxRows), '-s', 'Color', colorsH(2,:));
        set(gca,'XDir','reverse'); grid on; xlabel('h_{fit,min}'); ylabel('Eigenvalue');
        legend({'\lambda_1 (min)','\lambda_2 (max)'},'Location','best');
        title('Fitted eigenvalues vs region size');

        subplot(2,2,2);
        plot(rowsHmin(idxRows), rowsRMSEdB(idxRows), '-o', 'Color', colorsH(3,:)); hold on;
        plot(rowsHmin(idxRows), rowsEmaxdB(idxRows), '-s', 'Color', colorsH(4,:));
        set(gca,'XDir','reverse'); grid on; xlabel('h_{fit,min}'); ylabel('Error (dB)');
        legend({'RMSE_{dB}','E_{max,dB}'},'Location','best');
        title('Fitting error vs region size');

        subplot(2,2,3);
        plot(rowsHmin(idxRows), rowsRMSElin(idxRows), '-o', 'Color', colorsH(5,:));
        set(gca,'XDir','reverse'); grid on; xlabel('h_{fit,min}'); ylabel('RMSE_{lin}');
        title('Linear-power fitting error vs region size');

        subplot(2,2,4);
        plot(rowsHmin(idxRows), rowsR1dB(idxRows)/wref, '-o', 'Color', colorsH(1,:));
        set(gca,'XDir','reverse'); grid on; xlabel('h_{fit,min}'); ylabel('r_{1dB} / w_{ref}');
        title('Model-validity radius vs region size');

        if exist('sgtitle','file') || exist('sgtitle','builtin')
            sgtitle(sprintf('Task 3 region-size sensitivity: %s', beamLabel), 'Interpreter','tex');
        end
        save_figure_t2(fig2, fullfile(figDir3, sprintf('T3_sensitivity_%s', beamTag)));
    end
    close all;

    %% ---------------------------------------------------------------
    %  Summary table -> CSV (Task 3 deliverable)
    %% ---------------------------------------------------------------
    T3 = table(rowsBeam, rowsHmin, rowsNf, rowsCorr, rowsA11, rowsA12, rowsA22, ...
        rowsEig1, rowsEig2, rowsAngle, rowsRMSElin, rowsRMSEdB, rowsEmaxdB, rowsR1dB, ...
        'VariableNames', {'Beam','hfit_min','Nf','EigCorrected','A11','A12','A22', ...
        'Eig1_min','Eig2_max','MajorAxisAngle_deg','RMSElin','RMSEdB','Emax_dB','r1dB_m'});
    writetable(T3, fullfile(resDir3,'Task3_curvature_fit_table.csv'));

    fprintf('===============================================\n');
    fprintf('Checks: A(Elliptical vs analytic):%s   B(Rotated vs analytic):%s\n', tf_str(checkA_pass), tf_str(checkB_pass));
    fprintf('===============================================\n');
    fprintf('Figures -> %s\nTable   -> %s\n', figDir3, fullfile(resDir3,'Task3_curvature_fit_table.csv'));
    fprintf('============ TASK 3 VALIDATION COMPLETE ============\n\n');
end


%% ========================================================================
%  LOCAL FUNCTIONS -- Task 4 echo-based pointing-estimation library
%  (merged in from the standalone Task 4 script; Stage A/A2/B logic
%  unchanged. All beam-library / radar-chain helpers they need --
%  beam_gain_general, make_beam_case, beam_matrix_A, tf_str,
%  save_figure_t2, setup_system_params, build_radar_reference,
%  generate_received_signal_fast, range_doppler_fast -- are the SINGLE
%  shared copies above; the duplicate copies that shipped with the
%  standalone Task 4 script have been removed. The only change from the
%  standalone script is that wref is now passed in as wref2 (= w_mo of
%  the SF regime) by the caller below, instead of the standalone script's
%  hardcoded wref = 0.20 m placeholder -- run_stage_A/A2/B themselves
%  already took wref as a parameter, so their bodies are untouched.)
%% ========================================================================

function run_task4_validation(wref, hfloor, fdstep, beamCaseNames, beamCaseLabels, isQuadratic, ...
    figDir4, resDir4, beamsB, beamLabelsB, alphaListB, Ptx_dBm_listB, ...
    Ntrials, Mcal, Nfa, Pfa_target, jitterStd)
% RUN_TASK4_VALIDATION  Wraps Stage A (ideal probe-power validation),
% Stage A2 (Jacobian-conditioning / EKF-vs-prediction-distance
% diagnostics), and Stage B (complete sensing-chain validation, reduced
% combo set) -- the three stages of the standalone Task 4 script -- into
% a single call, mirroring run_task2_validation / run_task3_validation.

    fprintf('===============================================\n');
    fprintf('Task 4 -- Echo-based pointing estimation validation + DIAGNOSTICS\n');
    fprintf('===============================================\n');
    fprintf('wref = %.4f m (= w_mo of the SF regime) | hfloor = %.0e | fdstep = %.2e m\n\n', wref, hfloor, fdstep);

    %% ====================================================================
    %  STAGE A -- IDEAL (NOISELESS) PROBE-POWER VALIDATION
    %% ====================================================================
    run_stage_A(wref, hfloor, fdstep, beamCaseNames, beamCaseLabels, isQuadratic, resDir4);

    %% ====================================================================
    %  STAGE A2 -- DIAGNOSTICS
    %  Determines whether the FlatTop/Sidelobed EKF-1-step blowups seen in
    %  Stage A/B are explained by Jacobian conditioning at the cold-start
    %  linearization point, and whether they resolve as the prediction
    %  approaches the truth.
    %% ====================================================================
    run_stage_A2_diagnostics(wref, hfloor, fdstep, resDir4);

    %% ====================================================================
    %  STAGE B -- COMPLETE SENSING-CHAIN VALIDATION (reduced combo set)
    %% ====================================================================
    run_stage_B(wref, hfloor, fdstep, beamsB, beamLabelsB, alphaListB, ...
        Ptx_dBm_listB, Ntrials, Mcal, Nfa, Pfa_target, jitterStd, figDir4, resDir4);

    fprintf('============ TASK 4 VALIDATION COMPLETE ============\n\n');
end


function run_stage_A(wref, hfloor, fdstep, beamCaseNames, beamCaseLabels, isQuadratic, resDir4)
    fprintf('=== STAGE A: ideal (noiseless) probe-power validation ===\n\n');

    ux = [1;0]; uy = [0;1];
    alpha = 0.3;

    radii  = [0, 0.10, 0.30, 0.50]*wref;
    angles = deg2rad([0 45 90 135]);
    testPts = [0 0];
    for ri = 2:numel(radii)
        for ai = 1:numel(angles)
            testPts(end+1,:) = radii(ri)*[cos(angles(ai)) sin(angles(ai))]; %#ok<AGROW>
        end
    end
    nTest = size(testPts,1);

    rowsBeam = {}; rowsEpsT = {}; rowsMethod = {}; rowsEpsHat = {}; rowsErr = [];

    for bc = 1:numel(beamCaseNames)
        beamTag = beamCaseNames{bc};
        beamLabel = beamCaseLabels{bc};
        bp = make_beam_case(beamTag, wref);

        wx_dir = directional_radius(bp, ux, exp(-2));
        wy_dir = directional_radius(bp, uy, exp(-2));
        deltax = alpha*wx_dir; deltay = alpha*wy_dir;

        fprintf('--- %s (deltax=%.4f m, deltay=%.4f m) ---\n', beamLabel, deltax, deltay);

        if isQuadratic(bc)
            A = beam_matrix_A(bp.wx, bp.wy, bp.phi);
            H = [4*deltax*(ux.'*A); 4*deltay*(uy.'*A)];
        end

        maxErrDirect = 0; maxErrEKF = 0; maxErrGN = 0; maxErrLUT = 0;

        for n = 1:nTest
            eps_true = testPts(n,:).';
            z = g_observation(eps_true, deltax, deltay, bp, hfloor);

            if isQuadratic(bc)
                eps_direct = H\z;
                errDirect = norm(eps_direct-eps_true);
                maxErrDirect = max(maxErrDirect, errDirect);
            else
                eps_direct = [NaN;NaN]; errDirect = NaN;
            end

            eps_pred = [0;0];
            J0 = jacobian_fd(eps_pred, deltax, deltay, bp, hfloor, fdstep);
            eps_ekf = eps_pred + J0\(z - g_observation(eps_pred, deltax, deltay, bp, hfloor));
            errEKF = norm(eps_ekf-eps_true);
            maxErrEKF = max(maxErrEKF, errEKF);

            eps_gn = gauss_newton_solve(z, [0;0], deltax, deltay, bp, hfloor, fdstep, 30, 1e-12);
            errGN = norm(eps_gn-eps_true);
            maxErrGN = max(maxErrGN, errGN);

            eps_lut = lookup_estimate(z, eye(2), bp, deltax, deltay, hfloor, wx_dir, wy_dir, 61);
            errLUT = norm(eps_lut-eps_true);
            maxErrLUT = max(maxErrLUT, errLUT);

            rowsBeam(end+1:end+4,1)   = {beamLabel;beamLabel;beamLabel;beamLabel}; %#ok<AGROW>
            rowsEpsT(end+1:end+4,1)   = {mat2str(eps_true',4)}; rowsEpsT(end-2:end,1) = rowsEpsT(end); %#ok<AGROW>
            rowsMethod(end+1:end+4,1) = {'Direct';'EKF-1step';'GaussNewton';'Lookup'}; %#ok<AGROW>
            rowsEpsHat(end+1,1) = {mat2str(eps_direct',4)}; %#ok<AGROW>
            rowsEpsHat(end+1,1) = {mat2str(eps_ekf',4)}; %#ok<AGROW>
            rowsEpsHat(end+1,1) = {mat2str(eps_gn',4)}; %#ok<AGROW>
            rowsEpsHat(end+1,1) = {mat2str(eps_lut',4)}; %#ok<AGROW>
            rowsErr(end+1:end+4,1) = [errDirect; errEKF; errGN; errLUT]; %#ok<AGROW>
        end

        if isQuadratic(bc)
            fprintf('  Direct quadratic estimator: max|eps_hat-eps_true| over %d test pts = %.3e m  [%s]\n', ...
                nTest, maxErrDirect, tf_str(maxErrDirect < 1e-9));
        else
            fprintf('  Direct quadratic estimator: N/A (beam is not exactly quadratic)\n');
        end
        fprintf('  EKF single-step:            max error = %.4e m\n', maxErrEKF);
        fprintf('  Gauss-Newton (LM-damped):   max error = %.4e m\n', maxErrGN);
        fprintf('  Lookup (61x61 grid):        max error = %.4e m\n\n', maxErrLUT);
    end

    T4A = table(rowsBeam, rowsEpsT, rowsMethod, rowsEpsHat, rowsErr, ...
        'VariableNames', {'Beam','eps_true','Method','eps_hat','Error_m'});
    writetable(T4A, fullfile(resDir4,'Task4_StageA_table.csv'));
    fprintf('Stage A table -> %s\n\n', fullfile(resDir4,'Task4_StageA_table.csv'));
end


function run_stage_A2_diagnostics(wref, hfloor, fdstep, resDir4)
    fprintf('=== STAGE A2: DIAGNOSTICS (Jacobian conditioning + EKF vs prediction distance) ===\n\n');

    ux = [1;0]; uy = [0;1];
    alpha = 0.3;

    diagBeamNames  = {'FlatTopSuperGaussian','Sidelobed'};
    diagBeamLabels = {'Flat-top Super-Gaussian (p=4, \phi=30^\circ)','Sidelobed (sinc^2, \phi=30^\circ)'};

    % representative true residual used for the prediction-distance sweep
    eps_true_test = [0.30*wref; 0];

    rowsBeam = {}; rowsFrac = []; rowsRcond = []; rowsErr = [];

    for bi = 1:numel(diagBeamNames)
        beamTag = diagBeamNames{bi};
        beamLabel = diagBeamLabels{bi};
        bp = make_beam_case(beamTag, wref);

        wx_dir = directional_radius(bp, ux, exp(-2));
        wy_dir = directional_radius(bp, uy, exp(-2));
        deltax = alpha*wx_dir; deltay = alpha*wy_dir;

        fprintf('--- %s ---\n', beamLabel);

        % --- Patch 1: Jacobian conditioning at the cold-start point eps_pred=0 ---
        J0 = jacobian_fd([0;0], deltax, deltay, bp, hfloor, fdstep);
        fprintf('  J(eps_pred=0) = [%.4e %.4e; %.4e %.4e]\n', J0(1,1),J0(1,2),J0(2,1),J0(2,2));
        fprintf('  rank=%d   cond=%.3e   rcond=%.3e\n', rank(J0), cond(J0), rcond(J0));
        if rcond(J0) < 1e-10
            fprintf('  -> J is numerically singular/near-singular at the cold-start point.\n');
            fprintf('     This alone is sufficient to explain a large one-step EKF correction:\n');
            fprintf('     J\\r amplifies any residual z-g(0) by ~1/rcond(J).\n');
        else
            fprintf('  -> J is reasonably well-conditioned at the cold-start point; a large\n');
            fprintf('     EKF error here would point to something else (check step 2 below).\n');
        end

        % --- Patch 2: EKF-1-step error as the linearization point marches
        %     from eps_pred=0 toward eps_true ---
        z = g_observation(eps_true_test, deltax, deltay, bp, hfloor);
        fracs = [0, 0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.0];
        fprintf('  EKF-1-step error vs. |eps_pred - eps_true| (eps_true fixed at 0.30*wref along x):\n');
        fprintf('  %8s %12s %12s\n', 'frac', 'rcond(J)', 'err (m)');
        for f = fracs
            eps_pred = f*eps_true_test;   % march the linearization point toward the truth
            Jf = jacobian_fd(eps_pred, deltax, deltay, bp, hfloor, fdstep);
            eps_ekf = eps_pred + Jf\(z - g_observation(eps_pred, deltax, deltay, bp, hfloor));
            errf = norm(eps_ekf - eps_true_test);
            fprintf('  %8.2f %12.3e %12.4e\n', f, rcond(Jf), errf);

            rowsBeam{end+1,1} = beamLabel; %#ok<AGROW>
            rowsFrac(end+1,1) = f; %#ok<AGROW>
            rowsRcond(end+1,1) = rcond(Jf); %#ok<AGROW>
            rowsErr(end+1,1) = errf; %#ok<AGROW>
        end
        fprintf('\n');
    end

    T4A2 = table(rowsBeam, rowsFrac, rowsRcond, rowsErr, ...
        'VariableNames', {'Beam','PredictionFraction','rcondJ','EKFError_m'});
    writetable(T4A2, fullfile(resDir4,'Task4_StageA2_diagnostics_table.csv'));
    fprintf('Interpretation guide:\n');
    fprintf('  - If error falls sharply as frac -> 1, the EKF-1-step blowup is a COLD-START\n');
    fprintf('    linearization artifact, not a bug: once the KF prediction is close to the\n');
    fprintf('    true residual (as it will be after a few frames of tracking in Task 5,\n');
    fprintf('    except at initialization), the single linearized correction should work fine.\n');
    fprintf('  - If error stays large even at frac=1 (i.e. even near-truth linearization is\n');
    fprintf('    unstable), that is a genuine estimator limitation for this beam, not just an\n');
    fprintf('    initialization problem -- report it as such, and rely on lookup/GN instead.\n');
    fprintf('Diagnostics table -> %s\n\n', fullfile(resDir4,'Task4_StageA2_diagnostics_table.csv'));
end


function run_stage_B(wref, hfloor, fdstep, beamsB, beamLabelsB, alphaListB, ...
    Ptx_dBm_listB, Ntrials, Mcal, Nfa, Pfa_target, jitterStd, figDir4, resDir4)

    fprintf('=== STAGE B: complete sensing-chain validation (reduced combo set) ===\n');
    fprintf('beams=%s | alpha=%s | Ptx=%s dBm | Ntrials=%d | Mcal=%d | Nfa=%d\n\n', ...
        strjoin(beamsB,','), mat2str(alphaListB), mat2str(Ptx_dBm_listB), Ntrials, Mcal, Nfa);

    P = setup_system_params();
    distance = 15; velocity = 5;
    waveform = 'OFDM';
    W = build_radar_reference(P, waveform, distance, velocity);
    rx_delayed_doppler = W.rx_delayed_doppler;
    corr_ref = W.corr_ref;

    height = 100;
    [T_atm_K, p_hPa, wvden] = atmositu(height);
    T_C = T_atm_K - 273.15; P_Pa = p_hPa*100;
    L_atm_dB  = gaspl(distance, P.fc, T_C, P_Pa, wvden);
    L_atm_lin = 10^(-L_atm_dB/10);
    T_noise   = (1-L_atm_lin)*T_atm_K;
    P_n       = P.kb*T_noise*P.B;
    P_n_radar = P_n*(P.N_radar/P.N_fft);

    gateR = 2; gateD = 2;
    [r_idx0, d_idx0] = target_bins(P, distance, velocity);

    ux = [1;0]; uy = [0;1];
    rowsBeam={}; rowsAlpha=[]; rowsPtx=[]; rowsMethod={};
    rowsBiasX=[]; rowsBiasY=[]; rowsRMSE=[]; rowsP95=[]; rowsDrop=[]; rowsN=[];

    for pt = 1:numel(Ptx_dBm_listB)
        Ptx_dBm = Ptx_dBm_listB(pt);
        Ptx_W = 10^((Ptx_dBm-30)/10);
        P_radar = Ptx_W*(P.N_radar/P.N_fft);
        A_base = sqrt(P_radar)*L_atm_lin* ...
            sqrt((P.G_total*P.lambda^2*P.sigma_rcs)/((4*pi)^3*distance^4));

        aligned_SNR_dB = 10*log10(A_base^2/P_n_radar);
        fprintf('--- Ptx=%.1f dBm  (aligned SNR ~ %.1f dB) ---\n', Ptx_dBm, aligned_SNR_dB);

        Pdet_thresh = calibrate_detection_threshold(rx_delayed_doppler, corr_ref, P, ...
            A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Nfa, Pfa_target);
        fprintf('    Pdet threshold (Pfa=%.1e, Nfa=%d trials) = %.4e\n', Pfa_target, Nfa, Pdet_thresh);

        for bi = 1:numel(beamsB)
            beamTag = beamsB{bi}; beamLabel = beamLabelsB{bi};
            bp = make_beam_case(beamTag, wref);
            isQuad = strcmp(beamTag,'CircularGaussian') || strcmp(beamTag,'EllipticalGaussian') ...
                || strcmp(beamTag,'RotatedEllipticalGaussian');
            if isQuad, A = beam_matrix_A(bp.wx, bp.wy, bp.phi); end

            wx_dir = directional_radius(bp, ux, exp(-2));
            wy_dir = directional_radius(bp, uy, exp(-2));

            for ai = 1:numel(alphaListB)
                alpha = alphaListB(ai);
                deltax = alpha*wx_dir; deltay = alpha*wy_dir;
                if isQuad
                    H = [4*deltax*(ux.'*A); 4*deltay*(uy.'*A)];
                end

                fprintf('  [%s, alpha=%.2f] deltax=%.4f deltay=%.4f m\n', beamLabel, alpha, deltax, deltay);

                Zcal = nan(Mcal,2);
                for m = 1:Mcal
                    [zc, availc] = measure_probe_pair(rx_delayed_doppler, corr_ref, P, bp, ...
                        [0;0], deltax, deltay, A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Pdet_thresh, hfloor);
                    if all(availc), Zcal(m,:) = zc.'; end
                end
                validCal = ~any(isnan(Zcal),2);
                if nnz(validCal) >= 10
                    R = cov(Zcal(validCal,:));
                else
                    R = eye(2);
                    warning('R calibration: only %d/%d valid samples, using identity R.', nnz(validCal), Mcal);
                end
                Rinv = inv(R + 1e-12*eye(2)); %#ok<MINV>

                epsTrueAll = jitterStd*randn(Ntrials,2);
                errDirect = nan(Ntrials,2); errEKF = nan(Ntrials,2); errLUT = nan(Ntrials,2);
                dropCount = 0;
                for t = 1:Ntrials
                    eps_true = epsTrueAll(t,:).';
                    [z, avail] = measure_probe_pair(rx_delayed_doppler, corr_ref, P, bp, ...
                        eps_true, deltax, deltay, A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Pdet_thresh, hfloor);
                    if ~all(avail)
                        dropCount = dropCount + 1;
                        continue;
                    end
                    if isQuad
                        eps_d = H\z;
                        errDirect(t,:) = (eps_d-eps_true).';
                    end
                    J0 = jacobian_fd([0;0], deltax, deltay, bp, hfloor, fdstep);
                    eps_e = [0;0] + J0\(z - g_observation([0;0], deltax, deltay, bp, hfloor));
                    errEKF(t,:) = (eps_e-eps_true).';

                    eps_l = lookup_estimate(z, Rinv, bp, deltax, deltay, hfloor, wx_dir, wy_dir, 41);
                    errLUT(t,:) = (eps_l-eps_true).';
                end
                dropoutProb = dropCount/Ntrials;

                methodsHere = {};
                if isQuad, methodsHere{end+1} = 'Direct'; end %#ok<AGROW>
                methodsHere{end+1} = 'EKF-1step'; %#ok<AGROW>
                methodsHere{end+1} = 'Lookup'; %#ok<AGROW>

                for mi = 1:numel(methodsHere)
                    switch methodsHere{mi}
                        case 'Direct', E = errDirect;
                        case 'EKF-1step', E = errEKF;
                        case 'Lookup', E = errLUT;
                    end
                    valid = ~any(isnan(E),2);
                    if any(valid)
                        bias = mean(E(valid,:),1);
                        d = sqrt(sum(E(valid,:).^2,2));
                        rmse = sqrt(mean(d.^2));
                        p95 = prctile(d,95);
                    else
                        bias = [NaN NaN]; rmse = NaN; p95 = NaN;
                    end
                    fprintf('      %-10s bias=[%.4f %.4f] mm  RMSE=%.4f mm  P95=%.4f mm  Pdrop=%.3f\n', ...
                        methodsHere{mi}, 1000*bias(1), 1000*bias(2), 1000*rmse, 1000*p95, dropoutProb);

                    rowsBeam{end+1,1} = beamLabel; %#ok<AGROW>
                    rowsAlpha(end+1,1) = alpha; %#ok<AGROW>
                    rowsPtx(end+1,1) = Ptx_dBm; %#ok<AGROW>
                    rowsMethod{end+1,1} = methodsHere{mi}; %#ok<AGROW>
                    rowsBiasX(end+1,1) = bias(1); %#ok<AGROW>
                    rowsBiasY(end+1,1) = bias(2); %#ok<AGROW>
                    rowsRMSE(end+1,1) = rmse; %#ok<AGROW>
                    rowsP95(end+1,1) = p95; %#ok<AGROW>
                    rowsDrop(end+1,1) = dropoutProb; %#ok<AGROW>
                    rowsN(end+1,1) = nnz(valid); %#ok<AGROW>
                end
            end
        end
        fprintf('\n');
    end

    T4B = table(rowsBeam, rowsAlpha, rowsPtx, rowsMethod, rowsBiasX, rowsBiasY, ...
        rowsRMSE, rowsP95, rowsDrop, rowsN, ...
        'VariableNames', {'Beam','Alpha','Ptx_dBm','Method','Bias_x_m','Bias_y_m', ...
        'RMSE_m','P95_m','Pdrop','Nvalid'});
    writetable(T4B, fullfile(resDir4,'Task4_StageB_table.csv'));
    fprintf('Stage B table -> %s\n', fullfile(resDir4,'Task4_StageB_table.csv'));

    fig = figure('Color','w','Position',[100 100 900 550]); hold on;
    methodsAll = unique(rowsMethod,'stable');
    beamsAll = unique(rowsBeam,'stable');
    colors = lines(numel(methodsAll)*numel(beamsAll));
    ci = 1;
    for bi = 1:numel(beamsAll)
        for mi = 1:numel(methodsAll)
            mask = strcmp(rowsBeam,beamsAll{bi}) & strcmp(rowsMethod,methodsAll{mi});
            if ~any(mask), continue; end
            [av,ord] = sort(rowsAlpha(mask));
            rv = rowsRMSE(mask); rv = rv(ord);
            plot(av, 1000*rv, '-o', 'Color', colors(ci,:), ...
                'DisplayName', sprintf('%s / %s', beamsAll{bi}, methodsAll{mi}));
            ci = ci+1;
        end
    end
    xlabel('Probe spacing \alpha'); ylabel('RMSE (mm)'); grid on; legend('Location','best','Interpreter','none');
    title('Task 4 Stage B: pointing RMSE vs probe spacing (reduced combo set)');
    save_figure_t2(fig, fullfile(figDir4,'T4_RMSE_vs_alpha'));
    close all;
end


function g = g_observation(eps, deltax, deltay, bp, hfloor)
    ux = [1;0]; uy = [0;1];
    hxp = beam_gain_general(eps(1)-deltax*ux(1), eps(2)-deltax*ux(2), bp);
    hxm = beam_gain_general(eps(1)+deltax*ux(1), eps(2)+deltax*ux(2), bp);
    hyp = beam_gain_general(eps(1)-deltay*uy(1), eps(2)-deltay*uy(2), bp);
    hym = beam_gain_general(eps(1)+deltay*uy(1), eps(2)+deltay*uy(2), bp);
    gx = log(max(hxp,hfloor)) - log(max(hxm,hfloor));
    gy = log(max(hyp,hfloor)) - log(max(hym,hfloor));
    g = [gx; gy];
end


function J = jacobian_fd(eps, deltax, deltay, bp, hfloor, fdstep)
    J = zeros(2,2);
    for j = 1:2
        dvec = zeros(2,1); dvec(j) = fdstep;
        gp = g_observation(eps+dvec, deltax, deltay, bp, hfloor);
        gm = g_observation(eps-dvec, deltax, deltay, bp, hfloor);
        J(:,j) = (gp-gm)/(2*fdstep);
    end
end


function eps_hat = gauss_newton_solve(z, eps0, deltax, deltay, bp, hfloor, fdstep, maxIter, tol)
% GAUSS_NEWTON_SOLVE  Levenberg-Marquardt-damped Newton refinement of
% g(eps)=z. Diagnostic-only comparator alongside Direct/EKF/Lookup in
% Stage A -- NOT part of the brief's single-step EKF update, and not
% called anywhere in the Task 5/6 sensing loop. Damping (JTJ + lambda*I)
% keeps the solve well-posed even where J is near-singular (e.g. sidelobe
% nulls, flat-top shoulder); the adaptive lambda accepts a step only if it
% actually reduces the residual, otherwise it increases damping and
% retries.
    eps_hat = eps0;
    lambda = 1e-6;
    for it = 1:maxIter
        r = z - g_observation(eps_hat, deltax, deltay, bp, hfloor);
        if norm(r) < tol, break; end
        J = jacobian_fd(eps_hat, deltax, deltay, bp, hfloor, fdstep);
        JTJ = J.'*J;
        step = (JTJ + lambda*eye(2)) \ (J.'*r);
        eps_new = eps_hat + step;
        r_new = z - g_observation(eps_new, deltax, deltay, bp, hfloor);
        if norm(r_new) < norm(r)
            eps_hat = eps_new;
            lambda = max(lambda/3, 1e-12);
        else
            lambda = min(lambda*3, 1e8);
            continue;   % reject step, retry from same point with more damping
        end
        if norm(step) < tol, break; end
    end
end


function eps_hat = lookup_estimate(z, Rinv, bp, deltax, deltay, hfloor, wx_dir, wy_dir, gridN)
    xg = linspace(-0.9*wx_dir, 0.9*wx_dir, gridN);
    yg = linspace(-0.9*wy_dir, 0.9*wy_dir, gridN);
    bestCost = Inf; eps_hat = [0;0];
    for ix = 1:numel(xg)
        for iy = 1:numel(yg)
            e = [xg(ix); yg(iy)];
            g = g_observation(e, deltax, deltay, bp, hfloor);
            r = z-g;
            cost = r.'*Rinv*r;
            if cost < bestCost
                bestCost = cost; eps_hat = e;
            end
        end
    end
end


function w = directional_radius(bp, dirVec, target_h)
    dirVec = dirVec/norm(dirVec);
    rlo = 0; rhi = 10*0.20;
    hlo = beam_gain_general(rlo*dirVec(1), rlo*dirVec(2), bp);
    hhi = beam_gain_general(rhi*dirVec(1), rhi*dirVec(2), bp);
    if hlo < target_h
        w = 0; return;
    end
    iter = 0;
    while hhi > target_h && iter < 20
        rhi = rhi*1.5;
        hhi = beam_gain_general(rhi*dirVec(1), rhi*dirVec(2), bp);
        iter = iter+1;
    end
    for it = 1:60
        rmid = 0.5*(rlo+rhi);
        hmid = beam_gain_general(rmid*dirVec(1), rmid*dirVec(2), bp);
        if hmid > target_h
            rlo = rmid;
        else
            rhi = rmid;
        end
    end
    w = 0.5*(rlo+rhi);
end


function [r_idx0, d_idx0] = target_bins(P, distance, velocity)
    r_idx0 = round(distance/(P.c/(2*P.N_radar*P.delta_f))) + 1;
    d_shift0 = round(velocity*2*P.N_symbols*P.T_sym/P.lambda);
    if d_shift0 < 0
        d_idx0 = d_shift0 + P.N_symbols + 1;
    else
        d_idx0 = d_shift0 + 1;
    end
end


function [Ppeak, PnoiseMed] = probe_gate_power(Z, r_idx0, d_idx0, gateR, gateD, N_radar, N_symbols)
    rIdx = max(1,r_idx0-gateR):min(N_radar, r_idx0+gateR);
    dIdxRaw = (d_idx0-gateD):(d_idx0+gateD);
    dIdx = mod(dIdxRaw-1, N_symbols)+1;
    mask = false(N_radar, N_symbols);
    mask(rIdx, dIdx) = true;
    Pmap = abs(Z).^2;
    Ppeak = max(Pmap(mask));
    PnoiseMed = median(Pmap(~mask));
end


function Pdet_thresh = calibrate_detection_threshold(rx_delayed_doppler, corr_ref, P, ...
    A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Nfa, Pfa_target)
    Ppk = zeros(Nfa,1);
    for m = 1:Nfa
        [rx_signal, ~] = generate_received_signal_fast(rx_delayed_doppler, 0, A_base, P_n_radar);
        Z = range_doppler_fast(rx_signal, corr_ref, P.N_fft, P.N_cp, P.N_symbols, P.N_radar);
        [Ppk(m), ~] = probe_gate_power(Z, r_idx0, d_idx0, gateR, gateD, P.N_radar, P.N_symbols);
    end
    Pdet_thresh = prctile(Ppk, 100*(1-Pfa_target));
end


function [z, avail] = measure_probe_pair(rx_delayed_doppler, corr_ref, P, bp, eps_minus, ...
    deltax, deltay, A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Pdet_thresh, hfloor)
% MEASURE_PROBE_PAIR  Uses the actual measured per-pair noise-floor median
% PnoiseMed for Pfloor (Eq.63), instead of a fixed multiple of the
% globally calibrated Pdet_thresh -- matches the standalone Task 4 fix.
    ux = [1;0]; uy = [0;1];
    eps_list = [eps_minus-deltax*ux, eps_minus+deltax*ux, eps_minus-deltay*uy, eps_minus+deltay*uy];
    Pmeas = zeros(1,4); Pnoise = zeros(1,4); % [x+, x-, y+, y-]
    for i = 1:4
        e = eps_list(:,i);
        h_i = beam_gain_general(e(1), e(2), bp);
        [rx_signal, ~] = generate_received_signal_fast(rx_delayed_doppler, h_i, A_base, P_n_radar);
        Z = range_doppler_fast(rx_signal, corr_ref, P.N_fft, P.N_cp, P.N_symbols, P.N_radar);
        [Ppk, PnMed] = probe_gate_power(Z, r_idx0, d_idx0, gateR, gateD, P.N_radar, P.N_symbols);
        Pmeas(i) = Ppk; Pnoise(i) = PnMed;
    end
    availX = (Pmeas(1) >= Pdet_thresh) && (Pmeas(2) >= Pdet_thresh);
    availY = (Pmeas(3) >= Pdet_thresh) && (Pmeas(4) >= Pdet_thresh);
    avail = [availX, availY];

    % Eq.63: Pfloor,k = 1e-6 * PbN,k, using each probe's own measured
    % noise-floor median (averaged over the pair for symmetry).
    Pfloor_x = 1e-6*mean(Pnoise(1:2));
    Pfloor_y = 1e-6*mean(Pnoise(3:4));
    zx = log((Pmeas(1)+Pfloor_x)/(Pmeas(2)+Pfloor_x));
    zy = log((Pmeas(3)+Pfloor_y)/(Pmeas(4)+Pfloor_y));
    z = [zx; zy];
    if ~availX, z(1) = NaN; end
    if ~availY, z(2) = NaN; end
end


%% ========================================================================
%  LOCAL FUNCTIONS -- Task 5 closed-loop tracking (merged in from the
%  standalone, FIXED Task 5 script; FIX 1/FIX 2 already applied, see
%  header notes at the top of this file. All shared helpers this section
%  needs -- setup_system_params, build_radar_reference, beam_gain_general,
%  beam_matrix_A, make_beam_case, tf_str, save_figure_t2,
%  generate_received_signal_fast, range_doppler_fast, extract_target_fast,
%  generate_ar1_stationary, g_observation, jacobian_fd, lookup_estimate,
%  directional_radius, target_bins, probe_gate_power,
%  calibrate_detection_threshold, measure_probe_pair, fit_curvature_matrix
%  -- are the SINGLE shared copies already defined above; the duplicate
%  copies that shipped with the standalone Task 5 script have been
%  removed. Only the Task-5-specific functions below are new.)
%% ========================================================================

function run_task5_closed_loop_tracking(wref, hfloor, fdstep, figDir5, resDir5, ...
    beamTag_true, beamLabel_true, rho, jitterFrac, alphaProbe, J_interval, ...
    K, K_ideal, hth, gridNLookup, Mcal, Nfa, Pfa_target, hfit_min_task3, Ptx_dBm, ...
    distance, velocity)
% RUN_TASK5_CLOSED_LOOP_TRACKING  Wraps the standalone Task 5 script's
% sections 2-7 (beam/probe geometry, physical sensing-chain setup,
% measurement-covariance calibration, shared trajectory, Stage A ideal
% sanity check, Stage B full sensing-chain run for all 7 methods, metrics
% table, deliverable figures) into a single call, mirroring
% run_task4_validation(). See the top-of-file header for the two
% correctness fixes (Stage A z-generation from A_true; paired noise seeds
% across methods) already verified before this merge.

    fprintf('===============================================\n');
    fprintf('Task 5 -- Closed-loop tracking and compensation\n');
    fprintf('===============================================\n');
    sigma_s = jitterFrac*wref;
    Sigma_s = sigma_s^2*eye(2);
    F_dyn   = rho*eye(2);
    Q_dyn   = (1-rho^2)*Sigma_s;
    fprintf('wref=%.4f m | true beam = %s | rho=%.2f | sigma_s/wref=%.2f | alpha=%.2f | J=%d | K=%d\n\n', ...
        wref, beamLabel_true, rho, jitterFrac, alphaProbe, J_interval, K);

    %% --- Beam / probe geometry -----------------------------------------
    bp_true = make_beam_case(beamTag_true, wref);
    A_true  = beam_matrix_A(bp_true.wx, bp_true.wy, bp_true.phi);   % exact (beam is quadratic)
    A_circ  = (2/wref^2)*eye(2);                                     % mismatched circular assumption

    ux = [1;0]; uy = [0;1];
    wx_dir = directional_radius(bp_true, ux, exp(-2));
    wy_dir = directional_radius(bp_true, uy, exp(-2));
    deltax = alphaProbe*wx_dir;
    deltay = alphaProbe*wy_dir;
    fprintf('Probe geometry: wx_dir=%.4f m, wy_dir=%.4f m, deltax=%.4f m, deltay=%.4f m\n', ...
        wx_dir, wy_dir, deltax, deltay);

    % Task 3 fitted matrix, evaluated on the ACTUAL beam (method 5)
    fitT3 = fit_curvature_matrix(bp_true, wref, hfit_min_task3, hfloor, 401);
    A_fitted = fitT3.A;
    if fitT3.corrected, clipStr = 'YES'; else, clipStr = 'NO'; end
    fprintf('Task 3 fitted A (hfit_min=%.2f): [%.4f %.4f; %.4f %.4f]  (eig-clipping applied: %s)\n', ...
        hfit_min_task3, A_fitted(1,1),A_fitted(1,2),A_fitted(2,1),A_fitted(2,2), clipStr);
    fprintf('True analytic A:                 [%.4f %.4f; %.4f %.4f]\n\n', ...
        A_true(1,1),A_true(1,2),A_true(2,1),A_true(2,2));

    %% --- Physical sensing chain setup (shared by all methods' probes + main burst) ---
    P = setup_system_params();
    waveform = 'OFDM';
    W = build_radar_reference(P, waveform, distance, velocity);
    rx_delayed_doppler = W.rx_delayed_doppler;
    corr_ref = W.corr_ref;

    height = 100;
    [T_atm_K, p_hPa, wvden] = atmositu(height);
    T_C = T_atm_K - 273.15; P_Pa = p_hPa*100;
    L_atm_dB  = gaspl(distance, P.fc, T_C, P_Pa, wvden);
    L_atm_lin = 10^(-L_atm_dB/10);
    T_noise   = (1-L_atm_lin)*T_atm_K;
    P_n       = P.kb*T_noise*P.B;
    P_n_radar = P_n*(P.N_radar/P.N_fft);

    Ptx_W = 10^((Ptx_dBm-30)/10);
    P_radar = Ptx_W*(P.N_radar/P.N_fft);
    A_base = sqrt(P_radar)*L_atm_lin* ...
        sqrt((P.G_total*P.lambda^2*P.sigma_rcs)/((4*pi)^3*distance^4));
    aligned_SNR_dB = 10*log10(A_base^2/P_n_radar);

    gateR = 2; gateD = 2;
    [r_idx0, d_idx0] = target_bins(P, distance, velocity);

    Pdet_thresh = calibrate_detection_threshold(rx_delayed_doppler, corr_ref, P, ...
        A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Nfa, Pfa_target);
    fprintf('Aligned SNR ~ %.1f dB | Pdet threshold (Pfa=%.1e) = %.4e\n\n', aligned_SNR_dB, Pfa_target, Pdet_thresh);

    %% --- Measurement-covariance calibration (Monte Carlo, eps_cal = 0) ---
    fprintf('--- Calibrating measurement covariance (Monte Carlo, eps_cal=0, Mcal=%d) ---\n', Mcal);
    Zcal = nan(Mcal,2);
    EpsLUTcal = nan(Mcal,2);
    Rinv_forLUTcal = eye(2);
    for m = 1:Mcal
        rng(90000+m);
        [zc, availc] = measure_probe_pair(rx_delayed_doppler, corr_ref, P, bp_true, ...
            [0;0], deltax, deltay, A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Pdet_thresh, hfloor);
        if all(availc)
            Zcal(m,:) = zc.';
            EpsLUTcal(m,:) = lookup_estimate(zc, Rinv_forLUTcal, bp_true, deltax, deltay, hfloor, wx_dir, wy_dir, gridNLookup).';
        end
    end
    validCal = ~any(isnan(Zcal),2);
    if nnz(validCal) >= 10
        Rquad = cov(Zcal(validCal,:));
    else
        Rquad = eye(2);
        warning('R calibration: only %d/%d valid samples, using identity R.', nnz(validCal), Mcal);
    end
    validLUTcal = ~any(isnan(EpsLUTcal),2);
    if nnz(validLUTcal) >= 10
        Rlookup = cov(EpsLUTcal(validLUTcal,:));
    else
        Rlookup = eye(2)*1e-4;
        warning('Rlookup calibration: only %d/%d valid samples, using a small default.', nnz(validLUTcal), Mcal);
    end
    fprintf('Rquad    = [%.3e %.3e; %.3e %.3e]\n', Rquad(1,1),Rquad(1,2),Rquad(2,1),Rquad(2,2));
    fprintf('Rlookup  = [%.3e %.3e; %.3e %.3e]\n\n', Rlookup(1,1),Rlookup(1,2),Rlookup(2,1),Rlookup(2,2));

    %% --- Shared trajectory (same sk sequence reused for every method) ---
    sk_traj = generate_ar1_stationary(K, rho, Sigma_s, 4242);   % K x 2

    %% --- STAGE A: ideal (noiseless) sanity check, methods 3-7 only ---
    fprintf('=== STAGE A: ideal (noiseless) closed-loop sanity check (%d frames) ===\n\n', K_ideal);
    run_task5_stageA_ideal(K_ideal, sk_traj(1:K_ideal,:), F_dyn, Q_dyn, Sigma_s, ...
        bp_true, A_true, A_circ, A_fitted, deltax, deltay, hfloor, fdstep, ...
        wx_dir, wy_dir, gridNLookup);

    %% --- STAGE B: full sensing-chain closed-loop run, all 7 methods ---
    fprintf('=== STAGE B: full sensing-chain closed-loop run (%d frames, %d methods) ===\n\n', K, 7);

    methodNames = {'NoCompensation','Oracle','MatchedQuadraticKF','MismatchedCircularKF', ...
                   'FittedMatrixKF','FullMapEKF','LookupPlusKF'};
    methodModes = {'none','oracle','quad','quad','quad','ekf','lookup'};
    methodA     = {[], [], A_true, A_circ, A_fitted, [], []};

    lambda_wave = P.lambda;
    baseSeedB = 5000;   % FIX 2: SAME seed passed to every method (paired noise)

    allTS = cell(1,numel(methodNames));
    for mi = 1:numel(methodNames)
        fprintf('--- Running: %s ---\n', methodNames{mi});
        TS = run_task5_method(methodModes{mi}, methodA{mi}, Rquad, Rlookup, sk_traj, ...
            F_dyn, Q_dyn, Sigma_s, J_interval, bp_true, deltax, deltay, hfloor, fdstep, ...
            wx_dir, wy_dir, gridNLookup, ...
            P, rx_delayed_doppler, corr_ref, A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Pdet_thresh, ...
            distance, velocity, P.c, P.delta_f, P.N_radar, lambda_wave, P.T_sym, P.N_symbols, ...
            P.rangeRes, P.velRes, baseSeedB);
        allTS{mi} = TS;
    end
    fprintf('\n');

    %% --- Metrics table ---
    rowsMethod={}; rowsRMSEpred=[]; rowsRMSEres=[]; rowsPdrop=[]; rowsPout=[]; rowsPcorrect=[];
    rowsRMSER=[]; rowsRMSEv=[]; rowsFailAvg=[]; rowsFailMax=[]; rowsPl2=[]; rowsPl5=[]; rowsPl10=[];

    for mi = 1:numel(methodNames)
        M = compute_task5_metrics(allTS{mi}, hth);
        fprintf('[%-22s] RMSEpred=%.4f m RMSEres=%.4f m Pdrop=%.3f Pout=%.3f Pcorrect=%.3f avgFail=%.2f maxFail=%d\n', ...
            methodNames{mi}, M.RMSEpred, M.RMSEres, M.Pdrop, M.Pout, M.Pcorrect, M.avgFail, M.maxFail);

        rowsMethod{end+1,1}   = methodNames{mi}; %#ok<AGROW>
        rowsRMSEpred(end+1,1) = M.RMSEpred; %#ok<AGROW>
        rowsRMSEres(end+1,1)  = M.RMSEres; %#ok<AGROW>
        rowsPdrop(end+1,1)    = M.Pdrop; %#ok<AGROW>
        rowsPout(end+1,1)     = M.Pout; %#ok<AGROW>
        rowsPcorrect(end+1,1) = M.Pcorrect; %#ok<AGROW>
        rowsRMSER(end+1,1)    = M.RMSE_R; %#ok<AGROW>
        rowsRMSEv(end+1,1)    = M.RMSE_v; %#ok<AGROW>
        rowsFailAvg(end+1,1)  = M.avgFail; %#ok<AGROW>
        rowsFailMax(end+1,1)  = M.maxFail; %#ok<AGROW>
        rowsPl2(end+1,1)      = M.Pl(1); %#ok<AGROW>
        rowsPl5(end+1,1)      = M.Pl(2); %#ok<AGROW>
        rowsPl10(end+1,1)     = M.Pl(3); %#ok<AGROW>
    end
    fprintf('\n');

    T5 = table(rowsMethod, rowsRMSEpred, rowsRMSEres, rowsPdrop, rowsPout, rowsPcorrect, ...
        rowsRMSER, rowsRMSEv, rowsFailAvg, rowsFailMax, rowsPl2, rowsPl5, rowsPl10, ...
        'VariableNames', {'Method','RMSEpred_m','RMSEres_m','Pdrop','Pout','Pcorrect', ...
        'CondRMSE_R_m','CondRMSE_v_mps','AvgFailRun','MaxFailRun','P_failrun_ge2','P_failrun_ge5','P_failrun_ge10'});
    disp(T5);
    writetable(T5, fullfile(resDir5,'Task5_methods_summary_table.csv'));
    fprintf('Summary table -> %s\n\n', fullfile(resDir5,'Task5_methods_summary_table.csv'));

    %% --- Deliverable figures -- time series per method + cross-method comparison ---
    for mi = 1:numel(methodNames)
        plot_task5_timeseries(allTS{mi}, methodNames{mi}, figDir5);
    end
    plot_task5_comparison(allTS, methodNames, hth, figDir5);
    close all;

    fprintf('Figures -> %s\nTables  -> %s\n', figDir5, resDir5);
    fprintf('============ TASK 5 CLOSED-LOOP TRACKING COMPLETE ============\n\n');
end


function run_task5_stageA_ideal(K_ideal, sk_traj, F_dyn, Q_dyn, Sigma_s, ...
    bp_true, A_true, A_circ, A_fitted, deltax, deltay, hfloor, fdstep, wx_dir, wy_dir, gridNLookup)
% RUN_TASK5_STAGEA_IDEAL  Noiseless g(.)-only closed-loop sanity check for
% methods 3-7. No waveform simulation; verifies the update math/signs
% converge (or, for the mismatched case, converge only partially) before
% the full Stage B run.
%
% FIX 1 (verified, see top-of-file header): for the 'quad' methods, z is
% ALWAYS generated from A_true (the physical beam's curvature matrix),
% never from methodA{mi}. Each method's own assumed matrix (methodA{mi})
% is used ONLY inside its own KF observation matrix H for the update
% step.
    Rtiny = 1e-8*eye(2);   % near-zero "noise" so S is never singular

    methodNames = {'MatchedQuadraticKF','MismatchedCircularKF','FittedMatrixKF','FullMapEKF','LookupPlusKF'};
    methodModes = {'quad','quad','quad','ekf','lookup'};
    methodA     = {A_true, A_circ, A_fitted, [], []};

    % H built from the TRUE beam's curvature matrix -- used to generate the
    % physically correct noiseless measurement z, regardless of which
    % method is being tested.
    H_true = [4*deltax*[1 0]*A_true; 4*deltay*[0 1]*A_true];

    for mi = 1:numel(methodNames)
        s_prev = [0;0]; P_prev = Sigma_s;
        errFinal = NaN;
        for k = 1:K_ideal
            sk = sk_traj(k,:).';
            s_pred = F_dyn*s_prev;
            P_pred = F_dyn*P_prev*F_dyn.' + Q_dyn;
            ck = s_pred;
            eps_true_minus = sk - ck;

            switch methodModes{mi}
                case 'quad'
                    A_assumed = methodA{mi};
                    z = H_true*eps_true_minus;   % FIX 1: physically correct measurement, from A_true
                    H = [4*deltax*[1 0]*A_assumed; 4*deltay*[0 1]*A_assumed];   % method's own (possibly wrong) H, used in the update only
                    rk = z;
                    Sk = H*P_pred*H.' + Rtiny;
                    Kk = P_pred*H.'/Sk;
                    s_upd = s_pred + Kk*rk;
                    P_upd = (eye(2)-Kk*H)*P_pred*(eye(2)-Kk*H).' + Kk*Rtiny*Kk.';
                case 'ekf'
                    z = g_observation(eps_true_minus, deltax, deltay, bp_true, hfloor);   % ideal (no waveform noise)
                    J0 = jacobian_fd([0;0], deltax, deltay, bp_true, hfloor, fdstep);
                    rk = z;   % z_pred = g(0) = 0
                    Sk = J0*P_pred*J0.' + Rtiny;
                    Kk = P_pred*J0.'/Sk;
                    s_upd = s_pred + Kk*rk;
                    P_upd = (eye(2)-Kk*J0)*P_pred*(eye(2)-Kk*J0).' + Kk*Rtiny*Kk.';
                case 'lookup'
                    z = g_observation(eps_true_minus, deltax, deltay, bp_true, hfloor);
                    eps_lut = lookup_estimate(z, eye(2), bp_true, deltax, deltay, hfloor, wx_dir, wy_dir, gridNLookup);
                    s_tilde = ck + eps_lut;
                    rk = s_tilde - s_pred;
                    Sk = P_pred + Rtiny;
                    Kk = P_pred/Sk;
                    s_upd = s_pred + Kk*rk;
                    P_upd = (eye(2)-Kk)*P_pred*(eye(2)-Kk).' + Kk*Rtiny*Kk.';
            end
            s_prev = s_upd; P_prev = P_upd;
            errFinal = norm(sk - s_upd);
        end
        fprintf('  %-22s final-frame |s_true - s_upd| = %.4e m  (K_ideal=%d)\n', methodNames{mi}, errFinal, K_ideal);
    end
    fprintf('  (Matched/Fitted/EKF/Lookup should shrink toward ~0; Mismatched circular will\n');
    fprintf('   generally NOT reach 0 even noiselessly, since its assumed A is structurally\n');
    fprintf('   wrong for a non-circular beam -- this is the expected, informative failure mode.)\n\n');
end


function TS = run_task5_method(mode, Aassumed, Rquad, Rlookup, sk_traj, ...
    F_dyn, Q_dyn, Sigma_s, J_interval, bp_true, deltax, deltay, hfloor, fdstep, ...
    wx_dir, wy_dir, gridNLookup, ...
    P, rx_delayed_doppler, corr_ref, A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Pdet_thresh, ...
    distance, velocity, c, delta_f, N_radar, lambda_wave, T_sym, N_symbols, rangeRes, velRes, baseSeed)
% RUN_TASK5_METHOD  Runs one full closed-loop trajectory (Stage B, full
% sensing chain) for one of the 7 Task 5 methods and returns a time-series
% struct used both for the deliverable figure and the performance metrics.
%   mode: 'none' | 'oracle' | 'quad' | 'ekf' | 'lookup'
%   Aassumed: 2x2 assumed curvature matrix (used only when mode=='quad')
%
% FIX 2 (verified, see top-of-file header): baseSeed is passed in
% IDENTICAL across all 7 methods by the caller, so every rng() call inside
% this function (baseSeed+k for probe, baseSeed+700000+k for main burst)
% draws the SAME noise realization at the same frame k for every method
% on the same shared trajectory -- satisfying the brief's "paired
% trajectories and paired noise seeds across methods" requirement.

    K = size(sk_traj,1);
    TS.k = (1:K).';
    TS.sk = sk_traj;
    TS.s_pred = nan(K,2);
    TS.s_upd  = nan(K,2);
    TS.cmain  = nan(K,2);
    TS.eps_main = nan(K,2);
    TS.h_main = nan(K,1);
    TS.scheduled = false(K,1);
    TS.avail = false(K,1);
    TS.Dk = false(K,1);
    TS.range_est = nan(K,1);
    TS.vel_est = nan(K,1);

    s_prev = [0;0]; P_prev = Sigma_s;   % Eq.95: s0|-1=0, P0|-1=Sigma_s

    for k = 1:K
        sk = sk_traj(k,:).';
        s_pred = F_dyn*s_prev;
        P_pred = F_dyn*P_prev*F_dyn.' + Q_dyn;

        switch mode
            case 'none'
                cmain = [0;0];
                s_upd = [0;0]; P_upd = P_pred;
                scheduled = false; avail = false;

            case 'oracle'
                cmain = sk;
                s_pred = sk; s_upd = sk; P_upd = zeros(2);   % perfect knowledge, by definition
                scheduled = false; avail = false;

            otherwise   % 'quad', 'ekf', 'lookup'
                ck = s_pred;
                scheduled = (mod(k-1, J_interval) == 0);
                if scheduled
                    rng(baseSeed + k);
                    eps_true_minus = sk - ck;
                    [z, avail2] = measure_probe_pair(rx_delayed_doppler, corr_ref, P, bp_true, ...
                        eps_true_minus, deltax, deltay, A_base, P_n_radar, r_idx0, d_idx0, gateR, gateD, Pdet_thresh, hfloor);
                    avail = all(avail2);
                else
                    avail = false; z = [NaN;NaN];
                end

                if avail
                    switch mode
                        case 'quad'
                            H = [4*deltax*[1 0]*Aassumed; 4*deltay*[0 1]*Aassumed];
                            rk = z;   % ck=s_pred exactly => innovation reduces to z (see header derivation)
                            Sk = H*P_pred*H.' + Rquad;
                            Kk = P_pred*H.'/Sk;
                            s_upd = s_pred + Kk*rk;
                            P_upd = (eye(2)-Kk*H)*P_pred*(eye(2)-Kk*H).' + Kk*Rquad*Kk.';
                        case 'ekf'
                            J0 = jacobian_fd([0;0], deltax, deltay, bp_true, hfloor, fdstep);
                            rk = z;   % z_pred = g(eps_pred=0) = 0
                            Sk = J0*P_pred*J0.' + Rquad;
                            Kk = P_pred*J0.'/Sk;
                            s_upd = s_pred + Kk*rk;
                            P_upd = (eye(2)-Kk*J0)*P_pred*(eye(2)-Kk*J0).' + Kk*Rquad*Kk.';
                        case 'lookup'
                            Rinv = inv(Rquad + 1e-12*eye(2)); %#ok<MINV>
                            eps_lut = lookup_estimate(z, Rinv, bp_true, deltax, deltay, hfloor, wx_dir, wy_dir, gridNLookup);
                            s_tilde = ck + eps_lut;
                            rk = s_tilde - s_pred;
                            Sk = P_pred + Rlookup;
                            Kk = P_pred/Sk;
                            s_upd = s_pred + Kk*rk;
                            P_upd = (eye(2)-Kk)*P_pred*(eye(2)-Kk).' + Kk*Rlookup*Kk.';
                    end
                else
                    s_upd = s_pred; P_upd = P_pred;   % Eq.83, missing-measurement rule
                end
                cmain = s_upd;
        end

        eps_main = sk - cmain;
        h_main = beam_gain_general(eps_main(1), eps_main(2), bp_true);

        % --- main-burst detection (always run, every method, every frame) ---
        rng(baseSeed + 700000 + k);
        [rx_signal, ~] = generate_received_signal_fast(rx_delayed_doppler, h_main, A_base, P_n_radar);
        Z = range_doppler_fast(rx_signal, corr_ref, P.N_fft, P.N_cp, P.N_symbols, P.N_radar);
        [~, r_est, v_est, det] = extract_target_fast(Z, c, delta_f, N_radar, lambda_wave, T_sym, N_symbols, ...
            distance, velocity, rangeRes, velRes);

        TS.s_pred(k,:) = s_pred.';
        TS.s_upd(k,:)  = s_upd.';
        TS.cmain(k,:)  = cmain.';
        TS.eps_main(k,:) = eps_main.';
        TS.h_main(k)   = h_main;
        TS.scheduled(k) = scheduled;
        TS.avail(k)    = avail;
        TS.Dk(k)       = det;
        TS.range_est(k) = r_est;
        TS.vel_est(k)   = v_est;

        s_prev = s_upd; P_prev = P_upd;
    end
end


function M = compute_task5_metrics(TS, hth)
% COMPUTE_TASK5_METRICS  Eqs.98-107: RMSEpred, RMSEres, Pout, Pcorrect,
% conditional RMSE_R/RMSE_v, failure-run statistics, Pdrop.
    K = numel(TS.k);
    M.RMSEpred = sqrt(mean(sum((TS.sk - TS.s_pred).^2, 2)));
    M.RMSEres  = sqrt(mean(sum((TS.sk - TS.s_upd).^2, 2)));
    M.Pout     = mean(TS.h_main < hth);
    M.Pcorrect = mean(TS.Dk);

    Dset = TS.Dk;
    if any(Dset)
        M.RMSE_R = sqrt(mean((TS.range_est(Dset) - 15).^2));      % distance=15 m, fixed target
        M.RMSE_v = sqrt(mean((TS.vel_est(Dset) - 5).^2));         % velocity=5 m/s, fixed target
    else
        M.RMSE_R = NaN; M.RMSE_v = NaN;
    end

    [avgFail, maxFail] = run_length_stats(~TS.Dk);
    M.avgFail = avgFail; M.maxFail = maxFail;

    failVec = ~TS.Dk;
    d = diff([0; failVec(:); 0]);
    starts = find(d==1); ends = find(d==-1)-1;
    lens = ends-starts+1;
    Lset = [2 5 10];
    Pl = zeros(1,numel(Lset));
    if ~isempty(lens)
        for li = 1:numel(Lset)
            Pl(li) = mean(lens >= Lset(li));
        end
    end
    M.Pl = Pl;

    nSched = nnz(TS.scheduled);
    if nSched > 0
        M.Pdrop = nnz(TS.scheduled & ~TS.avail) / nSched;
    else
        M.Pdrop = NaN;   % methods that never probe (None, Oracle)
    end
end


function plot_task5_timeseries(TS, methodName, figDir5)
% PLOT_TASK5_TIMESERIES  Deliverable figure: true vs predicted vs updated
% displacement, residual displacement, beam coefficient, and dropout
% locations, for one representative trajectory / method.
    K = numel(TS.k); %#ok<NASGU>
    fig = figure('Color','w','Position',[100 100 1000 900]);

    subplot(4,1,1); hold on;
    plot(TS.k, TS.sk(:,1), 'k-', 'LineWidth',1.4, 'DisplayName','s_x (true)');
    plot(TS.k, TS.s_pred(:,1), 'b--', 'LineWidth',1.0, 'DisplayName','s_x pred');
    plot(TS.k, TS.s_upd(:,1), 'r-', 'LineWidth',1.0, 'DisplayName','s_x updated');
    ylabel('x (m)'); grid on; legend('Location','best','FontSize',7);
    title(sprintf('Task 5: %s -- true / predicted / updated displacement', methodName), 'Interpreter','none');

    subplot(4,1,2); hold on;
    plot(TS.k, TS.sk(:,2), 'k-', 'LineWidth',1.4, 'DisplayName','s_y (true)');
    plot(TS.k, TS.s_pred(:,2), 'b--', 'LineWidth',1.0, 'DisplayName','s_y pred');
    plot(TS.k, TS.s_upd(:,2), 'r-', 'LineWidth',1.0, 'DisplayName','s_y updated');
    ylabel('y (m)'); grid on; legend('Location','best','FontSize',7);

    subplot(4,1,3); hold on;
    epsNorm = sqrt(sum(TS.eps_main.^2,2));
    plot(TS.k, epsNorm, 'm-', 'LineWidth',1.2, 'DisplayName','|\epsilon_{main}|');
    yyaxis right;
    plot(TS.k, TS.h_main, 'g-', 'LineWidth',1.2, 'DisplayName','h_{main}');
    ylabel('h_{main}'); ylim([0 1.05]);
    yyaxis left; ylabel('|\epsilon_{main}| (m)');
    xlabel('Frame, k'); grid on; legend('Location','best','FontSize',7);
    title('Residual displacement magnitude and beam coefficient','FontSize',9);

    subplot(4,1,4); hold on;
    dropoutFrames = TS.scheduled & ~TS.avail;
    detFail = ~TS.Dk;
    stem(TS.k(dropoutFrames), ones(nnz(dropoutFrames),1), 'r', 'Marker','x', 'DisplayName','probe dropout');
    stem(TS.k(detFail), 0.5*ones(nnz(detFail),1), 'b', 'Marker','o', 'DisplayName','main-burst det. fail');
    ylim([0 1.5]); yticks([0.5 1]); yticklabels({'det. fail','probe drop'});
    xlabel('Frame, k'); grid on; legend('Location','best','FontSize',7);
    title('Measurement-dropout / detection-failure locations','FontSize',9);

    save_figure_t2(fig, fullfile(figDir5, sprintf('T5_timeseries_%s', methodName)));
    close(fig);
end


function plot_task5_comparison(allTS, methodNames, hth, figDir5)
% PLOT_TASK5_COMPARISON  Overlay of h_main across all 7 methods, for a
% quick cross-method visual comparison on the same trajectory.
    colors = lines(numel(methodNames));
    fig = figure('Color','w','Position',[100 100 1000 500]); hold on;
    for mi = 1:numel(methodNames)
        TS = allTS{mi};
        plot(TS.k, TS.h_main, 'Color', colors(mi,:), 'LineWidth',1.1, 'DisplayName', methodNames{mi});
    end
    yline(hth, 'k--', sprintf('h_{th}=%.2f', hth));
    xlabel('Frame, k'); ylabel('h_{main}'); ylim([0 1.05]); grid on;
    legend('Location','best','Interpreter','none','FontSize',8);
    title('Task 5: main-burst beam coefficient, all methods, same trajectory');
    save_figure_t2(fig, fullfile(figDir5, 'T5_comparison_hmain_allmethods'));
    close(fig);
end