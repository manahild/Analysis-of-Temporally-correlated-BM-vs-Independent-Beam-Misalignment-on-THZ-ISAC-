%% ========================================================================
%  TASKS 2, 3 -- GENERALIZED BEAM-MISALIGNMENT LIBRARY, CURVATURE
%  FIT/VALIDATION, AND THE FULL CORRELATED-MISALIGNMENT SENSING SIMULATOR
%  ------------------------------------------------------------------------
%  This is the "Tasks 2" integrated script with the standalone Task 3
%  curvature-fitting script merged in, exactly as flagged in Task 3's own
%  header: "This file is written to be dropped into the Task 2/5/6 master
%  script later: it reuses the Task 2 beam-library local functions
%  VERBATIM ... so there is nothing to reconcile except deleting the
%  duplicate copies when merging."
%
%  WHAT CHANGED VS. THE TWO SOURCE FILES (Tasks2/5/6 master + standalone
%  Task 3 script)
%  ----------------------------------------------------------------------
%  1) Task 3's duplicate copies of beam_gain_general, mysinc2,
%     beam_matrix_A, make_beam_case, tf_str, fmt_or_na, save_figure_t2
%     were DELETED -- this file now has exactly one copy of each, shared
%     by Task 2, Task 3, and the Task 5/6 sensing loop.
%  2) Task 3's remaining local functions (center_hessian_curvature,
%     fit_curvature_matrix, compute_fitting_errors, compute_r1dB,
%     run_task3_validation) are added verbatim (formulas/logic
%     unchanged) as new local functions of this file.
%  3) Task 3's placeholder "wref = 0.20 m (ASSUMPTION)" is replaced by the
%     SAME wref2 = w_mo of the SF regime that Task 2's validation section
%     already uses in this master script -- i.e. Task 3 now runs on the
%     ACTUAL beam-waist of the regime under test, not a hardcoded stand-in,
%     exactly mirroring how Task 2's own placeholder was resolved when it
%     was first wired into this master file.
%  4) run_task3_validation() is called once, immediately after
%     run_task2_validation(), using the same wref2, grid size, and the
%     5-beam-case list already defined for Task 2/5/6. Its outputs
%     (Figures_T3/, Results_T3/Task3_curvature_fit_table.csv) are
%     produced alongside Task 2's and Task 5/6's outputs, with no changes
%     to the Task 5/6 sensing loop, radar equation, range-Doppler
%     processing, or Section 10/11 stationary-covariance check.
%  5) hfitMinList, hfloor, and a dedicated Ngrid3 are declared once in
%     Section 2.5 (new) and passed into run_task3_validation(), matching
%     the sweep {0.8, 0.5, 0.25, 0.1} required by the brief (Eq.39).
%
%  Everything else about the radar waveform, range-Doppler processor,
%  AR(1) recursion, metric definitions, and Task 2/5/6 loop structure is
%  UNCHANGED from the Task 2/5/6 master script.
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
figDir  = fullfile(scriptDir,'Figures_T56');    % Task 5/6 sensing figures
resDir  = fullfile(scriptDir,'Results_T56');    % Task 5/6 sensing table
datDir  = fullfile(scriptDir,'Data_T56');       % Task 5/6 saved data
for f = {figDir2,resDir2,figDir3,resDir3,figDir,resDir,datDir}
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
%     validation section, the Task 3 curvature-fit section, AND the
%     Task 5/6 sensing loop below, via make_beam_case(name, wref).
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

%% ------------------------------------------------------------------
%  3. SENSING / MISALIGNMENT SCENARIO SETUP  (Task 5 requirements)
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
%  10. SAVE DATA (now includes the beam-case dimension)
%% ------------------------------------------------------------------
save(fullfile(datDir,'Task5_6_results.mat'), 'results', 'T6', ...
    'distance','theta_bm','sigma_s','Ptx_dBm_fixed','N_frames','h_th', ...
    'rangeRes','velRes','regimes','waveforms','scenarioNames','scenarioTagsSF', ...
    'beamCaseNames','beamCaseLabels', ...
    'w_0_SF','apert_T_SF','w_0_SL','apert_T_SL','apert_R','C_n_sq','fc', ...
    'hfitMinList','hfloor','Ngrid3');

fprintf('\nFigures -> %s\nTable   -> %s\nData    -> %s\n', figDir, resDir, datDir);
fprintf('============ TASKS 2, 3, 5 & 6 COMPLETE ============\n');

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
% ESTIMATE_POINTING_ERROR  RESERVED STUB (Task 4 scope). Unchanged.
    error('estimate_pointing_error: not implemented yet (Task 4 scope).');
end


function [xhat, Pcov] = filter_update(z, xpred, Ppred, model) %#ok<INUSD,DEFNU>
% FILTER_UPDATE  RESERVED STUB for the Kalman/EKF state update (Task 4/5 scope). Unchanged.
    if isempty(z)
        xhat = xpred;
        Pcov = Ppred;
        return;
    end
    error('filter_update: measurement update not implemented yet (Task 4/5 scope).');
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
%  Task 2 validation, Task 3 curvature fitting, AND the Task 5/6 sensing
%  loop above)
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
    targetFolder = fileparts(baseName);
    if ~exist(targetFolder, 'dir')
        error('save_figure_t2: target folder does not exist: %s', targetFolder);
    end
    try
        savefig(figHandle, [baseName '.fig']);
    catch
        % savefig not always available/needed under Octave batch mode
    end
    if exist('exportgraphics', 'file')
        exportgraphics(figHandle, [baseName '.png'], 'Resolution', 200);
    else
        print(figHandle, [baseName '.png'], '-dpng', '-r200');
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