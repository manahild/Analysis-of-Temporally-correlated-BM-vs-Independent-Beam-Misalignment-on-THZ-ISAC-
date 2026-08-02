%% ========================================================================
%  TASK 5 & 6 -- CORRELATED BEAM MISALIGNMENT IN THE FULL SENSING SIM
%  PURE 2026 Program
%
%  Plugs the Task 2-4 temporal misalignment model (Eq.5-6, Eq.7) into the
%  main THz-ISAC OFDM vs SC-FDMA range-Doppler simulation, per Eq.(8):
%       Pr,k = Pr,aligned * h_k        (LINEAR in h, NOT h^2)
%
%  Runs a sequence of sensing frames for 4 scenarios:
%       1. Aligned          (h_k = 1 for all k)
%       2. i.i.d. misalignment
%       3. Correlated, rho = 0.9
%       4. Correlated, rho = 0.99
%  for both regimes (SF, SL) and both waveforms (OFDM, SC-FDMA).
%
%  *** IMPORTANT FIX vs the provided main SIU script ***
%  The main script did:  A_echo = A_base * h_mo_mean;  P_rx = A_echo^2;
%  -> that makes RECEIVED POWER scale as h_mo_mean^2, which contradicts
%     Eq.(2)/(8) "Pr = Pr,aligned * h" (linear) and Task 1's explicit
%     instruction "Do not use h_k^2".
%  Fix used here: scale the ECHO AMPLITUDE by sqrt(h_k), so that
%     P_rx = A_base^2 * h_k = P_aligned * h_k   (linear, as required).
%  This is flagged again at the point it is used below, and should be
%  listed under "problems / unclear assumptions" in the write-up.
%
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
%  1. SYSTEM PARAMETERS (Table I, same as the main SIU script)
%% ------------------------------------------------------------------
c  = physconst('LightSpeed');
kb = physconst('Boltzmann');

fc     = 0.85e12;                  % 0.85 THz
lambda = c/fc;
k_wave = 2*pi/lambda;

Gt_dB = 20; Gr_dB = 20;
G_total = 10^((Gt_dB+Gr_dB)/10);

B         = 20e9;
N_fft     = 4096;
N_radar   = 2048;
N_cp      = 1024;
N_symbols = 1024;

delta_f = B/N_fft;
T_sym   = (N_fft+N_cp)/B;
Fs      = B;

% Regime configs -- SAME beam, target size sets SF/SL (matches main script)
w_0_SF = 0.20;  apert_T_SF = 1.5;
w_0_SL = 0.20;  apert_T_SL = 0.04;
apert_R = 0.075;
C_n_sq  = 1e-11;

sigma_rcs_dBsm = 10;
sigma_rcs = 10^(sigma_rcs_dBsm/10);
velocity  = 5;                      % m/s, target radial velocity

regimes      = {'SF','SL'};
w_0_list     = [w_0_SF, w_0_SL];
apert_T_list = [apert_T_SF, apert_T_SL];
waveforms    = {'OFDM','SCFDMA'};

%% ------------------------------------------------------------------
%  2. SENSING / MISALIGNMENT SCENARIO SETUP  (Task 5 requirements)
%% ------------------------------------------------------------------
distance = 15;                      % m  (matches Task 2-4 sigma_s calc)
theta_bm = 30e-3;                   % 30 mrad (matches Task 2-4)
sigma_s  = distance*tan(theta_bm);  % pointing-jitter std, m
Ptx_dBm_fixed = 20;                  % ASSUMPTION: single fixed operating
Ptx_W = 10^((Ptx_dBm_fixed-30)/10); % point (mid of the sweep used in the
                                     % main script), since T5 asks for a
                                     % *sequence of frames*, not a power
                                     % sweep. Change here if a different
                                     % operating power is preferred.

N_frames = 500;                     % ASSUMPTION: 500 sensing frames per
                                     % scenario (>=10,000 was for the
                                     % static histograms in Task 2-4; a
                                     % CPI-by-CPI tracking run of 500
                                     % frames keeps runtime reasonable
                                     % while still showing long
                                     % correlated fade events at rho=0.99)

h_th = 0.10;                        % same threshold as Task 4

rho_list       = [0, 0.9, 0.99];    % 0 stands in for "i.i.d." (rho=0
                                     % collapses AR(1) exactly to i.i.d.)
scenarioNames  = {'Aligned','i.i.d.','AR(1) rho=0.90','AR(1) rho=0.99'};
scenarioTagsSF = {'Aligned','iid','rho0p9','rho0p99'};
nScenario = 4;

fprintf('===============================================\n');
fprintf('Task 5-6 setup\n');
fprintf('===============================================\n');
fprintf('distance=%.1f m | theta_bm=%.0f mrad | sigma_s=%.4f m\n', distance, theta_bm*1000, sigma_s);
fprintf('Ptx = %.1f dBm | N_frames = %d | h_th = %.2f\n\n', Ptx_dBm_fixed, N_frames, h_th);

%% ------------------------------------------------------------------
%  3. ATMOSPHERIC MODEL (same as main script)
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
%  4. RADAR REFERENCE (Zadoff-Chu), built once, reused every frame
%% ------------------------------------------------------------------
n_zc = (0:N_radar-1)';
zc_root = 2;
radar_seq = exp(-1j*pi*zc_root*n_zc.*(n_zc+1)/N_radar);
radar_ref = repmat(radar_seq,1,N_symbols);

k_delay = round(2*distance/c*Fs);
f_D     = 2*velocity/lambda;
t_vec_full = (0:(N_fft+N_cp)*N_symbols-1)'/Fs;
doppler_phase_full = exp(1j*2*pi*f_D*t_vec_full);

rangeRes = c/(2*N_radar*delta_f);            % range bin size
velRes   = lambda/(2*N_symbols*T_sym);       % velocity bin size

%% ------------------------------------------------------------------
%  5. PARALLEL POOL
%% ------------------------------------------------------------------
pool = gcp('nocreate');
if isempty(pool), pool = parpool('local'); end
fprintf('Using parallel pool with %d workers\n\n', pool.NumWorkers);

%% ------------------------------------------------------------------
%  6. MAIN LOOP: regime x waveform x scenario x frame
%% ------------------------------------------------------------------
% Result containers, indexed results.(regime).(waveform).(scenarioTag)
results = struct();

for reg_idx = 1:numel(regimes)
    regime  = regimes{reg_idx};
    w_0     = w_0_list(reg_idx);
    apert_T = apert_T_list(reg_idx);

    fprintf('========== %s CONFIG (w_0=%.0f mm, a_T=%.2f m) ==========\n', regime, w_0*1000, apert_T);

    %% --- beam-waist / w_mo computation, IDENTICAL to main script ---
    rho_0 = (0.55*C_n_sq*k_wave^2*distance)^(-3/5);
    eps_turb = 1 + (2*w_0^2)/(rho_0^2);
    wz = w_0*sqrt(1 + eps_turb*((lambda*distance)/(pi*w_0^2))^2);
    eta = apert_T/wz;
    is_SF = (eta > 1);

    if is_SF
        w0_R = wz;
    else
        w0_R = apert_T*sqrt(2);
    end
    eps_R = 1 + (2*w0_R^2)/(rho_0^2);
    wz_R  = w0_R*sqrt(1 + eps_R*((lambda*distance)/(pi*w0_R^2))^2);
    v_R   = (sqrt(pi)*apert_R)/(sqrt(2)*wz_R);
    weq_R_sq = (wz_R^2*sqrt(pi)*erf(v_R))/(2*v_R*exp(-v_R^2));

    if is_SF
        w_mo_sq = weq_R_sq;
    else
        w_mo_sq = 1/(1/wz^2 + 1/weq_R_sq);
    end
    w_mo = sqrt(w_mo_sq);          % <-- Eq.(1)/(7) beam-width parameter
    fprintf('  eta=%.2f  ->  w_mo = %.4f m\n', eta, w_mo);

    %% --- Task 2/3 generators: build h_k(k) for every scenario --------
    ux = randn(N_frames,1);  uy = randn(N_frames,1);   % shared innovations
    hk_scenario = zeros(N_frames, nScenario);
    rk_scenario = zeros(N_frames, nScenario);

    % 1: aligned
    hk_scenario(:,1) = 1;
    rk_scenario(:,1) = 0;
    % 2: i.i.d. (rho=0 realization uses independent innovations, not the
    %    AR(1) recursion, to match Task 2 exactly)
    xk_iid = sigma_s*randn(N_frames,1);
    yk_iid = sigma_s*randn(N_frames,1);
    rk_scenario(:,2) = sqrt(xk_iid.^2+yk_iid.^2);
    hk_scenario(:,2) = exp(-2*rk_scenario(:,2).^2/w_mo^2);
    % 3 & 4: correlated AR(1), rho = 0.9, 0.99
    for si = 1:2
        rho = rho_list(si+1);
        [xk,yk,rk] = gen_correlated_xy(N_frames, rho, sigma_s, ux, uy);
        rk_scenario(:,si+2) = rk;
        hk_scenario(:,si+2) = exp(-2*rk.^2/w_mo^2);
    end

    %% --- radar equation amplitude (aligned reference) -----------------
    P_radar = Ptx_W*(N_radar/N_fft);
    A_base = sqrt(P_radar)*L_atm_lin* ...
             sqrt((G_total*lambda^2*sigma_rcs)/((4*pi)^3*distance^4));

    for wf_idx = 1:numel(waveforms)
        waveform = waveforms{wf_idx};
        is_SCFDMA = strcmp(waveform,'SCFDMA');

        if is_SCFDMA
            scfdma_input = ifft(radar_ref, N_radar, 1)*sqrt(N_radar);
            tx_symbols = fft(scfdma_input, N_radar, 1)/sqrt(N_radar);
        else
            tx_symbols = radar_ref;
        end
        corr_ref = radar_ref;

        tx_freq = zeros(N_fft, N_symbols);
        tx_freq(1:N_radar,:) = tx_symbols;
        tx_time = ifft(tx_freq, N_fft, 1)*sqrt(N_fft);
        tx_cp = [tx_time(end-N_cp+1:end,:); tx_time];
        tx_serial = tx_cp(:);
        rx_delayed = [zeros(k_delay,1); tx_serial(1:end-k_delay)];
        rx_delayed_doppler = rx_delayed .* doppler_phase_full;

        for sc_idx = 1:nScenario
            hk_vec = hk_scenario(:,sc_idx);

            range_est   = zeros(N_frames,1);
            vel_est     = zeros(N_frames,1);
            correct_det = false(N_frames,1);
            SNR_dB      = zeros(N_frames,1);
            Prx_lin     = zeros(N_frames,1);

            parfor kf = 1:N_frames
                h_k = hk_vec(kf); %#ok<PFBNS>

                % --- Eq.(8): Pr,k = Pr,aligned * h_k  --> LINEAR in h ---
                % Amplitude is scaled by sqrt(h_k) so that POWER scales
                % by h_k (NOT h_k^2). See header note.
                A_echo_k = A_base * sqrt(h_k);

                rx_echo = A_echo_k * rx_delayed_doppler;
                noise = sqrt(P_n_radar/2)*(randn(size(rx_echo))+1j*randn(size(rx_echo)));
                rx_signal = rx_echo + noise;

                rx_cp = reshape(rx_signal, N_fft+N_cp, N_symbols);
                rx_time = rx_cp(N_cp+1:end,:);
                rx_freq = fft(rx_time, N_fft, 1)/sqrt(N_fft);
                rx_radar = rx_freq(1:N_radar,:);

                corr_out = rx_radar .* conj(corr_ref);
                range_profile = ifft(corr_out, N_radar, 1);
                range_doppler = fft(range_profile, N_symbols, 2);

                [~, max_idx] = max(abs(range_doppler(:)));
                [r_idx, d_idx_peak] = ind2sub(size(range_doppler), max_idx);

                r_est = (r_idx-1)*c/(2*N_radar*delta_f);
                if d_idx_peak > N_symbols/2
                    d_shift = d_idx_peak - N_symbols - 1;
                else
                    d_shift = d_idx_peak - 1;
                end
                v_est = d_shift*lambda/(2*N_symbols*T_sym);

                range_est(kf) = r_est;
                vel_est(kf)   = v_est;
                correct_det(kf) = (abs(r_est-distance) <= rangeRes/2) && ...
                                  (abs(v_est-velocity)  <= velRes/2);

                Pr_k = A_echo_k^2;                      % = A_base^2*h_k
                Prx_lin(kf) = Pr_k;
                SNR_dB(kf) = 10*log10(Pr_k/P_n_radar);
            end

            failVec = ~correct_det;
            [avgFail, maxFail] = run_length_stats(failVec);

            tag = scenarioTagsSF{sc_idx};
            results.(regime).(waveform).(tag).h_k        = hk_vec;
            results.(regime).(waveform).(tag).r_k         = rk_scenario(:,sc_idx);
            results.(regime).(waveform).(tag).range_est   = range_est;
            results.(regime).(waveform).(tag).vel_est     = vel_est;
            results.(regime).(waveform).(tag).correct_det = correct_det;
            results.(regime).(waveform).(tag).SNR_dB      = SNR_dB;
            results.(regime).(waveform).(tag).Prx_lin      = Prx_lin;
            results.(regime).(waveform).(tag).avgRangeErr = mean(abs(range_est-distance));
            results.(regime).(waveform).(tag).avgVelErr   = mean(abs(vel_est-velocity));
            results.(regime).(waveform).(tag).Pd          = mean(correct_det);
            results.(regime).(waveform).(tag).avgFailRun  = avgFail;
            results.(regime).(waveform).(tag).maxFailRun  = maxFail;

            fprintf('  [%s,%s,%-14s] Pd=%.3f  avgRangeErr=%.4f m  avgVelErr=%.4f m/s  avgFailRun=%.2f  maxFailRun=%d\n', ...
                regime, waveform, scenarioNames{sc_idx}, mean(correct_det), ...
                mean(abs(range_est-distance)), mean(abs(vel_est-velocity)), avgFail, maxFail);
        end
    end
    fprintf('\n');
end

%% ------------------------------------------------------------------
%  7. TASK 6 -- RESULTS TABLE
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
%  8. TASK 6 -- FIGURES
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
    title(sprintf('Task 6 Fig 1: h_k vs Frame, %s regime', regime));
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
    title(sprintf('Task 6 Fig 2: SNR vs Frame, %s regime, %s', regime, wfShow));
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
            title(sprintf('Task 6 Fig 3: Detection Correct(1)/Incorrect(0) vs Frame, %s, %s', regime, wfShow));
        end
        if sc_idx == nScenario, xlabel('Frame Number, k'); end
    end
    save_figure(fig, fullfile(figDir, sprintf('T6_Fig3_detections_vs_frame_%s', regime)));
end

% --- Fig 4: number of consecutive failed frames vs rho (both regimes) ---
fig = figure('Color','w','Position',[100 100 900 550]); hold on;
barWidth = 0.35;
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
title(sprintf('Task 6 Fig 4: Consecutive Failed Frames vs \\rho (%s)', wfShow));
legend('Location','best'); grid on;
save_figure(fig, fullfile(figDir,'T6_Fig4_consec_fail_vs_rho'));

%% ------------------------------------------------------------------
%  9. SAVE DATA
%% ------------------------------------------------------------------
save(fullfile(datDir,'Task5_6_results.mat'), 'results', 'T6', ...
    'distance','theta_bm','sigma_s','Ptx_dBm_fixed','N_frames','h_th', ...
    'rangeRes','velRes','regimes','waveforms','scenarioNames','scenarioTagsSF');

fprintf('\nFigures -> %s\nTable   -> %s\nData    -> %s\n', figDir, resDir, datDir);
fprintf('============ TASKS 5-6 COMPLETE ============\n');

%% ========================================================================
%  LOCAL FUNCTIONS (Task 2-4 generator + helpers, reused verbatim)
%% ========================================================================

function [xk,yk,rk] = gen_correlated_xy(N, rho, sigma, ux, uy)
    % AR(1) correlated pointing-error generator, Eq.(5)-(6)
    xk = zeros(N,1); yk = zeros(N,1);
    xk(1) = sigma*ux(1); yk(1) = sigma*uy(1);
    for k = 2:N
        xk(k) = rho*xk(k-1) + sqrt(1-rho^2)*sigma*ux(k);
        yk(k) = rho*yk(k-1) + sqrt(1-rho^2)*sigma*uy(k);
    end
    rk = sqrt(xk.^2 + yk.^2);
end

function [avgRun, maxRun] = run_length_stats(boolVec)
    % Average & max length of consecutive TRUE runs in boolVec
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