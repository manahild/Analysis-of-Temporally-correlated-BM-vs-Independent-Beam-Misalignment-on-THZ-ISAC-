%% ========================================================================
%  CORRELATED BEAM MISALIGNMENT IN THz SENSING
%  PURE 2026 Program -- COMBINED SCRIPT: Tasks 2, 3, 4, 5, 6
%
%  Everything runs top-to-bottom in ONE script and ONE workspace -- no
%  .mat loading between tasks, so there is nothing to "run out of order."
%  Just press Run.
%
%  Equations used throughout (assignment sheet numbering):
%    Eq.(1)/(7): h_k = exp(-2*r_k^2 / w_mo^2)
%    Eq.(3)-(4): i.i.d. case, x_k~N(0,sigma^2), y_k~N(0,sigma^2), r_k=sqrt(x^2+y^2)
%    Eq.(5)-(6): correlated case (AR(1)):
%         x_k = rho*x_{k-1} + sqrt(1-rho^2)*sigma*ux,k
%         y_k = rho*y_{k-1} + sqrt(1-rho^2)*sigma*uy,k
%    Eq.(2)/(8): Pr,k = Pr,aligned * h_k   (LINEAR in h, NOT h^2)
%
%  Author: Manahil Ahmad (PURE 2026)
% =========================================================================

clear;
clc;
close all;

rng(1);                       % Single global seed -- reproducible top-to-bottom

%% ========================================================================
%  0. FOLDER SETUP (absolute paths, robust to current-folder drift)
%% ========================================================================

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

figDir = fullfile(scriptDir, 'Figures');
resDir = fullfile(scriptDir, 'Results');
datDir = fullfile(scriptDir, 'Data');

folderList = {figDir, resDir, datDir};
for f = 1:numel(folderList)
    if ~exist(folderList{f}, 'dir')
        [ok, msg] = mkdir(folderList{f});
        if ~ok
            error('Could not create folder "%s": %s', folderList{f}, msg);
        end
    end
end

fprintf('Script folder     = %s\n', scriptDir);
fprintf('Figures folder    = %s\n', figDir);
fprintf('Results folder    = %s\n', resDir);
fprintf('Data folder       = %s\n\n', datDir);

%% ========================================================================
%  0. SYSTEM / BEAM PARAMETERS (shared across all tasks)
%% ========================================================================

N        = 10000;             % samples/frames for Tasks 2-4 statistics
distance = 15;                 % m
theta_bm = 30e-3;              % 30 mrad
sigma_s  = distance*tan(theta_bm);   % pointing jitter std, m

w_mo_SF = 0.2158;              % m, Spot-Filling regime (from main THz-ISAC sim)
w_mo_SL = 0.0721;              % m, Spot-Limited regime (from main THz-ISAC sim)

rho_list = [0, 0.9, 0.99];     % correlation coefficients to test
h_th     = 0.10;               % detection/outage threshold (fixed for all comparisons)

fprintf('===============================================\n');
fprintf('Shared parameters\n');
fprintf('===============================================\n');
fprintf('Distance    = %.1f m | Beam angle = %.0f mrad | sigma_s = %.4f m\n', ...
    distance, theta_bm*1000, sigma_s);
fprintf('w_mo SF     = %.4f m | w_mo SL    = %.4f m\n', w_mo_SF, w_mo_SL);
fprintf('rho values  = %s | h_th = %.2f | N = %d\n\n', mat2str(rho_list), h_th, N);

%% ========================================================================
%  TASK 2 -- INDEPENDENT (i.i.d.) MISALIGNMENT SAMPLES
%  Eq.(3): xk~N(0,sigma^2), yk~N(0,sigma^2); Eq.(4): rk; Eq.(1): hk
% =========================================================================

fprintf('=== TASK 2: i.i.d. reference case ===\n');

xk_iid = sigma_s * randn(N,1);
yk_iid = sigma_s * randn(N,1);
rk_iid = sqrt(xk_iid.^2 + yk_iid.^2);
hk_iid_SF = exp(-2*rk_iid.^2 / w_mo_SF^2);
hk_iid_SL = exp(-2*rk_iid.^2 / w_mo_SL^2);

r_theory = linspace(0, max(rk_iid), 500);
pdf_rayleigh = (r_theory./sigma_s.^2).*exp(-(r_theory.^2)./(2*sigma_s.^2));

% Figure T2-1: x_k
fig = figure('Color','w','Position',[100 100 900 550]);
plot(xk_iid,'b'); xlabel('Sample Number'); ylabel('x_k (m)');
title('Task 2: Independent Horizontal Pointing Error, x_k'); grid on;
save_figure(fig, fullfile(figDir,'T2_Fig1_xk'));

% Figure T2-2: r_k
fig = figure('Color','w','Position',[100 100 900 550]);
plot(rk_iid,'k'); xlabel('Sample Number'); ylabel('r_k (m)');
title('Task 2: Radial Pointing Displacement, r_k'); grid on;
save_figure(fig, fullfile(figDir,'T2_Fig2_rk'));

% Figure T2-3: histogram + Rayleigh theory
fig = figure('Color','w','Position',[100 100 900 550]);
histogram(rk_iid,50,'Normalization','pdf'); hold on;
plot(r_theory,pdf_rayleigh,'r','LineWidth',2);
xlabel('r_k (m)'); ylabel('Probability Density');
title('Task 2: Rayleigh Distribution of Pointing Error, r_k');
legend('Simulation','Theory','Location','best'); grid on;
save_figure(fig, fullfile(figDir,'T2_Fig3_RayleighPDF'));

% Figure T2-4: h_k, SF
fig = figure('Color','w','Position',[100 100 900 550]);
plot(hk_iid_SF,'m'); xlabel('Sample Number'); ylabel('h_k');
title(sprintf('Task 2: Beam Coefficient, SF Regime (w_{mo}=%.4f m)', w_mo_SF));
ylim([0 1]); grid on;
save_figure(fig, fullfile(figDir,'T2_Fig4_hk_SF'));

% Figure T2-5: h_k, SL
fig = figure('Color','w','Position',[100 100 900 550]);
plot(hk_iid_SL,'Color',[0.20 0.55 0.20]); xlabel('Sample Number'); ylabel('h_k');
title(sprintf('Task 2: Beam Coefficient, SL Regime (w_{mo}=%.4f m)', w_mo_SL));
ylim([0 1]); grid on;
save_figure(fig, fullfile(figDir,'T2_Fig5_hk_SL'));

fprintf('Mean(x_k)=%.5f  Mean(r_k)=%.5f  Mean(h_k)SF=%.5f  Mean(h_k)SL=%.5f\n\n', ...
    mean(xk_iid), mean(rk_iid), mean(hk_iid_SF), mean(hk_iid_SL));

%% ========================================================================
%  TASK 3 -- CORRELATED (AR(1)) MISALIGNMENT SAMPLES
%  Eq.(5)-(6): AR(1) recursion; Eq.(7): rk, hk
% =========================================================================

fprintf('=== TASK 3: Correlated (AR(1)) cases, rho = %s ===\n', mat2str(rho_list));

% Shared innovations across all rho -- isolates rho as the ONLY variable
ux = randn(N,1);
uy = randn(N,1);

nRho = numel(rho_list);
xk_all   = zeros(N,nRho);
yk_all   = zeros(N,nRho);
rk_all   = zeros(N,nRho);
hkSF_all = zeros(N,nRho);
hkSL_all = zeros(N,nRho);

for ri = 1:nRho
    rho = rho_list(ri);
    [xk, yk, rk] = gen_correlated_xy(N, rho, sigma_s, ux, uy);
    xk_all(:,ri) = xk; yk_all(:,ri) = yk; rk_all(:,ri) = rk;
    hkSF_all(:,ri) = exp(-2*rk.^2/w_mo_SF^2);
    hkSL_all(:,ri) = exp(-2*rk.^2/w_mo_SL^2);
end

colors = lines(nRho);

% Figure T3-1: r_k vs frame, all rho
fig = figure('Color','w','Position',[100 100 1000 550]); hold on;
for ri = 1:nRho
    plot(rk_all(:,ri),'Color',colors(ri,:),'DisplayName',sprintf('\\rho=%.2f',rho_list(ri)));
end
xlabel('Frame Number, k'); ylabel('r_k (m)');
title('Task 3: r_k vs Frame Number, by \rho'); legend('Location','best'); grid on;
save_figure(fig, fullfile(figDir,'T3_Fig1_rk_vs_frame_allrho'));

% Figure T3-2: h_k SF vs frame, all rho
fig = figure('Color','w','Position',[100 100 1000 550]); hold on;
for ri = 1:nRho
    plot(hkSF_all(:,ri),'Color',colors(ri,:),'DisplayName',sprintf('\\rho=%.2f',rho_list(ri)));
end
xlabel('Frame Number, k'); ylabel('h_k');
title(sprintf('Task 3: h_k vs Frame, SF Regime (w_{mo}=%.4f m)', w_mo_SF));
ylim([0 1]); legend('Location','best'); grid on;
save_figure(fig, fullfile(figDir,'T3_Fig2_hk_SF_vs_frame_allrho'));

% Figure T3-3: h_k SL vs frame, all rho
fig = figure('Color','w','Position',[100 100 1000 550]); hold on;
for ri = 1:nRho
    plot(hkSL_all(:,ri),'Color',colors(ri,:),'DisplayName',sprintf('\\rho=%.2f',rho_list(ri)));
end
xlabel('Frame Number, k'); ylabel('h_k');
title(sprintf('Task 3: h_k vs Frame, SL Regime (w_{mo}=%.4f m)', w_mo_SL));
ylim([0 1]); legend('Location','best'); grid on;
save_figure(fig, fullfile(figDir,'T3_Fig3_hk_SL_vs_frame_allrho'));

% Figure T3-4: combined histogram comparison
fig = figure('Color','w','Position',[100 100 900 550]); hold on;
for ri = 1:nRho
    histogram(rk_all(:,ri),50,'Normalization','pdf', ...
        'DisplayName',sprintf('\\rho=%.2f',rho_list(ri)), ...
        'FaceAlpha',0.35,'FaceColor',colors(ri,:));
end
plot(r_theory,pdf_rayleigh,'k--','LineWidth',2,'DisplayName','Rayleigh theory');
xlabel('r_k (m)'); ylabel('Probability Density');
title('Task 3: Histogram of r_k, all \rho (marginal shape unchanged)');
legend('Location','best'); grid on;
save_figure(fig, fullfile(figDir,'T3_Fig4_rk_histogram_compare'));

% Figures T3-5,6,7: individual histograms per rho
histNames = {'T3_Fig5_rk_histogram_rho0','T3_Fig6_rk_histogram_rho0p9','T3_Fig7_rk_histogram_rho0p99'};
for ri = 1:nRho
    fig = figure('Color','w','Position',[100 100 900 550]);
    histogram(rk_all(:,ri),50,'Normalization','pdf','FaceColor',colors(ri,:),'FaceAlpha',0.6);
    hold on; plot(r_theory,pdf_rayleigh,'r','LineWidth',2);
    xlabel('r_k (m)'); ylabel('Probability Density');
    title(sprintf('Task 3: Histogram of r_k, \\rho=%.2f', rho_list(ri)));
    legend('Simulation','Rayleigh theory','Location','best'); grid on;
    save_figure(fig, fullfile(figDir, histNames{ri}));
end

for ri = 1:nRho
    xk_ri = xk_all(:,ri);
    ac = corrcoef(xk_ri(1:end-1), xk_ri(2:end));
    fprintf('rho=%.2f: measured lag-1 autocorr=%.4f | Mean(h_k)SF=%.5f | Mean(h_k)SL=%.5f\n', ...
        rho_list(ri), ac(1,2), mean(hkSF_all(:,ri)), mean(hkSL_all(:,ri)));
end
fprintf('\n');

%% ========================================================================
%  TASK 4 -- STATISTICAL COMPARISON TABLE
%  Uses the SAME in-memory variables from Task 2 & 3 -- no file I/O needed
% =========================================================================

fprintf('=== TASK 4: Statistical comparison table (h_th = %.2f) ===\n', h_th);

caseLabels = [{'i.i.d. (Task 2)'}, arrayfun(@(r) sprintf('AR(1), rho=%.2f',r), rho_list, 'UniformOutput', false)];
hkSF_cases = [hk_iid_SF, hkSF_all];
hkSL_cases = [hk_iid_SL, hkSL_all];
nCases = numel(caseLabels);

regimeNames = {'SF','SL'};
regimeData  = {hkSF_cases, hkSL_cases};
w_mo_vals   = [w_mo_SF, w_mo_SL];

T4 = table();
for regIdx = 1:2
    hkMat = regimeData{regIdx};
    for c = 1:nCases
        hk = hkMat(:,c);
        below = hk < h_th;
        [avgRun, maxRun] = run_length_stats(below);
        newRow = table({regimeNames{regIdx}}, {caseLabels{c}}, w_mo_vals(regIdx), ...
            mean(hk), min(hk), 100*mean(below), avgRun, maxRun, ...
            'VariableNames', {'Regime','Case','w_mo_m','AvgH','MinH', ...
            'PctFramesBelowThresh','AvgConsecRunLen','MaxConsecRunLen'});
        T4 = [T4; newRow]; %#ok<AGROW>
    end
end

disp(T4);
writetable(T4, fullfile(resDir,'Task4_comparison_table.csv'));

fprintf('\n');



%% ========================================================================
%  SAVE ALL DATA (Tasks 2-4)
%% ========================================================================

save(fullfile(datDir,'Combined_Task2to4_data.mat'), ...
    'xk_iid','yk_iid','rk_iid', ...
    'hk_iid_SF','hk_iid_SL', ...
    'rho_list', ...
    'xk_all','yk_all','rk_all', ...
    'hkSF_all','hkSL_all', ...
    'w_mo_SF','w_mo_SL', ...
    'sigma_s','h_th','T4');

fprintf('\nAll figures saved to  %s\n', figDir);
fprintf('All tables saved to   %s\n', resDir);
fprintf('All data saved to     %s\n', fullfile(datDir,'Combined_Task2to4_data.mat'));
fprintf('\n============ TASKS 2-4 COMPLETE ============\n');

%% ========================================================================
%  LOCAL FUNCTIONS
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
    % Saves a figure as FIG, PNG, and PDF at publication quality.
    % baseName should already be an ABSOLUTE path with no extension.
    % Uses exportgraphics (robust, R2020a+) rather than legacy print().
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