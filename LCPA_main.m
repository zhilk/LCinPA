%% LCPA MAIN ANALYSIS

% Paper: A latent state perspective on Placebo Analgesia
% Author: Alina Panzel
%
% Last Date of changes: 15.02.26
%
%% Set up & initiate

clearvars; clc; close all;

% Project root (one level up from Scripts/)
projRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projRoot, 'Scripts'));

% Load behavioral data
dataDir = fullfile(projRoot,'LCinPA', 'SubjectData');

csvFiles = {
    fullfile(dataDir, 'LCPA_STUDY1_behavioral.csv')
    fullfile(dataDir, 'LCPA_STUDY2_behavioral.csv')
    fullfile(dataDir, 'LCPA_STUDY3_behavioral.csv')
};
studyLabels = {'Study 1 (Behavioural)', 'Study 2 (Control)', 'Study 3 (fMRI)'};

% Read in all three study tables
behavioralData = cell(3,1);
for s = 1:3
    fprintf('Loading %s ...\n', csvFiles{s});
    behavioralData{s} = readtable(csvFiles{s});
    fprintf('  -> %d rows, %d columns\n', height(behavioralData{s}), width(behavioralData{s}));
end

% Load atlases

%% Behavioral Analysis: Linear Mixed Models

allResults = run_lmm(behavioralData, studyLabels);

% plot behavioral ratings
for s = 1:3
    plot_behavioral(behavioralData{s}, studyLabels{s});
end

%% Behavioral Analysis: Computational Modeling (RW vs BCC vs LCM)

modelResults = run_models(behavioralData, studyLabels);

% Diagnose cause structure (uses stored fits — no re-running)
causeDiag = diagnose_causes(modelResults, behavioralData, studyLabels);

% Plot computational modeling results (pass causeDiag for Panels E & F)
plot_modeling(modelResults, studyLabels, causeDiag);

% Dominant cause analysis: show that top-2 causes per block align with
% conditioning-phase categories despite participants inferring >2 causes
domCause = analyze_dominant_causes(modelResults, behavioralData, studyLabels);

% Plot dominant cause structure (Panel A: cause portrait, Panel B: posterior weight)
plot_dominant_causes(domCause, studyLabels);

%% Neural Analysis: ROI Analysis
%%
