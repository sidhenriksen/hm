% Full filtering pipeline:

clear all; close all; clc

%fileName = '/b/data/lem/M303/lemM303.rds.dxXmixacXdw.Cells.mat';
fileName = '/b/data/lem/M302/lemM302.rds.dxXmixacXdd.Cells.mat';

load(fileName);

% For some reason, the suppress plot thing doesn't work so we will get a
% tonne of plots. Just close these when we're done.
P = run_ANOVA(AllExpt);
close all; 

neuronsAfterANOVA = find(P < 0.1);

if exist('cExpt','var');
    array = 0;
else
    array = 1;
end

% We need different procedures for laminar recordings and multi-electrode
% arrays
if array
    [Dx,MeanResponse,SDResponse,SEM] = AnalyzeNewLaminar(AllExpt);
    Expt = AllExpt;
else
    [Dx,MeanResponse,SDResponse,SEM] = AnalyzeLaminar(cExpt);
    Expt = cExpt;
end

%spikeCountM(dx,dd,corr,dw,neuron) = currentMean;

dw=1;
MaxResp = squeeze(max(MeanResponse(:,2,1,dw,neuronsAfterANOVA)));

neuronsAfterMax = neuronsAfterANOVA(MaxResp > 4);

% Inspection of the tuning curves for lemM303 reveals that 9 and 11 are duplicates of 
% other cells. We remove these because there are fewer spikes in these than
% in the other cells.
%neuronsAfterMax = setdiff(neuronsAfterMax,[9,11]);

