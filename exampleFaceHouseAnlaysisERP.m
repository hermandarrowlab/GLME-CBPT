%Script imports and does basic analysis of example data for the CBPT with [G]LME paper.
%This script shows example of how to analyze Event-Related Potentials (ERPs)
%from individual channels using cluster-based permutation testing (CBPT) with LMEs.
%Sampling rate of the data is 1000 Hz.
%
%Code may take a few minutes to run.
%
%
% Requirements to run:
%   1) Need access to example data named "exampleGlmeFaceHouseAnalysisData.mat"
%   2) clusterBasedGlmeERD.m: main function for cluster-based permutation testing with GLMEs
%   3) findgaps.m: custom function for finding start and end points of clusters
%
%
% Example Data structure:
%   0) Parameter variables
%       a)erpParameters: time window parmaeters
%       b) glmParams: glme parameters
%           i) glmParams.glmeFormulaChannel: GLME formula for equations
%           ii) glmParams.groupModels: GLME group-level formulations
%           iii) glmParams.glmParams.altGlmeFormulas: alternative group-level models
%           iv) glmParams.numShuffs: number of bootstrapping shuffles (i.e. 1000)
%           v) glmParams.alphaIndividual: alpha-value for signficance testing (i.e. 0.04)
%           vi) glmParams.alphaGroup: alpha-value for group testing (i.e. 0.01/3)
%           vii) glmParams.minClusterSize: minimum cluster size (i.e. 55 ms)
%   1) Event related potential (ERP): 1st cell is sig channel, 2nd cell is not-signifciant
%       a) selectedErpData: LFPs (uV) algined to image onset
%       b) selectedErpEventTable: corresponding event table
%
%
% Event Table Structure:
%   *Holds all the information for each channel. Row structure parallels
%   event data row structure.
%   1) eventTable.patientName: original subject ID from downloaded datarepository
%   2) eventTable.blockNumber: image block number (1-3)
%   3) eventTable.eventTime: event time in seconds
%   4) eventTable.eventValue: unique image number from 1 to 100
%   5) eventTable.eventType: image type-1 for house, 2 face
%   6) eventTable.novelRepeat: whether image was novel (0), or repeated (1)
%   7) eventTable.eventNumber: ordinal event number from original event
%   data variable named "stim".
%   8) eventTable.novelRepeatReps: 0-2, number of times image has been viewed
%
%
% Interpretation of Beta Weights/Coefficient Estimates
%   1) The unit of the beta wieght varys across analyses and the type of [G]LME used.
%   2) Because house is event type 1 and face is event type 2, postive beta weights
%   for category (eventType) correspond to face-selective channels whereas
%   negative beta coefficents/weights correspond to house-selective channels.
%   3) Because novel images are encoded as 0 and repeated images as 1,
%   postive beta weights for novelty (novelRepeat) correspond to repeat-selective channels
%   whereas negative beta weights correspond to novelty-selective channels.
%
%
%Written by Seth Konig 7/28/2022
%Lasted updated 1/16/2023


%clear everything
clear, clc, close all
tic

%ranomdly seed so people can reproduce
rng(7282022,'twister')


%---What and Where---%
exampleDataFolder = 'C:\Users\koeni117\Documents\Matlab\CBPT-GLME-Code\';
exampleFileName = 'exampleGlmeFaceHouseAnalysisData';

%specify if you want to use parfor loops to speed up analysis
%use of parfor loops will tend to run faster!
useParfor = 1;%1/true or 0/false


%% Load the Data

%---Load Data---%
load([exampleDataFolder exampleFileName '.mat'],'erpParameters','glmParams',...
    'selectedErpData','selectedErpEventTable')


%---Define Time Window Of Interest---%
%we're going to analyze the full time window of the provided data
%from 200 ms before image onset to 600 ms after image onset
%note image turns off after 400 ms
timeWindowOfInterest = [erpParameters.timeWindow(1) erpParameters.timeWindow(end)];


%--Fix Dummy Variable Encoding for Image Category---%
%technically image type 1 and 2 work just fine, but not the traditional
%method for categorical variables; similarly, use categorical() will turn
%image type into categorical dummy variables. 
%Centering is not usually used with dummy variables!
%Centering of all predictors is necessary for interpretation of intercepts,
%which we do not care about usually. 

%for example ERP data
selectedErpEventTable{1}.eventType = selectedErpEventTable{1}.eventType-1;
selectedErpEventTable{2}.eventType = selectedErpEventTable{2}.eventType-1;


%% Run Example Channel-Level ERP Anlaysis
%---Analyze the ERP Data---%
%run LME analysis on ERP data using channel formula
%Assumes ERPs/LFPs are normally distributed at each time point
%and also normally distributed at the cluster level!
%ERPs/LFPs tend to have ryhtmic/fluctating responses that are bidrectional
%and sometimes multi-phasic!
%Units of LFPs are in microvolts (uV).
disp('Running Example LFP/ERP Analysis...')
erpLme = clusterBasedGlmeERD(glmParams.glmeFormulaChannel,glmParams.erpDataType,glmParams.minClusterSize,...
    erpParameters.timeWindow,selectedErpData,selectedErpEventTable,glmParams.numShuffs,...
    glmParams.alphaIndividual,useParfor,timeWindowOfInterest);


%% Plot Channel-Level Results

%---Plot the ERP Data---%
figure('units','normalized','outerposition',[0 0 1 1]);%full screen

%plot significant channel's average response divided by image type
subplot(2,3,1)
plot(erpParameters.timeWindow,mean(selectedErpData{1}),'k','LineWidth',2); %All images
hold on
plot(erpParameters.timeWindow,mean(selectedErpData{1}(selectedErpEventTable{1}.eventType == 1,:))) %house images
plot(erpParameters.timeWindow,mean(selectedErpData{1}(selectedErpEventTable{1}.eventType == 2,:))) %face images
plot(erpParameters.timeWindow,mean(selectedErpData{1}(selectedErpEventTable{1}.novelRepeat == 0,:))) %novel images
plot(erpParameters.timeWindow,mean(selectedErpData{1}(selectedErpEventTable{1}.novelRepeat == 1,:))) %repeated images
hold off
xlabel('Time from Image Onset (ms)')
ylabel('LFP (uV)')
legend({'All','House','Face','Novel','Repeat'},'location','northWest','FontSize',8);
title('Significant Channel: Average Reponse')

%plot significant channel's beta weights
subplot(2,3,2)
plot(erpParameters.timeWindow,erpLme.betaWeightsTimePoints{1})
hold on
yl = ylim;
plot([erpParameters.timeWindow(1) erpParameters.timeWindow(end)],[0 0],'k--')
plot([0 0],[yl(1) yl(2)],'k--')
if erpLme.significantGlmeFit(1) %has a significant cluster
    %then add shading to indicate significant times
    [~,sigPeriodStart,sigPeriodEnd] = findgaps(find(erpLme.significantGlmeTimes(1,:) > 0)); %find start/end of significant clusters
    sigPeriodStart = sigPeriodStart-erpParameters.tWin2; %correct for image onset time
    sigPeriodEnd = sigPeriodEnd-erpParameters.tWin2;%correct for image onset time
    for sPS = 1:length(sigPeriodStart)
        patch([sigPeriodStart(sPS) sigPeriodEnd(sPS) sigPeriodEnd(sPS) sigPeriodStart(sPS)],[yl(1) yl(1) yl(2) yl(2)],'k','FaceAlpha',0.2)
    end
end
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Beta Weight (uV)')
legend(erpLme.fixedEffectNames,'location','northwest')
title('Significant Channel: Beta Weight Time Course')

%plot significant channel's cluster-level beta weights
subplot(2,3,3)
b = bar(1:erpLme.numFixedEffects,erpLme.sigGlmes{1}{1}.Coefficients.Estimate(2:end),'FaceAlpha',0.4); %ignore intercept
b(1).FaceColor = [.2 .6 .5];
xticks(1:erpLme.numFixedEffects)
xticklabels(erpLme.fixedEffectNames)
hold on
for nFE = 1:erpLme.numFixedEffects
    text(nFE,erpLme.sigGlmes{1}{1}.Coefficients.Estimate(1+nFE)+25,...
        ['\beta = ' num2str(erpLme.sigGlmes{1}{1}.Coefficients.Estimate(1+nFE),3) ...
        ' uV, p = ' num2str(erpLme.sigGlmes{1}{1}.Coefficients.pValue(1+nFE),3)])
end
hold off
ylabel('Beta Weight (uV)')
box off
ylim([0 850])
title('Significant Channel: Cluster-Level Beta Weights')


%non-significant channel's average response divided by image type
subplot(2,3,4)
plot(erpParameters.timeWindow,mean(selectedErpData{2}),'k','LineWidth',2); %All images
hold on
plot(erpParameters.timeWindow,mean(selectedErpData{2}(selectedErpEventTable{2}.eventType == 1,:))) %house images
plot(erpParameters.timeWindow,mean(selectedErpData{2}(selectedErpEventTable{2}.eventType == 2,:))) %face images
plot(erpParameters.timeWindow,mean(selectedErpData{2}(selectedErpEventTable{2}.novelRepeat == 0,:))) %novel images
plot(erpParameters.timeWindow,mean(selectedErpData{2}(selectedErpEventTable{2}.novelRepeat == 1,:))) %repeated images
hold off
xlabel('Time from Image Onset (ms)')
ylabel('LFP (uV)')
legend({'All','House','Face','Novel','Repeat'},'location','northWest','FontSize',8);
title('Non-significant Channel: Average Reponse')

%plot non significant channel's beta weights
subplot(2,3,5)
plot(erpParameters.timeWindow,erpLme.betaWeightsTimePoints{2})
hold on
yl = ylim;
plot([erpParameters.timeWindow(1) erpParameters.timeWindow(end)],[0 0],'k--')
plot([0 0],[yl(1) yl(2)],'k--')
if erpLme.significantGlmeFit(2) %has a significant cluster
    %then add shading to indicate significant times
    [~,sigPeriodStart,sigPeriodEnd] = findgaps(find(erpLme.significantGlmeTimes(2,:) > 0)); %find start/end of significant clusters
    sigPeriodStart = sigPeriodStart-erpParameters.tWin2; %correct for image onset time
    sigPeriodEnd = sigPeriodEnd-erpParameters.tWin2;%correct for image onset time
    for sPS = 1:length(sigPeriodStart)
        patch([sigPeriodStart(sPS) sigPeriodEnd(sPS) sigPeriodEnd(sPS) sigPeriodStart(sPS)],[yl(1) yl(1) yl(2) yl(2)],'k','FaceAlpha',0.2)
    end
end
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Beta Weight (uV)')
legend(erpLme.fixedEffectNames,'location','northwest')
title('Non-significant Channel: Beta Weight Time Course')

%plot non significant channel's p-values over time
potentialSigTimes = any(erpLme.sigIndividualTimePoints{2} < glmParams.alphaIndividual/erpLme.numFixedEffects,2);%find potential significant times
[~,potentialSigStart,potentialSigEnd] = findgaps(find(potentialSigTimes == 1)); %find potential clusters
largestPotentialCluster = max(potentialSigEnd-potentialSigStart+1);
subplot(2,3,6)
plot(erpParameters.timeWindow,log(erpLme.sigIndividualTimePoints{2}))
hold on
plot([erpParameters.timeWindow(1) erpParameters.timeWindow(end)],[log(glmParams.alphaIndividual/erpLme.numFixedEffects) log(glmParams.alphaIndividual/erpLme.numFixedEffects)],'k--')
hold off
text(-180,-7,['Largest Potential Cluster was only ' num2str(largestPotentialCluster) ' ms'],'FontSize',10,'HorizontalAlignment','left')
xlabel('Time from Image Onset (ms)')
ylabel('Log(p-value)')
legend(erpLme.fixedEffectNames,'location','southwest')
title('Non-significant Channel: Beta Weight P-values')

sgtitle('Example ERP Analysis using Cluster-Based Statistics with a LME')

%%
toc