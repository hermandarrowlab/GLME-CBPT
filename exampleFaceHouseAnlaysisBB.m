%Script imports and does basic analysis of example data for the CBPT with [G]LME paper.
%This script shows example of how to analyze Broadband (BB) power signals
%from individual channels using cluster-based permutation testing (CBPT) with LMEs.
%Code also shows how to anlayze group-level broadband power signals!
%
%Sampling rate of the data is 1000 Hz.
%
%Code takes at least few minutes to run since anlayzes data from individual
%channels and group-level data.
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
%   2) Broadband (BB): 1st cell is sig channel, 2nd cell is not-signifciant
%       a) selectedBBData: broadband power (a.u.) aligned to image onset.
%       b) selectedBBEventTable: corresponding event table
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
    'selectedBBData','selectedBBEventTable',...
    'groupedFaceBBData','groupedFaceBBEventTable')


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

%for example broadband power data
selectedBBEventTable{1}.eventType = selectedBBEventTable{1}.eventType-1;
selectedBBEventTable{2}.eventType = selectedBBEventTable{2}.eventType-1;

%for example broadband power group data
groupedFaceBBEventTable{1}.eventType = groupedFaceBBEventTable{1}.eventType-1;


%% Run Channel-Level Example BroadBand Power (BB) Anlaysis
%---Analyze the BB Data---%
%run LME analysis on BB data using channel formula
%Assumes broadband power is normally distributed at each time point
%and also normally distributed at the cluster level!
%Broadband repsonse tend to be unidirectional as channels tend to have
%increased broadband power when they're selective!
%Units of braodband power are arbitrary (a.u.) as they are calculated using
%PCA of the broadband power spectrum and then z-scored.
disp('Running Example Broadband Power (BB) Analysis...')
bbLme = clusterBasedGlmeERD(glmParams.glmeFormulaChannel,glmParams.bbDataType,glmParams.minClusterSize,...
    erpParameters.timeWindow,selectedBBData,selectedBBEventTable,glmParams.numShuffs,...
    glmParams.alphaIndividual,useParfor,timeWindowOfInterest);


%% Plot Channel-Level Results
%---Plot the BB Data---%
figure('units','normalized','outerPosition',[0 0 1 1]);%full screen

%plot significant channel's average response divided by image type
subplot(2,3,1)
plot(erpParameters.timeWindow,mean(selectedBBData{1}),'k','LineWidth',2); %All images
hold on
plot(erpParameters.timeWindow,mean(selectedBBData{1}(selectedBBEventTable{1}.eventType == 1,:))) %house images
plot(erpParameters.timeWindow,mean(selectedBBData{1}(selectedBBEventTable{1}.eventType == 2,:))) %face images
plot(erpParameters.timeWindow,mean(selectedBBData{1}(selectedBBEventTable{1}.novelRepeat == 0,:))) %novel images
plot(erpParameters.timeWindow,mean(selectedBBData{1}(selectedBBEventTable{1}.novelRepeat == 1,:))) %repeated images
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Broadband Power (a.u.)')
legend({'All','House','Face','Novel','Repeat'},'location','northWest','FontSize',8);
title('Significant Channel: Average Reponse')

%plot significant channel's beta weights
subplot(2,3,2)
plot(erpParameters.timeWindow,bbLme.betaWeightsTimePoints{1})
hold on
yl = ylim;
plot([erpParameters.timeWindow(1) erpParameters.timeWindow(end)],[0 0],'k--')
plot([0 0],[yl(1) yl(2)],'k--')
if bbLme.significantGlmeFit(1) %has a significant cluster
    %then add shading to indicate significant times
    [~,sigPeriodStart,sigPeriodEnd] = findgaps(find(bbLme.significantGlmeTimes(1,:) > 0)); %find start/end of significant clusters
    sigPeriodStart = sigPeriodStart-erpParameters.tWin2; %correct for image onset time
    sigPeriodEnd = sigPeriodEnd-erpParameters.tWin2;%correct for image onset time
    for sPS = 1:length(sigPeriodStart)
        patch([sigPeriodStart(sPS) sigPeriodEnd(sPS) sigPeriodEnd(sPS) sigPeriodStart(sPS)],[yl(1) yl(1) yl(2) yl(2)],'k','FaceAlpha',0.2)
    end
end
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Beta Weight (a.u.)')
legend(bbLme.fixedEffectNames,'location','northwest')
title('Significant Channel: Beta Weight Time Course')

%plot significant channel's cluster-level beta weights
subplot(2,3,3)
b = bar(1:bbLme.numFixedEffects,[bbLme.sigGlmes{1}{1}.Coefficients.Estimate(2:end) bbLme.sigGlmes{1}{2}.Coefficients.Estimate(2:end)],'FaceAlpha',0.5); %ignore intercept
b(1).FaceColor = [.2 .6 .5];
b(2).FaceColor = [0.8 0.2 0.5];
xticks(1:bbLme.numFixedEffects)
xticklabels(bbLme.fixedEffectNames)
hold on
for clust = 1:2
    for nFE = 1:bbLme.numFixedEffects
        if bbLme.sigGlmes{1}{clust}.Coefficients.Estimate(1+nFE) > 0
            offset = 0.1;
        else
            offset = -0.1;
        end
        text(nFE,bbLme.sigGlmes{1}{clust}.Coefficients.Estimate(1+nFE)+offset,...
            ['\beta = ' num2str(bbLme.sigGlmes{1}{clust}.Coefficients.Estimate(1+nFE),3) ...
            ' a.u., p = ' num2str(bbLme.sigGlmes{1}{clust}.Coefficients.pValue(1+nFE),3)])
    end
end
hold off
ylim([-1 1])
ylabel('Beta Weight (a.u.)')
box off
legend({'Cluster 1','Cluster 2'},'location','northwest')
title('Significant Channel: Cluster-Level Beta Weights')


%non-significant channel's average response divided by image type
subplot(2,3,4)
plot(erpParameters.timeWindow,mean(selectedBBData{2}),'k','LineWidth',2); %All images
hold on
plot(erpParameters.timeWindow,mean(selectedBBData{2}(selectedBBEventTable{2}.eventType == 1,:))) %house images
plot(erpParameters.timeWindow,mean(selectedBBData{2}(selectedBBEventTable{2}.eventType == 2,:))) %face images
plot(erpParameters.timeWindow,mean(selectedBBData{2}(selectedBBEventTable{2}.novelRepeat == 0,:))) %novel images
plot(erpParameters.timeWindow,mean(selectedBBData{2}(selectedBBEventTable{2}.novelRepeat == 1,:))) %repeated images
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Broadband Power (a.u.)')
legend({'All','House','Face','Novel','Repeat'},'location','northWest','FontSize',8);
title('Non-significant Channel: Average Reponse')

%plot non significant channel's beta weights
subplot(2,3,5)
plot(erpParameters.timeWindow,bbLme.betaWeightsTimePoints{2})
hold on
yl = ylim;
plot([erpParameters.timeWindow(1) erpParameters.timeWindow(end)],[0 0],'k--')
plot([0 0],[yl(1) yl(2)],'k--')
if bbLme.significantGlmeFit(2) %has a significant cluster
    %then add shading to indicate significant times
    [~,sigPeriodStart,sigPeriodEnd] = findgaps(find(bbLme.significantGlmeTimes(2,:) > 0)); %find start/end of significant clusters
    sigPeriodStart = sigPeriodStart-erpParameters.tWin2; %correct for image onset time
    sigPeriodEnd = sigPeriodEnd-erpParameters.tWin2;%correct for image onset time
    for sPS = 1:length(sigPeriodStart)
        patch([sigPeriodStart(sPS) sigPeriodEnd(sPS) sigPeriodEnd(sPS) sigPeriodStart(sPS)],[yl(1) yl(1) yl(2) yl(2)],'k','FaceAlpha',0.2)
    end
end
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Beta Weight (a.u.)')
legend(bbLme.fixedEffectNames,'location','northwest')
title('Non-significant Channel: Beta Weight Time Course')

%plot non significant channel's p-values over time
potentialSigTimes = any(bbLme.sigIndividualTimePoints{2} < glmParams.alphaIndividual/bbLme.numFixedEffects,2);%find potential significant times
[~,potentialSigStart,potentialSigEnd] = findgaps(find(potentialSigTimes == 1)); %find potential clusters
largestPotentialCluster = max(potentialSigEnd-potentialSigStart+1);
subplot(2,3,6)
plot(erpParameters.timeWindow,log(bbLme.sigIndividualTimePoints{2}))
hold on
plot([erpParameters.timeWindow(1) erpParameters.timeWindow(end)],[log(glmParams.alphaIndividual/bbLme.numFixedEffects) log(glmParams.alphaIndividual/bbLme.numFixedEffects)],'k--')
hold off
text(-180,-7,['Largest Potential Cluster was only ' num2str(largestPotentialCluster) ' ms'],'FontSize',10,'HorizontalAlignment','left')
xlabel('Time from Image Onset (ms)')
ylabel('Log(p-value)')
legend(bbLme.fixedEffectNames,'location','southwest')
title('Non-significant Channel: Beta Weight P-Values')

sgtitle('Example BB Analysis using Cluster-Based Statistics with a LME')



%% Example Grouped-Face-Selective BroadBand Power (BB) Data Anlaysis
%run LME analysis on BB data using group-level formulas one at a time so we
%can compare which model is best using a log-liklihood ratio test.
%Assumes broadband power is normally distributed at each time point
%and also normally distributed at the cluster level!
%Broadband repsonse tend to be unidirectional and since all face-selective
%channels then the category selectivity should also be in the same
%direction.
%Group-level models also allow us to look at the random effects across
%channels.
disp('Running Example Grouped-Face-Selective Broadband Power (BB) Analysis...')

%---Analyze the Group-level BB Data---%
%run through each group model
groupLmes = cell(1,length(glmParams.groupModels));
for gM = 1:length(glmParams.groupModels) %run through each group model
    groupLmes{gM} = clusterBasedGlmeERD(glmParams.groupModels{gM},glmParams.bbDataType,glmParams.minClusterSize,...
        erpParameters.timeWindow,groupedFaceBBData,groupedFaceBBEventTable,glmParams.numShuffs,...
        glmParams.alphaGroup,useParfor,timeWindowOfInterest);
end

%compare models using log-liklihood ratio test
clustersToUse = [1 1 1 1];%here there's only 1 significant clsuter to compare so this is easy
isMoreComplicatedModelBetter = NaN(1,length(glmParams.groupModels));
groupGlmRSquared = NaN(1,length(glmParams.groupModels)); %adjusted
groupGlmAIC = NaN(1,length(glmParams.groupModels));
for gM = 1:length(glmParams.groupModels) %run through each group model
    groupGlmRSquared(gM) = groupLmes{gM}.sigGlmes{1}{clustersToUse(gM)}.Rsquared.Adjusted;
    groupGlmAIC(gM) =  groupLmes{gM}.sigGlmes{1}{clustersToUse(gM)}.ModelCriterion.AIC;
    if gM > 1
        if groupLmes{gM-1}.sigGlmes{1}{clustersToUse(gM)}.LogLikelihood < groupLmes{gM}.sigGlmes{1}{clustersToUse(gM)}.LogLikelihood
            llrCompare = compare(groupLmes{gM-1}.sigGlmes{1}{clustersToUse(gM)},groupLmes{gM}.sigGlmes{1}{clustersToUse(gM)});
            if llrCompare.pValue(2) < glmParams.alphaGroup
                isMoreComplicatedModelBetter(gM) = 1;
            else
                isMoreComplicatedModelBetter(gM) = 0;
            end
        else %conisder the model not better
            isMoreComplicatedModelBetter(gM) = 0;
        end
    end
end


%get the best group model, should be 3, AIC alone would say this too
bestGroupModel = find(isMoreComplicatedModelBetter);
bestGroupModel = bestGroupModel(end);
disp(['Best Group Model is '  glmParams.groupModelNames{bestGroupModel}])


%---Test alternative LME formulas---%
%run through each group model
altGroupLmes = cell(1,length(glmParams.altGlmeFormulas));
for altGM = 1:length(glmParams.altGlmeFormulas) %run through each alternative group model
    altGroupLmes{altGM} = clusterBasedGlmeERD(glmParams.altGlmeFormulas{altGM},glmParams.bbDataType,glmParams.minClusterSize,...
        erpParameters.timeWindow,groupedFaceBBData,groupedFaceBBEventTable,glmParams.numShuffs,...
        glmParams.alphaGroup,useParfor,timeWindowOfInterest);
end

%compare alternative using log-liklihood ratio test
altClustersToUse = [1 1 1];%here there's only 1 significant clsuter to compare so this is easy
altIsMoreComplicatedModelBetter = NaN(1,length(glmParams.altGlmeFormulas));
altGroupGlmRSquared = NaN(1,length(glmParams.altGlmeFormulas)); %adjusted
altGroupGlmAIC = NaN(1,length(glmParams.altGlmeFormulas));
for altGM = 1:length(glmParams.altGlmeFormulas) %run through each group model
    altGroupGlmRSquared(altGM) = altGroupLmes{altGM}.sigGlmes{1}{altClustersToUse(altGM)}.Rsquared.Adjusted;
    altGroupGlmAIC(altGM) =  altGroupLmes{altGM}.sigGlmes{1}{altClustersToUse(altGM)}.ModelCriterion.AIC;
    if altGM > 1
        if altGroupLmes{altGM-1}.sigGlmes{1}{altClustersToUse(altGM)}.LogLikelihood < altGroupLmes{altGM}.sigGlmes{1}{altClustersToUse(altGM)}.LogLikelihood
            llrCompare = compare(altGroupLmes{altGM-1}.sigGlmes{1}{altClustersToUse(altGM)},altGroupLmes{altGM}.sigGlmes{1}{altClustersToUse(altGM)});
            if llrCompare.pValue(2) < glmParams.alphaGroup
                altIsMoreComplicatedModelBetter(altGM) = 1;
            else
                altIsMoreComplicatedModelBetter(altGM) = 0;
            end
        else %alt model is not better
            altIsMoreComplicatedModelBetter(altGM) = 0;
        end
    end
end

%make sure interaction + imag Id is actually better than basic
if altGroupLmes{1}.sigGlmes{1}{altClustersToUse(1)}.LogLikelihood < altGroupLmes{3}.sigGlmes{1}{altClustersToUse(3)}.LogLikelihood
    llrCompare = compare(altGroupLmes{1}.sigGlmes{1}{altClustersToUse(1)},altGroupLmes{3}.sigGlmes{1}{altClustersToUse(3)});
    if llrCompare.pValue(2) > glmParams.alphaGroup
        altIsMoreComplicatedModelBetter(3) = 0;
    end
else
    altIsMoreComplicatedModelBetter(3) = 0;
end
 

%get the best alternative group model, should be 3, AIC alone would say this too
bestAltGroupModel = find(altIsMoreComplicatedModelBetter);
bestAltGroupModel = bestAltGroupModel(end);
disp(['Best <strong> Alternative </strong> Group Model is '  glmParams.altGlmeNames{bestAltGroupModel}])


%% Plot Group-Level Results
%---Plot the Group-Level Face BB Data---%
figure('units','normalized','outerPosition',[0 0 1 1]);%full screen

subplot(2,3,1)
plot(erpParameters.timeWindow,mean(groupedFaceBBData{1}),'k','LineWidth',2); %All images
hold on
plot(erpParameters.timeWindow,mean(groupedFaceBBData{1}(groupedFaceBBEventTable{1}.eventType == 1,:))) %house images
plot(erpParameters.timeWindow,mean(groupedFaceBBData{1}(groupedFaceBBEventTable{1}.eventType == 2,:))) %face images
plot(erpParameters.timeWindow,mean(groupedFaceBBData{1}(groupedFaceBBEventTable{1}.novelRepeat == 0,:))) %novel images
plot(erpParameters.timeWindow,mean(groupedFaceBBData{1}(groupedFaceBBEventTable{1}.novelRepeat == 1,:))) %repeated images
hold off
xlabel('Time from Image Onset (ms)')
ylabel('Broadband Power (a.u.)')
legend({'All','House','Face','Novel','Repeat'},'location','northWest','FontSize',8);
title('Average Reponse')

subplot(2,3,2)
plot(erpParameters.timeWindow,groupLmes{bestGroupModel}.betaWeightsTimePoints{1})
hold on
yl = ylim;
plot([erpParameters.timeWindow(1) erpParameters.timeWindow(end)],[0 0],'k--')
plot([0 0],[yl(1) yl(2)],'k--')
if groupLmes{bestGroupModel}.significantGlmeFit %has a significant cluster
    %then add shading to indicate significant times
    [~,sigPeriodStart,sigPeriodEnd] = findgaps(find(groupLmes{bestGroupModel}.significantGlmeTimes > 0)); %find start/end of significant clusters
    sigPeriodStart = sigPeriodStart-erpParameters.tWin2; %correct for image onset time
    sigPeriodEnd = sigPeriodEnd-erpParameters.tWin2;%correct for image onset time
    for sPS = 1:length(sigPeriodStart)
        patch([sigPeriodStart(sPS) sigPeriodEnd(sPS) sigPeriodEnd(sPS) sigPeriodStart(sPS)],[yl(1) yl(1) yl(2) yl(2)],'k','FaceAlpha',0.2)
    end
end
hold off
xlim([erpParameters.timeWindow(1) erpParameters.timeWindow(end)])
xlabel('Time from Image Onset (ms)')
ylabel('Beta Weight (a.u.)')
legend(groupLmes{bestGroupModel}.fixedEffectNames,'location','northwest')
title(['Beta Weight Time Course for ' glmParams.groupModelNames{bestGroupModel}])

subplot(2,3,3)
b = bar(1:groupLmes{bestGroupModel}.numFixedEffects,[groupLmes{bestGroupModel}.sigGlmes{1}{clustersToUse(bestGroupModel)}.Coefficients.Estimate(2:end)],'FaceAlpha',0.5); %ignore intercept
b(1).FaceColor = [.2 .6 .5];
xticks(1:groupLmes{bestGroupModel}.numFixedEffects)
xticklabels(groupLmes{bestGroupModel}.fixedEffectNames)
hold on
for nFE = 1:groupLmes{bestGroupModel}.numFixedEffects
    if groupLmes{bestGroupModel}.sigGlmes{1}{clustersToUse(bestGroupModel)}.Coefficients.Estimate(1+nFE) < 0
        offset = -0.05;
    else
        offset = 0.05;
    end
    text(nFE,groupLmes{bestGroupModel}.sigGlmes{1}{clustersToUse(bestGroupModel)}.Coefficients.Estimate(1+nFE)+offset,...
        ['\beta = ' num2str(groupLmes{bestGroupModel}.sigGlmes{1}{clustersToUse(bestGroupModel)}.Coefficients.Estimate(1+nFE),3) ...
        ' a.u., p = ' num2str(groupLmes{bestGroupModel}.sigGlmes{1}{clustersToUse(bestGroupModel)}.Coefficients.pValue(1+nFE),3)])
end
hold off
ylim([-0.2 1])
ylabel('Beta Weight (a.u.)')
box off
title(['Cluster-Level Beta Weights for ' glmParams.groupModelNames{bestGroupModel}])

subplot(2,3,4)
bar(1:length(glmParams.groupModels),groupGlmRSquared)
hold on
marks = groupGlmRSquared.*isMoreComplicatedModelBetter;
marks(marks == 0) = NaN; %not using zero baseline
plot(1:length(glmParams.groupModels),marks,'k*') %mark if more complicated model is better
plot(bestGroupModel,groupGlmRSquared(bestGroupModel),'r*') %mark best model
hold off
xticks(1:length(glmParams.groupModels))
xticklabels(glmParams.groupModelNames)
legend({'Data','Better','best'},'location','northwest','FontSize',8)
ylabel('Adjusted R^2')
title('Gomparison of Group Modles: R^2')

subplot(2,3,5)
bar(1:length(glmParams.groupModels),groupGlmAIC)
hold on
marks = groupGlmAIC.*isMoreComplicatedModelBetter;
marks(marks == 0) = NaN; %not using zero baseline
plot(1:length(glmParams.groupModels),marks,'k*') %mark if more complicated model is better
plot(bestGroupModel,groupGlmAIC(bestGroupModel),'r*') %mark best model
hold off
xticks(1:length(glmParams.groupModels))
xticklabels(glmParams.groupModelNames)
ylabel('AIC')
legend({'Data','Better','best'},'location','northeast','FontSize',8)
ylim([3e4 4e4])
title('Gomparison of Group Modles: AIC')

subplot(2,3,6)
plot(erpParameters.timeWindow,groupLmes{bestGroupModel}.randomEffectTimePoints{1})
xlabel('Time from Image Onset (ms)')
ylabel('Random Effects (std)')
title(['Random Effects Time Course for ' glmParams.groupModelNames{bestGroupModel}])
legend(erase(groupLmes{bestGroupModel}.randomEffectNames,': (Intercept) vs (Intercept)'),'location','northwest','FontSize',8)

sgtitle('Grouped-Face-Selective Broadband Power (BB) Analysis')


%---Plot the Alternative Group-Level Face BB Data---%
figure('units','normalized','outerPosition',[0 0 1 1]);%full screen

subplot(2,2,1)
plot(erpParameters.timeWindow,altGroupLmes{bestAltGroupModel}.betaWeightsTimePoints{1})
hold on
yl = ylim;
plot([erpParameters.timeWindow(1) erpParameters.timeWindow(end)],[0 0],'k--')
plot([0 0],[yl(1) yl(2)],'k--')
if altGroupLmes{bestAltGroupModel}.significantGlmeFit %has a significant cluster
    %then add shading to indicate significant times
    [~,sigPeriodStart,sigPeriodEnd] = findgaps(find(altGroupLmes{bestAltGroupModel}.significantGlmeTimes > 0)); %find start/end of significant clusters
    sigPeriodStart = sigPeriodStart-erpParameters.tWin2; %correct for image onset time
    sigPeriodEnd = sigPeriodEnd-erpParameters.tWin2;%correct for image onset time
    for sPS = 1:length(sigPeriodStart)
        patch([sigPeriodStart(sPS) sigPeriodEnd(sPS) sigPeriodEnd(sPS) sigPeriodStart(sPS)],[yl(1) yl(1) yl(2) yl(2)],'k','FaceAlpha',0.2)
    end
end
hold off
xlim([erpParameters.timeWindow(1) erpParameters.timeWindow(end)])
xlabel('Time from Image Onset (ms)')
ylabel('Beta Weight (a.u.)')
legend(altGroupLmes{bestAltGroupModel}.fixedEffectNames,'location','northwest')
title(['Beta Weight Time Course for ' glmParams.altGlmeNames{bestAltGroupModel}])

subplot(2,2,2)
b = bar(1:altGroupLmes{bestAltGroupModel}.numFixedEffects,[altGroupLmes{bestAltGroupModel}.sigGlmes{1}{altClustersToUse(bestGroupModel)}.Coefficients.Estimate(2:end)],'FaceAlpha',0.5); %ignore intercept
b(1).FaceColor = [.2 .6 .5];
xticks(1:altGroupLmes{bestAltGroupModel}.numFixedEffects)
xticklabels(altGroupLmes{bestAltGroupModel}.fixedEffectNames)
hold on
for nFE = 1:altGroupLmes{bestAltGroupModel}.numFixedEffects
    if altGroupLmes{bestAltGroupModel}.sigGlmes{1}{altClustersToUse(bestGroupModel)}.Coefficients.Estimate(1+nFE) < 0
        offset = -0.05;
    else
        offset = 0.05;
    end
    text(nFE,altGroupLmes{bestAltGroupModel}.sigGlmes{1}{altClustersToUse(bestGroupModel)}.Coefficients.Estimate(1+nFE)+offset,...
        ['\beta = ' num2str(altGroupLmes{bestAltGroupModel}.sigGlmes{1}{altClustersToUse(bestGroupModel)}.Coefficients.Estimate(1+nFE),3) ...
        ' a.u., p = ' num2str(altGroupLmes{bestAltGroupModel}.sigGlmes{1}{altClustersToUse(bestGroupModel)}.Coefficients.pValue(1+nFE),3)])
end
hold off
ylim([-0.2 1])
ylabel('Beta Weight (a.u.)')
box off
title(['Cluster-Level Beta Weights for ' glmParams.altGlmeNames{bestAltGroupModel}])


subplot(2,2,3)
bar(1:length(glmParams.altGlmeFormulas),altGroupGlmRSquared)
hold on
marks = altGroupGlmRSquared.*altIsMoreComplicatedModelBetter;
marks(marks == 0) = NaN; %not using zero baseline
plot(1:length(glmParams.altGlmeFormulas),marks,'k*') %mark if more complicated model is better
plot(bestAltGroupModel,altGroupGlmRSquared(bestAltGroupModel),'r*') %mark best model
hold off
xticks(1:length(glmParams.altGlmeNames))
xticklabels(glmParams.altGlmeNames)
legend({'Data','Better','best'},'location','northwest','FontSize',8)
ylabel('Adjusted R^2')
ylim([0 0.35])
title('Gomparison of Group Modles: R^2')

subplot(2,2,4)
plot(erpParameters.timeWindow,altGroupLmes{bestAltGroupModel}.randomEffectTimePoints{1})
xlabel('Time from Image Onset (ms)')
ylabel('Random Effects (std)')
title(['Random Effects Time Course for ' glmParams.altGlmeNames{bestAltGroupModel}])
legend(erase(altGroupLmes{bestAltGroupModel}.randomEffectNames,': (Intercept) vs (Intercept)'),'location','northwest','FontSize',8)

sgtitle('Alternative Models: Grouped-Face-Selective Broadband Power (BB)')

%%
toc