function glmeClusterStats = clusterBasedGlmeERD(glmeFormula,glmeDataType,minClusterSize,...
    timeWindow,allEventData,allEventTables,numShuffs,alphaValue,useParfor,...
    timeWindowOfInterest,clusterDefinitionType)
%function runs cluster based statistics with glmes to determine if ERP, burst rate, or
%other event related data (ERD) with multiple conditioned based responses show
%significant time-locked locked responses. Different glmes are used based on the
%type of data provided e.g. for ERPs a LME is used based on the assumption of
%normality in ERP voltages while for burst rate a logistic regression is used based on it's
%binary distribution of bursting or not.
%
% Code does not handle 2D data (e.g. TF analysis). SEE clusterBasedGlmeERD2D.m
%
%Important note: if you are only interested in one type of event (e.g.
%stimulus onset) then you should only include the event related data for
%that/for these event(s) in allEventData & allEventTables. For example for stimulus
%onset in N-back you may want to use 51 & 52 but you don't want 17 for
%response times. If you include those extra events as inputs to this
%function (e.g. 17) the GLME will inappropriately run the GLME with 17's data!
%
% Processing steps (done channel by channel):
%   1) Concatenate all the channel's data across all events.
%   2) Fit GLME to each time point of interest.
%   3) Find time points that have significant beta weights for any fixed effect (any p < alpha/numFixEffects).
%   4) Find clusters of significant time points >= minClusterSize
%   5) For each cluster, calculate cluster-level statistic (sum all t-stat squared)
%   6) Find the largest cluster-level statistic and calculate the
%   permutated cluster-level statistic for that largest cluster.
%   7) Compare all clusters to the permutated cluster-level statistic from
%   the largest cluster. Any clusters whose cluster-level statistic is
%   larger than the permutated values is considered significant.
%
%
% Inputs:
%   1) glmeFormula: GLME formula in Wilkinson Notation (see https://www.mathworks.com/help/stats/wilkinson-notation.html)
%       *) note the default dependent variable name is 'erd' for all
%       event-related data types!
%   2) glmeDataType: data type being inputs as this effects the type of glme to use.
%       a) 'ERP': uses LME (normal distrbuted glme)
%       b) 'burstRate': uses logistical GLME (binary data)
%   3) minClusterSize: minimum cluster size in milliseconds (e.g 50 ms)
%   4) timeWindow: full time window corresponding to timestamps in allEventData (in milliseconds)
%   5) allEventData: cell array of all event-related data organized by
%   event (rows) and channel (columns) within each cell organized by trial
%   (row) and time (columns); produced when running ERP-type scripts.
%   6) allEventTables: cell array corresponding to allEventData; this is
%   used for fixed and random effects in the GLME.
%   7) numShuffs: number of permutation shuffles (default is 1000)
%   8) alphaValue: significance value (default is 0.05)
%   9) useParfor: true/false flag to use parfor. If value is > 1 then this
%   is the number of parpools to use. If value is 0 then will use
%   for-loops. Parfor is only applied to permutation testing!
%   10) timeWindowOfInterest: vector of time points of interest to run analysis over.
%   Useful if time window in allEventData is much larger than time window
%   of interest (e.g. when created for TF analysis and wanted buffers for
%   edge effects).
%   11) clusterDefinitionType: determines how clusters are defined
%       a) 'allTogether' (default): combines all fixed effects together and 
%       identifies clusters by p-values of ANY fixed effect. This method is
%       will likely identify fewer clusters but may wash out smaller
%       effects. This method is good for simpler signals and tasks.
%       b) 'uniqueCombos': defines cluster for unique combinations of significant 
%       effects. For example, if there are 2 fixed effects, clusters are
%       defined as no sig effects, sig for FE 1, sig for F2, or sig for FE1
%       & F2. This method is more likely to find a greater number of clusters 
%       and clusters are more likely to be smaller. This method is good for
%       more complicate signals and tasks where partial overlap of significant
%       time courses is likely. 
%
% Outputs (stored in glmeClusterStats):
%   1) numFixedEffects: number of fixed effects
%   2) numRandomEffects: number of random effects
%   3) fixedEffectNames: fixed effect names
%   4) randomEffectNames: random effect names
%   5) significantGlmeTimes: 0's and #'s for significant cluster-corrected
%   times organized by channel (row) and time (column)
%       a) The number(#) corresponds to the cluster number.
%   6) sigGlmes: significant cluster-level GLMEs organized by {channel} then
%   {cluster}
%   7) significantGlmeFit: vector of 0's and 1's representing whether a
%   channel has a significant cluster anywhere.
%   8) betaWeightsTimePoints: beta weights for each time point and fixed
%   effect organized by {channel}, time (row), and fixed effect (col).
%   9) sigIndividualTimePoints: significant of each time point for each
%   fixed effected organized by {channel}, time (row), and fixed effect (col).
%       *)These values are uncorrected p-values.
%   10) randomEffectTimePoints: random effect values for each time point
%   and random effect organized by {channel}, time (row), and random effect (col).
%   11) seIndividualTimePoints: standard error of the beta weights at
%   individual time points following the same organization as betaWeightsTimePoints.
%   12) tStatIndividualTimePoints: t-statistic of the beta weights at
%   individual time points following the same organization as betaWeightsTimePoints.
%   13) ciLowerIndividualTimePoints: lower bound of confidennce interval of
%   the beta weights at individual time points following the same organization as betaWeightsTimePoints.
%   14) ciUpperIndividualTimePoints: upper bound of confidennce interval of
%   the beta weights at individual time points following the same organization as betaWeightsTimePoints.
%
%
%Example Variable Usage For N-back Stimulus Onset ERPs:
% glmeFormula = 'erd ~ neuralEventValue + correctReponse + nBackLevel + (1|sessionNumber)'; %neuralEventValue - > target (51)/nontarget(52)
% glmeDataType = 'ERP';
% minClusterSize = 80;%ms
% timeWindow = erpParameters.timeWindow;
% timeWindowOfInterest = [timeWindow(1) timeWindow(end)]; %full time window
% useParfor = 1;
% alphaValue = 0.05;
% numShuffs = 1000;
%
%
%Written by Seth Konig 3/14/2022, last updated 12/16/2022


%---Parse Inputs---%
%check number of inputs
if nargin < 6
    error('There are 6 parameters inputs that must be speicfied!')
end
if nargin < 11
    clusterDefinitionType = 'allTogether';
end
if nargin < 10
    timeWindowOfInterest = [timeWindow(1) timeWindow(end)];
end
if nargin < 9
    useParfor = 1;
end
if nargin < 8
    alphaValue = 0.05;
end
if nargin < 7
    numShuffs = 1000;
end

%time window parameters
samplingRate = 1e3/mean(diff(timeWindow)); %in Hz
minClusterSizeSamples = minClusterSize/1e3 * samplingRate; %in samples
timeWindowOfInterestStart =  find(timeWindow == timeWindowOfInterest(1));
timeWindowOfInterestEnd =  find(timeWindow == timeWindowOfInterest(2));

%get settings for from glmeDataType (distribution 1 for time by time,
%distrbution 2 is for cluster level)
switch lower(glmeDataType)
    case {'erp'} %essentially LME
        distribution = 'Normal';
        link = 'identity'; %defualt for normal
        estimateDispersion = false; %unnecessary
        distribution2 = 'Normal';
        link2 = 'identity'; %defualt for normal
        estimateDispersion2 = false; %unnecessary
    case {'burstrate'} %logistical GLME
        distribution = 'Binomial';
        link = 'logit'; %default for binomial distrbution
        estimateDispersion = true; %may be necessary but likely around 1
        distribution2 = 'Poisson';
        link2 = 'log'; %default for Poisson distrbution
        estimateDispersion2 = true; %may be necessary
    otherwise
        error('glmeDataType type not recognized/supported!')
end

%check parpool parameters
if useParfor >= 1
    usingParfor = true;
    if isempty(gcp('nocreate')) %no pool is open
        try
            p = parallel.defaultProfile; %for 2022b and later
        catch
            p = parallel.defaultClusterProfile; %get default pool name for earlier versions
        end
        if useParfor == 1 %use default number of parallel pools
            parpool(p); %creates local parpool
        else
            parpool(p,useParfor) %creates a local parpool with the number of works specified in useParfor
        end
    end
else
    usingParfor = false;
end



%---Run GLME on Each Channels---%
numEvents = size(allEventData,1);
numChannels = size(allEventData,2);
numTimePoints = timeWindowOfInterestEnd-timeWindowOfInterestStart+1;
significantGlmeTimes = zeros(numChannels,numTimePoints);
seIndividualTimePoints = cell(1,numChannels);
tStatIndividualTimePoints = cell(1,numChannels);
ciLowerIndividualTimePoints = cell(1,numChannels);
ciUpperIndividualTimePoints= cell(1,numChannels);
sigIndividualTimePoints = cell(1,numChannels);%sig times before clustering for each time point for each channel
betaWeightsTimePoints = cell(1,numChannels); %beta weigths for each time point for each channel
glmeFixedEffectNames = cell(1,numChannels); %beta weight names
randomEffectTimePoints = cell(1,numChannels); %random effect sizes for each time point for each channel
glmeRandomEffectNames = cell(1,numChannels);
glmeRandomEffectTypes = cell(1,numChannels); %random types for each time point for each channel
sigGlmes = cell(1,numChannels);
significantGlmeFit = NaN(1,numChannels);
for chan = 1:numChannels
    disp(['Running GLME analysis on channel# ' num2str(chan)])
    
    %---Grab This Channels Data---%
    %concatenate all the data within a channel
    thisChannelTable = [];
    thisChannelData = [];
    for event = 1:numEvents
        if ~isempty(allEventTables{event,chan})
            thisChannelTable = [thisChannelTable; allEventTables{event,chan}];
            thisChannelData = [thisChannelData; allEventData{event,chan}(:,timeWindowOfInterestStart:timeWindowOfInterestEnd)];
        end
    end
    
    %check if channel is empty
    if isempty(thisChannelTable)
        disp(['No data for channel# ' num2str(chan)])
        continue
    end
    
    
    %---Fig GLMEs To Each Time Point (no shuffling)---%
    %get time time point by time point glmes
    %run once so we can figure out the size of everything cuz I don't
    %know how to get number of fixed effects and this will screw with parfor
    nTP = 1;
    thisDataTable = [thisChannelTable array2table(thisChannelData(:,nTP),'VariableNames',{'erd'})];
    if strcmpi(distribution,'normal')
        glme = fitlme(thisDataTable,glmeFormula);
    else
        glme = fitglme(thisDataTable,glmeFormula,'Distribution',distribution,'Link',link,'DispersionFlag',estimateDispersion);
    end
    
    %get number of effects and preallocate space
    numFixedEffects = glme.NumCoefficients-1;%ignore intercept
    allGlmePvalues = NaN(numTimePoints,numFixedEffects);
    thisChanBetaWeightsTimePoints = NaN(numTimePoints,numFixedEffects);
    thisChanSE =  NaN(numTimePoints,numFixedEffects);
    thisChantStat =  NaN(numTimePoints,numFixedEffects);
    thisChan95CILower =  NaN(numTimePoints,numFixedEffects);
    thisChan95CIUpper =  NaN(numTimePoints,numFixedEffects);
    [~,~,randomStats] = covarianceParameters(glme);
    if size(randomStats,1) == 1%error term only
        numRandomEffects = 0;
    else
        numRandomEffects = size(randomStats,1)-1;
    end
    thisRandomEffects = NaN(numTimePoints,numRandomEffects);
    
    %store estimates
    allGlmePvalues(nTP,:) =  glme.Coefficients.pValue(2:end)'; %ignore intercept
    thisChanBetaWeightsTimePoints(nTP,:) = glme.Coefficients.Estimate(2:end)'; %ignore intercept
    thisChanSE(nTP,:)= glme.Coefficients.SE(2:end)'; %ignore intercept
    thisChantStat(nTP,:)= glme.Coefficients.tStat(2:end)'; %ignore intercept
    thisChan95CILower(nTP,:)= glme.Coefficients.Lower(2:end)'; %ignore intercept
    thisChan95CIUpper(nTP,:)= glme.Coefficients.Upper(2:end)'; %ignore intercept
    [thisRandomNames,thisRandomEffectType,thisRandomEffects(nTP,:)] = getRandomEffectStats(glme);
    
    
    %run for rest of time points
    if usingParfor
        parfor nTP = 2:numTimePoints
            try
                thisDataTable2 = thisDataTable; %save new variable for parfor
                thisDataTable2.erd = thisChannelData(:,nTP); %update ERD
                if strcmpi(distribution,'normal')
                    glme = fitlme(thisDataTable2,glmeFormula);
                else
                    glme = fitglme(thisDataTable2,glmeFormula,'Distribution',distribution,'Link',link,'DispersionFlag',estimateDispersion);
                end
                allGlmePvalues(nTP,:) =  glme.Coefficients.pValue(2:end)'; %ignore intercept
                thisChanBetaWeightsTimePoints(nTP,:) = glme.Coefficients.Estimate(2:end)'; %ignore intercept
                thisChanSE(nTP,:)= glme.Coefficients.SE(2:end)'; %ignore intercept
                thisChantStat(nTP,:)= glme.Coefficients.tStat(2:end)'; %ignore intercept
                thisChan95CILower(nTP,:)= glme.Coefficients.Lower(2:end)'; %ignore intercept
                thisChan95CIUpper(nTP,:)= glme.Coefficients.Upper(2:end)'; %ignore intercept
                [~,~,thisRandomEffects(nTP,:)] = getRandomEffectStats(glme);
            catch ME
                if ~strcmpi(ME.identifier,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed')
                    error('encountered an uknown error type!')
                else %probably an error caused by too few non-zero values (e.g. low burst rates)
                    disp(['Channel #' num2str(chan) ', unable to run GLME for time point ' num2str(nTP) ' so skipping and saying not significant!'])
                    allGlmePvalues(nTP,:) = 1; %not significant;
                end
            end
        end
    else
        for nTP = 2:numTimePoints
            try
                thisDataTable2 = thisDataTable; %save new variable for parfor
                thisDataTable2.erd = thisChannelData(:,nTP); %update ERD
                if strcmpi(distribution,'normal')
                    glme = fitlme(thisDataTable2,glmeFormula);
                else
                    glme = fitglme(thisDataTable2,glmeFormula,'Distribution',distribution,'Link',link,'DispersionFlag',estimateDispersion);
                end
                allGlmePvalues(nTP,:) =  glme.Coefficients.pValue(2:end)'; %ignore intercept
                thisChanBetaWeightsTimePoints(nTP,:) = glme.Coefficients.Estimate(2:end)'; %ignore intercept
                thisChanSE(nTP,:)= glme.Coefficients.SE(2:end)'; %ignore intercept
                thisChantStat(nTP,:)= glme.Coefficients.tStat(2:end)'; %ignore intercept
                thisChan95CILower(nTP,:)= glme.Coefficients.Lower(2:end)'; %ignore intercept
                thisChan95CIUpper(nTP,:)= glme.Coefficients.Upper(2:end)'; %ignore intercept
                [~,~,thisRandomEffects(nTP,:)] = getRandomEffectStats(glme);
            catch ME
                if ~strcmpi(ME.identifier,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed')
                    error('encountered an uknown error type!')
                else %probably an error caused by too few non-zero values (e.g. low burst rates)
                    disp(['Channel #' num2str(chan) ', unable to run GLME for time point ' num2str(nTP) ' so skipping and saying not significant!'])
                    allGlmePvalues(nTP,:) = 1; %not significant;
                end
            end
        end
    end
    
    %store beta weights (must do outside parfor)
    betaWeightsTimePoints{chan} = thisChanBetaWeightsTimePoints;
    sigIndividualTimePoints{chan} = allGlmePvalues;
    seIndividualTimePoints{chan} = thisChanSE;
    tStatIndividualTimePoints{chan} = thisChantStat;
    ciLowerIndividualTimePoints{chan} =  thisChan95CILower;
    ciUpperIndividualTimePoints{chan} =  thisChan95CIUpper;
    glmeFixedEffectNames{chan} = glme.CoefficientNames(2:end);
    randomEffectTimePoints{chan} = thisRandomEffects;
    glmeRandomEffectNames{chan} = thisRandomNames;
    glmeRandomEffectTypes{chan} = thisRandomEffectType;
    
    
    %----Define Clusters---%
    switch clusterDefinitionType
        case {'uniqueCombos','uniquecombos'}
            %get significant time points and their combinations (bin2dec)
            sigPValues = allGlmePvalues < alphaValue/numFixedEffects;
            sigCombos = NaN(size(sigPValues,1),1);
            for sPV = 1:size(sigPValues,1)
                sigCombos(sPV) = bin2dec(strrep(num2str(sigPValues(sPV,:)),' ',''));
            end
            uniquePValueCombos = unique(sigCombos);
            uniquePValueCombos(uniquePValueCombos == 0) = [];
            
            %find clusters 
            clusterStart = [];
            clusterEnd = [];
            for uPVC = 1:length(uniquePValueCombos)
                theseInd = find(sigCombos == uniquePValueCombos(uPVC));
                [~,theseClusterStarts,theseClusterEnds] = findgaps(theseInd);
                clusterStart = [clusterStart theseClusterStarts];
                clusterEnd = [clusterEnd theseClusterEnds];
            end
            
            %resort
            [~,si] = sort(clusterStart);
            clusterStart = clusterStart(si);
            clusterEnd = clusterEnd(si);
            
            %remove clusters that are too short
             clusterDuration = clusterEnd-clusterStart+1;
            tooSmall = clusterDuration < minClusterSizeSamples;
            clusterDuration(tooSmall) = [];
            clusterStart(tooSmall) = [];
            clusterEnd(tooSmall) = [];
            
            
        case {'allTogether','alltogether'}
            %find time points that are signficiant for any
            if numFixedEffects == 1
                combinedPvalues = allGlmePvalues < alphaValue;
            else
                combinedPvalues = any(allGlmePvalues' < alphaValue/numFixedEffects);
            end
            
            %find clusters that are big enough
            [~,clusterStart,clusterEnd] = findgaps(find(combinedPvalues == 1));
            clusterDuration = clusterEnd-clusterStart+1;
            tooSmall = clusterDuration < minClusterSizeSamples;
            clusterDuration(tooSmall) = [];
            clusterStart(tooSmall) = [];
            clusterEnd(tooSmall) = [];
            
        otherwise
            error('cluster definition type unknown!')
    end
    
    
    %---Do Cluster Level Stats (shuffling)---%
    %get cluster-level stats using sum square t-stat for cluster-level statistic
    if ~isempty(clusterStart)
        clusterGlmes = cell(1,length(clusterStart));
        sumSquaredTstats = NaN(1,length(clusterStart));
        for cS = 1:length(clusterStart)
            %get average data for this cluster for cluster-level measure
            %varies by method
            switch lower(glmeDataType)
                case {'erp'} %essentially LME
                    thisClusterData = mean(thisChannelData(:,clusterStart(cS):clusterEnd(cS)),2);
                case {'burstrate'} %logistical GLME
                    thisClusterData = sum(thisChannelData(:,clusterStart(cS):clusterEnd(cS)),2);
                otherwise
                    error('glmeDataType type not recognized/supported!')
            end
            
            %fit glme to cluster and get t-stat
            thisDataTableCluster = [thisChannelTable array2table(thisClusterData,'VariableNames',{'erd'})];
            try
                if strcmpi(distribution,'normal')
                    clusterGlmes{cS} = fitlme(thisDataTableCluster,glmeFormula);
                else
                    clusterGlmes{cS} = fitglme(thisDataTableCluster,glmeFormula,'Distribution',distribution2,'Link',link2,'DispersionFlag',estimateDispersion2);
                end
                sumSquaredTstats(cS) = sum((clusterGlmes{cS}.Coefficients.tStat(2:end)).^2); %ignore intercept
            catch
                continue
            end
        end
        
        %find biggest cluster stat and get cluster-level data (to be used
        %below to fill table below)
        biggestClust = find(sumSquaredTstats == max(sumSquaredTstats));
        switch lower(glmeDataType)
            case {'erp'} %essentially LME
                biggestClusterData = mean(thisChannelData(:,clusterStart(biggestClust):clusterEnd(biggestClust)),2);
            case {'burstrate'} %logistical GLME
                biggestClusterData = sum(thisChannelData(:,clusterStart(biggestClust):clusterEnd(biggestClust)),2);
            otherwise
                error('glmeDataType type not recognized/supported!')
        end
        
        %check if clusters are significant by getting permutated/shuffled stats for this cluster
        numTrials = size(thisClusterData,1);
        shuffLmeTStat = NaN(1,numShuffs);
        if usingParfor
            parfor shuff = 1:numShuffs
                try
                    randInd = randperm(numTrials);
                    thisDataTableShuff = thisDataTableCluster;
                    thisBiggestClusterData = biggestClusterData; %reduces over head in parfor according to matlab
                    thisDataTableShuff.erd = thisBiggestClusterData(randInd);
                    if strcmpi(distribution,'normal')
                        shuffLME = fitlme(thisDataTableShuff,glmeFormula);
                    else
                        shuffLME = fitglme(thisDataTableShuff,glmeFormula,'Distribution',distribution2,'Link',link2,'DispersionFlag',estimateDispersion2);
                    end
                    shuffLmeTStat(shuff) = sum((shuffLME.Coefficients.tStat(2:end)).^2); %ignore intercept
                catch
                    continue
                end
            end
        else
            for shuff = 1:numShuffs
                randInd = randperm(numTrials);
                thisDataTableShuff = thisDataTableCluster;
                thisBiggestClusterData = biggestClusterData; %reduces over head in parfor according to matlab
                thisDataTableShuff.erd = thisBiggestClusterData(randInd);
                if strcmpi(distribution,'normal')
                    shuffLME = fitlme(thisDataTableShuff,glmeFormula);
                else
                    shuffLME = fitglme(thisDataTableShuff,glmeFormula,'Distribution',distribution2,'Link',link2,'DispersionFlag',estimateDispersion2);
                end
                shuffLmeTStat(shuff) = sum((shuffLME.Coefficients.tStat(2:end)).^2); %ignore intercept
            end
        end
        
        %determine significant clusters
        sigClusters = NaN(1,length(clusterStart));
        for cS = 1:length(clusterStart)
            if sum(sumSquaredTstats(cS) > shuffLmeTStat) > (1-alphaValue)*numShuffs
                sigClusters(cS) = 1;
            else
                sigClusters(cS) = 0;
            end
        end
        
        %---Store Data---%
        %store times and lmes
        if any(sigClusters == 1)
            significantGlmeFit(chan) = 1;
            for cS = 1:length(clusterStart)
                if sigClusters(cS) == 1
                    significantGlmeTimes(chan,clusterStart(cS):clusterEnd(cS)) = cS;
                    sigGlmes{chan} = [sigGlmes{chan} clusterGlmes(cS)];
                end
            end
        else
            significantGlmeFit(chan) = 0;
        end
    else
        significantGlmeFit(chan) = 0;
    end
end



%---Store Output---%
glmeClusterStats = [];
glmeClusterStats.numFixedEffects = numFixedEffects;
glmeClusterStats.numRandomEffects = numRandomEffects;
glmeClusterStats.fixedEffectNames = glmeFixedEffectNames{1};
glmeClusterStats.randomEffectNames = glmeRandomEffectNames{1};
glmeClusterStats.significantGlmeTimes = significantGlmeTimes;
glmeClusterStats.sigGlmes = sigGlmes;
glmeClusterStats.significantGlmeFit = significantGlmeFit;
glmeClusterStats.betaWeightsTimePoints = betaWeightsTimePoints;
glmeClusterStats.sigIndividualTimePoints = sigIndividualTimePoints;
glmeClusterStats.randomEffectTimePoints = randomEffectTimePoints;
glmeClusterStats.seIndividualTimePoints = seIndividualTimePoints;
glmeClusterStats.tStatIndividualTimePoints = tStatIndividualTimePoints;
glmeClusterStats.ciLowerIndividualTimePoints = ciLowerIndividualTimePoints;
glmeClusterStats.ciUpperIndividualTimePoints = ciUpperIndividualTimePoints;

end


function [thisRandomNames,thisRandomEffectType,thisRandomEffects] = getRandomEffectStats(glme)
%own subfunction to keep it clean
%if more than random effect then includes corr

%get random effects
[~,~,randomStats] = covarianceParameters(glme);
if size(randomStats,1) == 1 %then error term only
    thisRandomNames = [];
    thisRandomEffectType = [];
    thisRandomEffects = [];
    return
end
numRandomEffects = size(randomStats,1)-1;

%store in more convient format
thisRandomEffectType = cell(1,numRandomEffects);
thisRandomNames = cell(1,numRandomEffects);
thisRandomEffects = NaN(1,numRandomEffects);
for nRE = 1:numRandomEffects
    thisRandomEffectType{nRE} = randomStats{nRE}.Type{1};
    thisRandomNames{nRE} = [randomStats{nRE}.Group ': ' randomStats{nRE}.Name1{1} ' vs ' randomStats{nRE}.Name2{1}];
    thisRandomEffects(nRE) =  randomStats{nRE}.Estimate(1);
end

end