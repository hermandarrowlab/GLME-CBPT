function glmeClusterStats2D = clusterBasedGlmeERD2D(glmeFormula,glmeDataType,...
    clusterDefType,minClusterSize,timeWindow,allEventData,allEventTables,...
    numShuffs,alphaValue,useParfor,timeWindowOfInterest,clusterDefinitionType)
%function runs cluster based statistics with glmes to determine if 2D data
%(e.g. spectrograms/time-frequency data) with multiple conditioned based responses
%show signficant time-locked locked responses. Different glmes are used based on the
%type of data provided e.g. for spectograms a LME is used based on the assumption of
%normality in noramlized power. Please also check your data
%
%This code is based on the 1D code version for time series data (e.g. ERPs) 
%called clusterBasedGlmeERD!
%
%Important note: if you are only interested in one type of event (e.g.
%stimulus onset) then you should only include the event related data for
%that/for these event(s) in allEventData & allEventTables. For example for stimulus
%onset in N-back you may want to use 51 & 52 but you don't want 17 for
%repsonse times. If you include those extra events as inputs to this
%function (e.g. 17) the GLME will inappropriately run the GLME with 17's data!
%
% Processing steps (done channel by channel):
%   1) Concatenate all the channel's data across all events.
%   2) Fit GLME to each time/freq pixel of interest.
%   3) Find pixels that have significant beta weights for any fixed effect (any p < alpha/numFixEffects).
%   4) Find clusters of significant pixels.
%   5) For each cluster, calculate cluster-level statistic (sum t-stat squared)
%   6) Find the largest cluster-level statistic and calculate the
%   permuatated cluster-level stastic for that largest cluster.
%   7) Compare all clusters to the permuatated cluster-level statistic from
%   the largest cluster. Any clusters whos cluster-level statistic is
%   larger than the permutated values is considered significant.
%
% Inputs:
%   1) glmeFormula: GLME formula in Wilkinson Notation (see https://www.mathworks.com/help/stats/wilkinson-notation.html)
%       *) note the default depedent variable name is 'erd' for event-related data
%   2) glmeDataType: data type being inputs as this effects the type of glme to use.
%       a) 'spectral': uses LME (normal distrbuted glme)
%   3) clusterDefType: cluster definition type
%       a) 'area': cluster based on area, assuming "square" pixels
%       b) 'mag': cluster based on sum of the sum-squared tStat of the effects
%       c) 'areaMag':cluster based on area x mag
%   4) minClusterSize: minimum cluster size 
%       *) can pre-calculate minimum cluster size using findMinimum2DClusterThreshold.m
%       a) for 'area' this is in time samples * freq samples. Since sampling rate
%       varies and freq can be in log-freq units, area doesn't necessarily
%       correspond to ms x Hz! 
%       b) for 'mag' this is in sum-squared-t-stat
%   5) timeWindow: full time window corresponding to timestamps in allEventData (in milliseconds)
%   6) allEventData: cell array of all event-related data organzied by
%   event (rows) and channel (columns) within each cell organized by freq
%   (row), time (columns), and trial (z-col); produced when running ERP-type scripts.
%   7) allEventTables: cell array corresponding to allEventData; this is
%   used for fixed and random effects in the GLME.
%   8) numShuffs: number of permutation shuffles (default is 1000)
%   9) alphaValue: significance value (default is 0.05)
%   10) useParfor: true/false flag to use parfor. If value is > 1 then this
%   is the number of parpools to use. If value is 0 then will use
%   for-loops. Parfor is only applied to permutation testing!
%   11) timeWindowOfInterest: vector of time points of interest to run analysis over.
%   Useful if time window in allEventData is much larger than time window
%   of enterest (e.g. when created for TF analysis and wanted buffers for
%   edge effects).
%   12) clusterDefinitionType: determines how clusters are defined
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
% Outputs (stored in glmeClusterStats2D):
%   1) significantGlmeTimes: binary matrix of times at which GLMEs were
%   significant (1) or not (0) organized by channel (row) and time (column).
%   2) sigGlmes: binary vector on whether a channel had a significant GLME fit (1) or not (0)
%   3) significantGlmeFit: cell array organized by channel with significant
%   GLMEs stored within for each significant time window.
%
%Example Variable Usage For N-back Stimulus Onset ERPs:
% glmeFormula = 'erd ~ neuralEventValue + correctReponse + nBackLevel + (1|sessionNumber)'; %neuralEventValue - > target (51)/nontarget(52)
% glmeDataType = 'spectral';
% minClusterSize = 900;%time/freq samples
% timeWindow = erpParameters.timeWindow;
% timeWindowOfInterest = [timeWindow(1) timeWindow(end)]; %full time window
% useParfor = 1;
% alphaValue = 0.05;
% numShuffs = 1000;
%
%Written by Seth Konig 8/8/2022 based on clusterBasedGlmeERD.m and Blair Vail's pacman_ClusterGLME.m


%---Processing Parameters---%
bwConn = 4;%black-white label connectivity. 4 makes edges have to touch 
%and defualt was 8 which corners only which is weird and not realistic.


%---Parse Inputs---%
%check number of inputs
if nargin < 7
    error('There are 6 parameters inputs that must be speicfied!')
end
if nargin < 12
    clusterDefinitionType = 'allTogether';
end
if nargin < 11
    timeWindowOfInterest = [timeWindow(1) timeWindow(end)];
end
if nargin < 10
    useParfor = 1;
end
if nargin < 9
    alphaValue = 0.05;
end
if nargin < 8
    numShuffs = 1000;
end

%time window parameters
samplingRate = 1e3/mean(diff(timeWindow)); %in Hz
timeWindowOfInterestStart =  find(timeWindow == timeWindowOfInterest(1));
timeWindowOfInterestEnd =  find(timeWindow == timeWindowOfInterest(2));


%get settings for from glmeDataType (distribution 1 for time by time,
%distrbution 2 is for cluster level)
switch lower(glmeDataType)
    case {'spectral'} %essentially LME
        distribution = 'Normal';
        link = 'identity'; %defualt for normal
        estimateDispersion = false; %unnecessary
        distribution2 = 'Normal';
        link2 = 'identity'; %defualt for normal
        estimateDispersion2 = false; %unnecessary
        %     case {'burstrate'} %logistical GLME
        %         distribution = 'Binomial';
        %         link = 'logit'; %default for binomial distrbution
        %         estimateDispersion = true; %may be necessary but likley around 1
        %         distribution2 = 'Poisson';
        %         link2 = 'log'; %default for Poisson distrbution
        %         estimateDispersion2 = true; %may be necessary
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
idx = find(~cellfun(@isempty,allEventData),1);
numFreqencies = size(allEventData{idx},1);
significantGlmeTimes = zeros(numFreqencies,numTimePoints,numChannels); %0s & 1s
significantGlmeTimesL = zeros(numFreqencies,numTimePoints,numChannels);%cluster numbers
sigIndividualTimePoints = cell(1,numChannels);%sig times before clustering for each time point for each channel
tStatIndividualTimePoints = cell(1,numChannels);%tstat for individual times points
betaWeightsTimePoints = cell(1,numChannels); %beta weigths for each time point for each channel
glmeFixedEffectNames = cell(1,numChannels); %beta weight names
randomEffectTimePoints = cell(1,numChannels); %random effect sizes for each time point for each channel
interceptTimePoints = cell(1,numChannels); %intercepts for each time point for each channel (only valid with centered predictors)
glmeRandomEffectNames = cell(1,numChannels);
glmeRandomEffectTypes = cell(1,numChannels); %random types for each time point for each channel
sigGlmes = cell(1,numChannels);
significantGlmeFit = NaN(1,numChannels);
for chan = 1:numChannels
    disp(['Running GLME analysis on channel# ' num2str(chan)])
    
    %---Grab This Channels Data---%
    %concatenate all the data within a channel
    %thisChannelData is a 3D matrix of freq x time x trial
    thisChannelTable = [];
    thisChannelData = [];
    if numEvents == 1 %can have memory issue with large grouped data
        thisChannelTable = allEventTables{1,chan};
        thisChannelData = allEventData{1,chan};
    else
        for event = 1:numEvents
            if ~isempty(allEventTables{event,chan})
                thisChannelTable = [thisChannelTable; allEventTables{event,chan}];
                thisChannelData = cat(3,thisChannelData,allEventData{event,chan}(:,timeWindowOfInterestStart:timeWindowOfInterestEnd,:));
            end
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
    nF = 1;
    thisDataTable = [thisChannelTable array2table(squeeze(thisChannelData(nF,nTP,:)),'VariableNames',{'erd'})];
    if strcmpi(distribution,'normal')
        glme = fitlme(thisDataTable,glmeFormula);
    else
        glme = fitglme(thisDataTable,glmeFormula,'Distribution',distribution,'Link',link,'DispersionFlag',estimateDispersion);
    end
    
    %get number of effects and preallocate space, need some of these
    %parameters too here if we use a parfor since they won't be available
    %outside of the parfor loop
    numFixedEffects = glme.NumCoefficients-1;%ignore intercept
    allGlmePvalues = NaN(numFreqencies,numTimePoints,numFixedEffects);
    thisChanTstatTimePoints = NaN(numFreqencies,numTimePoints,numFixedEffects);
    thisChanBetaWeightsTimePoints = NaN(numFreqencies,numTimePoints,numFixedEffects);
    [thisRandomNames,thisRandomEffectType,~] = getRandomEffectStats(glme);
    [~,~,randomStats] = covarianceParameters(glme);
   
    
    if size(randomStats,1) == 1%error term only
        numRandomEffects = 0;
    else
        numRandomEffects = size(randomStats,1)-1;
    end
    thisRandomEffects = NaN(numFreqencies,numTimePoints,numRandomEffects);
    thisIntercepts = NaN(numFreqencies,numTimePoints);

    %run for rest of time points
    if usingParfor
        for nF = 1:numFreqencies
            parfor nTP = 1:numTimePoints
                try
                    thisDataTable2 = thisDataTable; %save new variable for parfor
                    thisDataTable2.erd = squeeze(thisChannelData(nF,nTP,:)); %update ERD
                    if strcmpi(distribution,'normal')
                        glme = fitlme(thisDataTable2,glmeFormula);
                    else
                        glme = fitglme(thisDataTable2,glmeFormula,'Distribution',distribution,'Link',link,'DispersionFlag',estimateDispersion);
                    end
                    allGlmePvalues(nF,nTP,:) =  glme.Coefficients.pValue(2:end)'; %ignore intercept
                    thisChanTstatTimePoints(nF,nTP,:) = glme.Coefficients.tStat(2:end)'; %ignore intercept
                    thisChanBetaWeightsTimePoints(nF,nTP,:) = glme.Coefficients.Estimate(2:end)'; %ignore intercept
                    [~,~,thisRandomEffects(nF,nTP,:)] = getRandomEffectStats(glme);
                    thisIntercepts(nF,nTP) = glme.Coefficients.Estimate(1);%the intercept
                catch ME
                    if ~strcmpi(ME.identifier,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed')
                        error('encountered an uknown error type!')
                    else %probably an error caused by too few non-zero values (e.g. low burst rates)
                        disp(['Channel #' num2str(chan) ', unable to run GLME for time point ' num2str(nTP) ' so skipping and saying not significant!'])
                        allGlmePvalues(nF,nTP,:) = 1; %not significant;
                    end
                end
            end
        end
    else
        for nF = 1:numFreqencies
            for nTP = 2:numTimePoints
                try
                    thisDataTable2 = thisDataTable; %save new variable for parfor
                    thisDataTable2.erd = squeeze(thisChannelData(nF,nTP,:)); %update ERD
                    if strcmpi(distribution,'normal')
                        glme = fitlme(thisDataTable2,glmeFormula);
                    else
                        glme = fitglme(thisDataTable2,glmeFormula,'Distribution',distribution,'Link',link,'DispersionFlag',estimateDispersion);
                    end
                    allGlmePvalues(nF,nTP,:) =  glme.Coefficients.pValue(2:end)'; %ignore intercept
                    thisChanTstatTimePoints(nF,nTP,:) = glme.Coefficients.tStat(2:end)'; %ignore intercept
                    thisChanBetaWeightsTimePoints(nF,nTP,:) = glme.Coefficients.Estimate(2:end)'; %ignore intercept
                    [~,~,thisRandomEffects(nF,nTP,:)] = getRandomEffectStats(glme);
                    thisIntercepts(nF,nTP) = glme.Coefficients.Estimate(1);%the intercept                    
                catch ME
                    if ~strcmpi(ME.identifier,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed')
                        error('encountered an uknown error type!')
                    else %probably an error caused by too few non-zero values (e.g. low burst rates)
                        disp(['Channel #' num2str(chan) ', unable to run GLME for time point ' num2str(nTP) ' so skipping and saying not significant!'])
                        allGlmePvalues(nF,nTP,:) = 1; %not significant;
                    end
                end
            end
        end
    end
    
    %store beta weights (must do outside parfor)
    tStatIndividualTimePoints{chan} = thisChanTstatTimePoints;
    betaWeightsTimePoints{chan} = thisChanBetaWeightsTimePoints;
    sigIndividualTimePoints{chan} = allGlmePvalues;
    glmeFixedEffectNames{chan} = glme.CoefficientNames(2:end);
    randomEffectTimePoints{chan} = thisRandomEffects;
    interceptTimePoints{chan} = thisIntercepts;
    glmeRandomEffectNames{chan} = thisRandomNames;
    glmeRandomEffectTypes{chan} = thisRandomEffectType;
        
    
    %----Define Clusters---%
    switch clusterDefinitionType
        case {'uniqueCombos','uniquecombos'}
            %get significant time-frequency points and their combinations (bin2dec)
            sigPValues = allGlmePvalues < alphaValue/numFixedEffects;
            sigCombos = NaN(size(sigPValues,1),size(sigPValues,2));
            for sPVr = 1:size(sigPValues,1)
                for sPVc = 1:size(sigPValues,2)
                    sigCombos(sPVr,sPVc) = bin2dec(strrep(num2str(squeeze(sigPValues(sPVr,sPVc,:))'),' ',''));
                end
            end
            uniquePValueCombos = unique(sigCombos(:));
            uniquePValueCombos(uniquePValueCombos == 0) = [];
            
            %define clusters for each combo
            clustersLs = NaN(size(sigCombos,1),size(sigCombos,2),length(uniquePValueCombos));
            allNumClusts = NaN(1,length(uniquePValueCombos));
            for uPVC = 1:length(uniquePValueCombos)
                %find clusters
                theseClusters = sigCombos == uniquePValueCombos(uPVC);
                L = bwlabel(theseClusters,bwConn);
                numClusts = max(L(:));
                clusterArea = NaN(1,numClusts);
                for nC = 1:numClusts
                    clusterArea(nC) = sum(L(:) == nC);
                end
                
                %remove clusters that are too small based on area
                switch clusterDefType
                    case 'area'
                        tooSmall = find(clusterArea < minClusterSize);
                        for tS = 1:length(tooSmall)
                            L(L == tooSmall(tS)) = 0;
                        end
                        
                        %recompute L & clsuter areas
                        L(L ~= 0) = 1;
                        L = bwlabel(L,bwConn);
                        numClusts = max(L(:));
                        clusterArea = NaN(1,numClusts);
                        for nC = 1:numClusts
                            clusterArea(nC) = sum(L(:) == nC);
                        end
                    otherwise
                        error('need code here!')
                end
                clustersLs(:,:,uPVC) = L;
                allNumClusts(uPVC) = numClusts;
            end
            
            %combine clusters across fixed effects
            L = zeros(size(sigCombos));
            newClusterNumber = 1;
            for uPVC = 1:length(uniquePValueCombos)
                for aNC = 1:allNumClusts(uPVC)
                    theseInd = clustersLs(:,:,uPVC) == aNC;
                    L(theseInd == 1) = newClusterNumber;
                    newClusterNumber = newClusterNumber+1;
                end
            end
            numClusts = newClusterNumber-1;
            
            
        case {'allTogether','alltogether'}
            %find time points that are signficiant for any
            combinedPvalues = any(allGlmePvalues < alphaValue/numFixedEffects,3);
            
            %find clusters that are big enough
            L = bwlabel(combinedPvalues,bwConn);
            numClusts = max(L(:));
            clusterArea = NaN(1,numClusts);
            for nC = 1:numClusts
                clusterArea(nC) = sum(L(:) == nC);
            end
            
            %remove clusters that are too small based on area
            switch clusterDefType
                case 'area'
                    tooSmall = find(clusterArea < minClusterSize);
                    for tS = 1:length(tooSmall)
                        L(L == tooSmall(tS)) = 0;
                    end
                    
                    %recompute L & clsuter areas
                    L(L ~= 0) = 1;
                    L = bwlabel(L,bwConn);
                    numClusts = max(L(:));
                    clusterArea = NaN(1,numClusts);
                    for nC = 1:numClusts
                        clusterArea(nC) = sum(L(:) == nC);
                    end
                otherwise
                    error('need code here!')
            end
            
        otherwise
            error('cluster definition type not recognized')
    end
    
    
    %---Do Cluster Level Stats (shuffling)---%
    %get cluster-level stats using sum square t-stat for cluster-level statistic
    if numClusts > 0
        thisClusterData = cell(1,numClusts);
        clusterGlmes = cell(1,numClusts);
        sumSquaredTstats = NaN(1,numClusts);
        for nC = 1:numClusts
            %get average data for this cluster for cluster-level measure
            %varies by method
            indexMatrix = L == nC;%2D matrix of this clusters indices
            switch lower(glmeDataType)
                case {'spectral'} %essentially LME, take the mean
                    thisClusterData{nC} = NaN(size(thisChannelData,3),1);
                    for trial = 1:size(thisChannelData,3)
                        thisTrial = thisChannelData(:,:,trial);
                        thisClusterData{nC}(trial) = mean(thisTrial(indexMatrix));
                    end
                %case {'burstrate'} %logistical GLME
                    %thisClusterData = sum(thisChannelData(:,clusterStart(cS):clusterEnd(cS)),2);
                otherwise
                    error('glmeDataType type not recognized/supported!')
            end
            
            %fit glme to cluster and get t-stat
            thisDataTableCluster = [thisChannelTable array2table(thisClusterData{nC},'VariableNames',{'erd'})];
            try
                if strcmpi(distribution,'normal')
                    clusterGlmes{nC} = fitlme(thisDataTableCluster,glmeFormula);
                else
                    clusterGlmes{nC} = fitglme(thisDataTableCluster,glmeFormula,'Distribution',distribution2,'Link',link2,'DispersionFlag',estimateDispersion2);
                end
                sumSquaredTstats(nC) = sum((clusterGlmes{nC}.Coefficients.tStat(2:end)).^2); %ignore intercept
            catch
                continue
            end
        end
        
        %find biggest cluster stat and get cluster-level data (to be used
        %below to fill table below)
        biggestClust = find(sumSquaredTstats == max(sumSquaredTstats));
        biggestClusterData = thisClusterData{biggestClust};

        %check if clusters are significant by getting permutated/shuffled stats for this cluster
        numTrials = length(biggestClusterData);
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
        sigClusters = NaN(1,numClusts);
        for nC = 1:numClusts
            if sum(sumSquaredTstats(nC) > shuffLmeTStat) > (1-alphaValue)*numShuffs
                sigClusters(nC) = 1;
            else
                sigClusters(nC) = 0;
            end
        end
        
        %create new matrix of sig times
        sigClusterTimes = zeros(size(L));
        for nC = 1:numClusts
            if sigClusters(nC) == 1
                sigClusterTimes(L == nC) = nC;
            else
                L(L == nC) = 0;  %remove the clusters from L
            end
        end
        
        %---Store Data---%
        %store sig times and glmes
        if any(sigClusters == 1)
            significantGlmeFit(chan) = 1;
            significantGlmeTimes(:,:,chan) = L ~= 0;%general indices
            significantGlmeTimesL(:,:,chan) = L;%cluster numbers
            for nC = 1:numClusts
                if sigClusters(nC) == 1
                    sigGlmes{chan} = [sigGlmes{chan} clusterGlmes(nC)];
                else
                    sigGlmes{chan} = [sigGlmes{chan} {}];%empty so that we can keep L's indexing!
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
glmeClusterStats2D = [];
glmeClusterStats2D.numFixedEffects = numFixedEffects;
glmeClusterStats2D.numRandomEffects = numRandomEffects;
glmeClusterStats2D.fixedEffectNames = glmeFixedEffectNames{1};
glmeClusterStats2D.randomEffectNames = glmeRandomEffectNames{1};
glmeClusterStats2D.significantGlmeTimes = significantGlmeTimes;
glmeClusterStats2D.significantGlmeTimesL = significantGlmeTimesL;
glmeClusterStats2D.sigGlmes = sigGlmes;
glmeClusterStats2D.significantGlmeFit = significantGlmeFit;
glmeClusterStats2D.tStatIndividualTimePoints = tStatIndividualTimePoints;
glmeClusterStats2D.betaWeightsTimePoints = betaWeightsTimePoints;
glmeClusterStats2D.interceptTimePoints = interceptTimePoints;
glmeClusterStats2D.sigIndividualTimePoints = sigIndividualTimePoints;
glmeClusterStats2D.randomEffectTimePoints = randomEffectTimePoints;

end