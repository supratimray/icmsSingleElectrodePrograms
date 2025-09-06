% Display All Channels after grouping them based on electrode distance from
% the stimulation electrode.

% In these protocols
% azimuth is mapped to amplitude of the stimulation
% elevation is mapped to the frequency of stimulation (in the old version,
% it was the number of pulses)
% temporal frequency is mapped to duration (in the old version, the
% frequency of stimulation)

function displayAllChannelsICMSSingleElectrode(subjectName,expDate,protocolName,folderSourceString,stimulationElectrode,badTrialNameStr,useCommonBadTrialsFlag)

if ~exist('folderSourceString','var');   folderSourceString='E:';       end
if ~exist('badTrialNameStr','var');     badTrialNameStr = '_v5';        end
if ~exist('useCommonBadTrialsFlag','var'); useCommonBadTrialsFlag = 1;  end

gridType = 'Microelectrode';

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');

% Get Combinations
[~,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);

% Guess the type of protocol
if ~isscalar(aValsUnique) && isscalar(eValsUnique) && isscalar(tValsUnique)
    protocolType = 1; 
    condVals = aValsUnique;
elseif isscalar(aValsUnique) && ~isscalar(eValsUnique) && isscalar(tValsUnique)
    protocolType = 2;
    condVals = eValsUnique;
else
    error('Parameter combinations not in valid format');
end
numConditions = length(condVals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.25; panelStartHeight = 0.7;
staticPanelWidth = 0.2; staticStartPos = 0.1;
dynamicPanelWidth = 0.2; dynamicStartPos = 0.3;
timingPanelWidth = 0.2; timingStartPos = 0.5;
plotOptionsPanelWidth = 0.2; plotOptionsStartPos = 0.7;
backgroundColor = 'w';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.12; dynamicGap=0.02; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Amplitude
amplitudeString = getStringFromValues(aValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Amplitude (microAmps)','FontSize',fontSizeSmall);

if protocolType == 1
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','text','String','variable','FontSize',fontSizeSmall);
else
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','text','String',amplitudeString,'FontSize',fontSizeSmall);
end

% Frequency
freqString = getStringFromValues(eValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Freq / #Pulses','FontSize',fontSizeSmall);

if protocolType == 2
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','text','String','variable','FontSize',fontSizeSmall);
else
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','text','String',freqString,'FontSize',fontSizeSmall);
end

% Sigma
sigmaString = getStringFromValues(sValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Sigma (Deg)','FontSize',fontSizeSmall);
hSigma = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-3*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',sigmaString,'FontSize',fontSizeSmall);

% Spatial Frequency
spatialFreqString = getStringFromValues(fValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeSmall);
hSpatialFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-4*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeSmall);

% Orientation
orientationString = getStringFromValues(oValsUnique,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeSmall);
hOrientation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeSmall);

% Contrast
if ~isempty(cValsUnique)
    contrastString = getStringFromValues(cValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Contrast (%)','FontSize',fontSizeSmall);
    hContrast = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-6*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',contrastString,'FontSize',fontSizeSmall);
end

% Duration or frequency
if ~isempty(tValsUnique)
    temporalFreqString = getStringFromValues(tValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-7*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Dur / Freq','FontSize',fontSizeSmall);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-7*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','text','String',temporalFreqString,'FontSize',fontSizeSmall);
end

% Reference scheme
referenceChannelString = 'None';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.11; timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

signalRange = [-0.5 1.5];
freqRange = [0 100];
baseline = [-0.7 -0.2];
stimPeriod = [0.75 1.25];

% Signal Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Parameter','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Min','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeMedium);

% Stim Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);

% Freq Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-4*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Freq Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-4*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(freqRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-4*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(freqRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Z Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-9*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Z Range','FontSize',fontSizeSmall);
hZMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','-5','FontSize',fontSizeSmall);
hZMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','10','FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.15;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale Z','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleZ_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale XY','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show electrode array and bad channels
% Get Bad channels from the main impedance file

impedanceFileName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,'impedanceValues.mat');
badImpedanceCutoff = 2500;

if exist(impedanceFileName,'file')
    impedanceValues = getImpedanceValues(impedanceFileName);
    badChannels = [find(impedanceValues>badImpedanceCutoff) find(isnan(impedanceValues))];
else
    disp('Could not find impedance values');
    badChannels=[];
end

%%%%%%%%%%%%%%%%%%%%%%%% Get good electrode lists %%%%%%%%%%%%%%%%%%%%%%%%%
electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
[~,~,electrodeArray] = electrodePositionOnGrid(1,gridType,subjectName);
[electrodeGroupList,groupNameList,goodElectrodes] = getElectrodeGroups(subjectName,gridType,electrodeArray,stimulationElectrode,badChannels);

numElectrodeGroups = length(electrodeGroupList);
colorNamesElectrodeGroups = copper(numElectrodeGroups);
for iG=1:numElectrodeGroups
    hElectrodes = showElectrodeLocations(electrodeGridPos,electrodeGroupList{iG},colorNamesElectrodeGroups(iG,:),[],1,0,gridType,subjectName);
    text(hElectrodes,-0.35,iG/10,groupNameList{iG},'color',colorNamesElectrodeGroups(iG,:),'unit','normalized');
end

if ~isempty(badChannels)
    showElectrodeLocations(electrodeGridPos,badChannels,'r',[],1,0,gridType,subjectName);
    text(hElectrodes,-0.35,1,'Bad','color','r','unit','normalized');
end

% Get main plots and message handles

hERP = getPlotHandles(numElectrodeGroups,1,[0.05 0.05 0.1 0.6]);
hFR  = getPlotHandles(numElectrodeGroups,1,[0.175 0.05 0.1 0.6]);
hDeltaPSD = getPlotHandles(numElectrodeGroups,1,[0.3 0.05 0.1 0.6]);
hDeltaTF  = getPlotHandles(numElectrodeGroups,numConditions,[0.425 0.05 0.55 0.6]);

uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],'Style','text',...
    'String',[subjectName expDate protocolName],'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%% Get data from  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numGoodElectrodes = length(goodElectrodes);
allData = cell(1,numGoodElectrodes);

for i=1:numGoodElectrodes
    channelString = ['elec' num2str(goodElectrodes(i))];
    disp(['Getting data: ' num2str(i) ' of ' num2str(numGoodElectrodes) ', ' channelString]);

    allData{i} = getSpikeLFPDataSingleChannel(subjectName,expDate,protocolName,folderSourceString,channelString,0,gridType,[],referenceChannelString,badTrialNameStr,useCommonBadTrialsFlag);
end        

colormap jet;
colorNames = jet(numConditions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)

        s=get(hSigma,'val');
        f=get(hSpatialFreq,'val');
        o=get(hOrientation,'val');
        c=get(hContrast,'val');

        signalRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        freqRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];

        blRange = [str2double(get(hBaselineMin,'String')) str2double(get(hBaselineMax,'String'))];
        stRange = [str2double(get(hStimPeriodMin,'String')) str2double(get(hStimPeriodMax,'String'))];
        zRange = [str2double(get(hZMin,'String')) str2double(get(hZMax,'String'))];

        removeERPFlag = 1;

        for iCond = 1:numConditions

            if protocolType == 1
                a = iCond; e = 1; t = 1;
            elseif protocolType == 2
                a = 1; e = iCond; t = 1;
            end

            for iGroup = 1:numElectrodeGroups

                % Get data from electrodes
                tmpElectrodes = electrodeGroupList{iGroup};
                
                if ~isempty(tmpElectrodes)

                    numTmpElectrodes = length(tmpElectrodes);

                    if numTmpElectrodes > 1 % Combine across electrodes

                        tmpERPData = []; tmpFRData = []; tmpDeltaPSD = []; tmpDeltaTF = [];

                        for k = 1:numTmpElectrodes
                            tmpData = getDataGRF(allData{tmpElectrodes(k)==goodElectrodes},a,e,s,f,o,c,t,blRange,stRange,removeERPFlag);

                            tmpERPData = cat(1,tmpERPData,tmpData.erp);
                            tmpFRData = cat(1,tmpFRData,tmpData.frVals);
                            tmpDeltaPSD = cat(1,tmpDeltaPSD,tmpData.deltaPSD');
                            tmpDeltaTF = cat(3,tmpDeltaTF,tmpData.deltaTF);
                        end

                        erpData = mean(tmpERPData,1);
                        frData = mean(tmpFRData,1);
                        deltaPSD = mean(tmpDeltaPSD,1);
                        deltaTF = squeeze(mean(tmpDeltaTF,3));
                        
                    else
                        tmpData = getDataGRF(allData{tmpElectrodes==goodElectrodes},a,e,s,f,o,c,t,blRange,stRange,removeERPFlag);
                        erpData = tmpData.erp;
                        frData = tmpData.frVals;
                        deltaPSD = tmpData.deltaPSD';
                        deltaTF = tmpData.deltaTF;
                    end

                    % Plot data
                    plot(hERP(iGroup),tmpData.timeVals,erpData,'color',colorNames(iCond,:)); hold(hERP(iGroup),'on');
                    plot(hFR(iGroup),tmpData.frTimeVals,frData,'color',colorNames(iCond,:)); hold(hFR(iGroup),'on');

                    plot(hDeltaPSD(iGroup),tmpData.freqST,deltaPSD,'color',colorNames(iCond,:)); hold(hDeltaPSD(iGroup),'on');
                    plot(hDeltaPSD(iGroup),tmpData.freqST,zeros(1,length(deltaPSD)),'color','k');
                    
                    pcolor(hDeltaTF(iGroup,iCond),tmpData.timeTF,tmpData.freqTF,deltaTF'); shading(hDeltaTF(iGroup,iCond),'interp');
                    clim(hDeltaTF(iGroup,iCond),zRange); axis(hDeltaTF(iGroup,iCond),[signalRange freqRange]);
                end
            end
        end

        % Rescale plots to same scale
        rescalePlots(hERP,[signalRange getYLims(hERP)]);
        rescalePlots(hFR,[signalRange getYLims(hFR)]);
        rescalePlots(hDeltaPSD,[freqRange getYLims(hDeltaPSD)]);

        % Titles
        title(hERP(1),'ERPs');
        title(hFR(1),'Firing Rates');
        title(hDeltaPSD(1),'\DeltaPSD (dB)');
        for iCond = 1:numConditions
            title(hDeltaTF(1,iCond),num2str(condVals(iCond)),'color',colorNames(iCond,:));
        end

        for iGroup = 1:numElectrodeGroups
            ylabel(hERP(iGroup),groupNameList{iGroup},'color',colorNamesElectrodeGroups(iGroup,:));
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleZ_Callback(~,~)

        zRange = [str2double(get(hZMin,'String')) str2double(get(hZMax,'String'))];
        rescaleZPlots(hDeltaTF,zRange);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function rescaleData_Callback(~,~)
        
        signalRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        freqRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];

        % Rescale plots to same scale
        rescalePlots(hERP,[signalRange getYLims(hERP)]);
        rescalePlots(hFR,[signalRange getYLims(hFR)]);
        rescalePlots(hDeltaPSD,[freqRange getYLims(hDeltaPSD)]);
        rescalePlots(hDeltaTF,[signalRange freqRange]);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        claPlots(hERP);
        claPlots(hFR);
        claPlots(hDeltaPSD);
        claPlots(hDeltaTF);

        function claPlots(plotHandles)
            [numRow,numCol] = size(plotHandles);
            for ii=1:numRow
                for jj=1:numCol
                    cla(plotHandles(ii,jj));
                end
            end
        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;

for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims=[yMin yMax];
end
function rescalePlots(plotHandles,axisLims)
[numRow,numCol] = size(plotHandles);

for i=1:numRow
    for j=1:numCol
        axis(plotHandles(i,j),axisLims);
    end
end
end
function rescaleZPlots(plotHandles,caxisLims)
[numRow,numCol] = size(plotHandles);

for i=1:numRow
    for j=1:numCol
        clim(plotHandles(i,j),caxisLims);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [electrodeGroupList,groupNameList,goodElectrodes] = getElectrodeGroups(subjectName,gridType,electrodeArray,stimulationElectrode,badChannels)

% get highRMSelectrodes
tmp = load([subjectName gridType 'RFData.mat']); % Get RF data
if strcmp(subjectName,'dona')
    highRMSElectrodes = tmp.highRMSElectrodes;
    highRMSElectrodes = highRMSElectrodes(highRMSElectrodes<=48); % Only V1
end
goodElectrodes = setdiff(highRMSElectrodes,badChannels);

[r1,c1] = find(electrodeArray==stimulationElectrode);

distances = [];
for i=1:length(goodElectrodes)
    [r2,c2] = find(electrodeArray==goodElectrodes(i));
    distances = cat(2,distances,0.4*sqrt((r1-r2).^2 + (c1-c2).^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Now pool by distance %%%%%%%%%%%%%%%%%%%%%%%%%%%
distanceRangeList = 0:0.5:2.5;
numDistanceRangeList = length(distanceRangeList);
electrodeGroupList = cell(1,numDistanceRangeList);
groupNameList = cell(1,numDistanceRangeList);
for i=1:numDistanceRangeList-1
    electrodeGroupList{i} = goodElectrodes(intersect(find(distances>=distanceRangeList(i)),find(distances<distanceRangeList(i+1))));
    groupNameList{i} = [num2str(distanceRangeList(i)) '<=d<' num2str(distanceRangeList(i+1))];
end
electrodeGroupList{numDistanceRangeList} = goodElectrodes(distances>=distanceRangeList(numDistanceRangeList));
groupNameList{numDistanceRangeList} = ['d>=' num2str(distanceRangeList(numDistanceRangeList))];

end
function outString = getStringFromValues(valsUnique,decimationFactor)

if isscalar(valsUnique)
    outString = convertNumToStr(valsUnique(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique)
        outString = cat(2,outString,[convertNumToStr(valsUnique(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%
% load Data
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract)

x=load(fullfile(folderExtract,'parameterCombinations.mat'));
parameterCombinations=x.parameterCombinations;
aValsUnique=x.aValsUnique;
eValsUnique=x.eValsUnique;

if ~isfield(x,'sValsUnique')
    sValsUnique = x.rValsUnique/3;         
else
    sValsUnique=x.sValsUnique;
end

fValsUnique=x.fValsUnique;
oValsUnique=x.oValsUnique;

if ~isfield(x,'cValsUnique')
    cValsUnique=[];
else
    cValsUnique=x.cValsUnique;
end

if ~isfield(x,'tValsUnique')
    tValsUnique=[];
else
    tValsUnique=x.tValsUnique;
end
end
function impedanceValues = getImpedanceValues(fileName)
x=load(fileName);
if isfield(x,'impedanceValues')
    impedanceValues = x.impedanceValues;
elseif isfield(x,'electrodeImpedances')
    impedanceValues = x.electrodeImpedances;
else
    disp('Impedance information is not available');
    impedanceValues = [];
end
end