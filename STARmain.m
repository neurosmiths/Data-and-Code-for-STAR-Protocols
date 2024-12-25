% Load downsampled LFP data and MUA times from STARmain_data.mat **
load('STARmain_data_filtered.mat'); % Load only the LFP and Fs from the data file

% Smoothing LFP data: convolve LFP with a Gaussian kernel using the "smoothdata" matlab function
LFP = -smoothdata(-LFP','gaussian',Fs./200)';

%% Action Potentials

% determine how many units are in the data
APtimes = {};
chans = unique(ChanUnitTimestamp(:,1));
for ch = 1:length(chans)
    nChanUnits = max(unique(ChanUnitTimestamp(ChanUnitTimestamp(:,1)==chans(ch),2)));
    for un = 1:nChanUnits
        % number of units per channel
        nUnits = size(APtimes,1);
        % put action potential times per unit into an array
        unitTimes = {ChanUnitTimestamp(ChanUnitTimestamp(:,1)==chans(ch) & ChanUnitTimestamp(:,2)==un,3)'};
        spikes.channel(ch).unit(un).times = cell2mat(unitTimes);
        % put unit times per unit in a cell array that is nUnits long
        APtimes(nUnits+1,1) = unitTimes;
    end
end
% update final number of units.
nUnits = length(APtimes);
timeSize = max(ChanUnitTimestamp(:,3)); %  finds the maximum timestamp in the third column of ChanUnitTimestamp, used to determine the time duration of the recording

%% MUA and AP time detection

% Smoothing MUA and AP data to produce mean firing rate
[rAP,rSmoAP,tAP] = psthBins(cell2mat(APtimes').*3e4, 25, 3e4, length(APtimes),tEnd);
[rMUA,rSmoMUA,tMUA] = psthBins(cell2mat(MUAtimes)*3e4, 25, 3e4, length(MUAtimes),tEnd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Neural Event Detection (Bursts) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[peaks, locs] = findpeaks(rSmoMUA, 'MinPeakProminence', 5); % Find peaks (bursts) in the mean firing rate with min peak prominence threshold criteria
length(peaks); % Gives the # of bursts detected

% time base (creating different time vectors using the linspace function)
tSec = linspace(0,tEndSecs,length(LFP)); % LFP
tSecAP = linspace(0,tEndSecs,length(tAP)); % sorted units
tSecMUA = linspace(0,tEndSecs,length(tMUA)); % MUA
peakTimes = tSecMUA(locs); % Get the burst times in samples

winwidth = min(diff(locs)); % Define the window size, % approximately half of the median minimum difference in time between bursts. The medium min diff was used to calculate the window around each burst/peak = winwidth
% Looping over detected peaks
whichAnalysis = 'point'; % 'phase' or 'point'
%for pk = 26   % For when you want to pick a single burst to look at
for pk = length(locs):-1:1 % For when you want to analyze all of the detected bursts
    % Check if the window goes out of bounds
    if locs(pk) - winwidth < 1 || locs(pk) + winwidth > length(tSecMUA)
        % Skip this peak if it goes out of bounds
        continue;
    end

    % Define the window start and end
    startIdx = locs(pk) - winwidth;
    endIdx = locs(pk) + winwidth;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Traveling Wave Classification %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch whichAnalysis
        case {'point'}   % point = multilinear regression
            % Making a cell array with all of the burst-associated
            % multiunit action potential times in order to fit them to
            % a plane using the multilinear regression methods outlined
            % in Liou et al. (2017)
            for chs = length(MUAtimes):-1:1
                % Ensure MUAtimes{chs} is within valid bounds of tSecMUA
                muaR{chs} = MUAtimes{chs}(MUAtimes{chs} > tSecMUA(startIdx) & MUAtimes{chs} < tSecMUA(endIdx));
            end

            % Looking at a time period between this window [tSecMUA(startIdx) tSecMUA(endIdx)]
            LFPtimeIdcs = (tSec >= tSecMUA(startIdx) & tSec <= tSecMUA(endIdx));
            burstTime = tSec(LFPtimeIdcs);

            % Looping over LFP channels.
            for chs = size(LFP, 1):-1:1
                whichLFPfeat = 'negPeak'; % 'maxDescent' or 'negPeak'
                switch whichLFPfeat
                    case {'maxDescent'}
                        % Maximal descent
                        [amps, times] = min(diff(LFP(:, LFPtimeIdcs), 1, 2), [], 2);
                    case {'negPeak'}
                        % Negative peak
                        [amps, times] = min(LFP(:, LFPtimeIdcs), [], 2);
                end
                lfpR = num2cell(burstTime(times'));
            end
    end

    % converting the 2D array map into a 96 X 2 matrix with cm units, P.
    for kgl=length(pt)
        if pt==1
            map = electrodepinout; %location of the channels, corners are empty
            for el = numel(chans):-1:1
                [P(el,1),P(el,2)] = find(map == chans(el));

            end
            % original unscaled data is preserved in POL, while Plfp contains the scaled values ready for further use. This approach ensures that the original data is not lost and can be referenced if needed later on.
            POU=P;
            P = P*0.04; % cm  P= electrode locations for APs, .04 = space between each channel in CM

            LFPchans = 1:96;  % for all chans   Again now for LFPs, everychannels has an LFP while the APs may not be there for each channel
            for elt = size(LFP,1):-1:1
                [Plfp(elt,1),Plfp(elt,2)] = find(map == LFPchans(elt));

            end
            % original unscaled data is preserved in POL, while Plfp contains the scaled values ready for further use. This approach ensures that the original data is not lost and can be referenced if needed later on.
            POL = Plfp;
            Plfp = Plfp*0.04;

        elseif pt==2

            map = electrodepinout;

            chans = sort(map(map>0)); % for all chans, finding any electrode map that has a number

            for el = length(MUAtimes):-1:1
                [P(el,1),P(el,2)] = find(map == chans(el));
            end
            POU=P;

            P = P*0.04; % cm

            LFPchans = find(~badChannels);  % for all chans
            for elt = size(LFP,1):-1:1
                [Plfp(elt,1),Plfp(elt,2)] = find(map == LFPchans(elt));

            end
            POL = Plfp;
            Plfp = Plfp*0.04;

        end
    end

    %%%%%%%%% Figure 5 %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    whichLossFunction = 'L1'; % can be 'L1' or 'L2'
    % calculate traveling wave velocity vector from units.
    plotregressions = true; %sets a flag to enable plotting of regression results
    if sum(~cellfun(@isempty,muaR),2) >= 30
        [unitWaves(pk).beta, unitWaves(pk).V, unitWaves(pk).p] = ...0
            SpatialLinearRegression(muaR',P,'switch_plot',plotregressions,'Lossfun',whichLossFunction,'Fs',Fs);
        if plotregressions
            saveas(gcf,[whichLossFunction 'regularized_regressionplotUNITS_burst' num2str(pk) '_patient' num2str(pt)  '.pdf']);
            close(gcf);
        end

        % calculate traveling wave velocity vector from LFP times.
        [lfpWaves(pk).beta, lfpWaves(pk).V, lfpWaves(pk).p] = ...
            SpatialLinearRegression(lfpR,Plfp,'switch_plot',plotregressions,'Lossfun',whichLossFunction,'Fs',Fs);
        if plotregressions
            saveas(gcf, ([whichLossFunction 'regularized_regressionplotLFP_burst' num2str(pk) '_patient' num2str(pt) '.pdf']));
            close (gcf)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Traveling wave feature characterization %%%%%
        %%%%%%%%%%%%%% Figures 3, 6 and 8 %%%%%%%%%%%%%%%%%

        % per-peak LFP velocities (i.e., pseudoinverse of the first two linear regression model betas)
        unitV = pinv(unitWaves(pk).beta(1:2));
        lfpV = pinv(lfpWaves(pk).beta(1:2));

        %% saving data
        % saving traveling wave data for each peak as POLAR coordinates:'V' defines the velocity and 'angle' defines the direction
        lfpWaves(pk).V = hypot(lfpV(1),lfpV(2));
        unitWaves(pk).V = hypot(unitV(1),unitV(2));
        lfpWaves(pk).angle = atan2(lfpV(2),lfpV(1));
        unitWaves(pk).angle = atan2(unitV(2),unitV(1));

        % saving traveling wave data for each peak as CARTESIAN cooridnates
        lfpWaves(pk).carts = lfpV;
        unitWaves(pk).carts = unitV;

    end
    % updating user
    updateUser('cacluclated traveling waves for detected peaks:',pk,10,length(locs))

    % plotting individual burst features here.
    burstPLT = false;
    if burstPLT
        cMap = copper(size(LFP,1));
        %Units
        for ix = 1:length(POU) % 1x96 cell array now put into a matrix instread of a cell array, taking the location of each row of P = 71x2
            respMatunit(POU(ix,1), POU(ix,2)) = median(muaR{ix});
        end

        for ix = 1:length(lfpR) % 1x96 cell array now put into a matrix instead of a cell array, taking the location of each row of P = 71x2
            respMatLFP(POL(ix,1), POL(ix,2)) = lfpR{ix};
        end

        %Change the dynamic range to visualize a narrower time window
        unitsclim= [min(cell2mat(muaR)) max(cell2mat(muaR))];

        % Making the corners NaNs
        respMatunit([1, 1, 10, 10], [1, 10, 1, 10]) = NaN;
        respMatLFP([1, 1, 10, 10], [1, 10, 1, 10]) = NaN;

        % Figure 3A: Visualization of detected bursts (peaks)
        subplot(3,3,1)
        data = rSmoMUA;
        plot(data) % Plot the original data
        hold on % Hold the current plot to add peaks
        plot(locs, peaks, 'ro') % Plot the peaks as red circles on the same plot
        hold off % Release the hold
        title('Detected Peaks in rSmoMUA')
        xlabel('Time (samples)')
        ylabel('Amplitude')
        legend('rSmoMUA', 'Peaks')

        % Figure 3B: visualize LFP during burst (and mean across chans in black)
        subplot(3,3,2)
        hold on
        for chs = 1:size(LFP,1)
            plot(tSec(LFPtimeIdcs), LFP(chs,LFPtimeIdcs),'color',cMap(chs,:),'linewidth',0.5)
        end
        plot(tSec(LFPtimeIdcs), mean(LFP(:,LFPtimeIdcs)),'color','b')
        scatter(burstTime(times),amps,5,[0 1 0],'filled')
        hold off
        axis tight square
        xs = xlim;
        ylabel('LFP (uV)')

        % Figure 3C: raster plot during burst.
        % organizing spike times for each burst.
        for hc = 1:length(MUAtimes)
            bstTimes{hc} = MUAtimes{hc}(MUAtimes{hc}>xs(1) & MUAtimes{hc}<xs(2)) - xs(1);
        end
        subplot(3,3,4)
        [~,I]=sort(cellfun(@median, bstTimes))
        plotSpikeRaster(bstTimes(I)','PlotType','vertline')
        axis square
        xlim(xs-xs(1))
        xlabel('time (s)')
        ylabel('microelectrodes')
        clear bstTimes

        % Figure 3E: Heaviside Function
        binFactor = 40;
        smoothWindowSize = Fs./binFactor;
        smoothedLFPdata = LFP(:, LFPtimeIdcs);
        nChans = size(smoothedLFPdata,1);

        % finding mins...
        [localMinTimeIdcs,localMinAmps] = islocalmin(smoothedLFPdata,2,'MinProminence',median(std(abs(smoothedLFPdata),[],2)),'MaxNumExtrema',nChans*3);
        [absoluteMinAmps,absoluteMinTimeIdcs] = min(smoothedLFPdata,[],2);

        % minimum detection parameters.
        edges = 0:0.02:1;
        weightingFunction = [0.1*ones(1,(floor(length(edges)/2))-1) linspace(1,0.5,ceil(length(edges)/2))];
        [N,edges,bindices] = histcounts(tSec(LFPtimeIdcs),edges);
        centers = (edges(1:end-1) + edges(2:end))/2;

        % finding the local minima in the largest weighted histogram bin
        [~,maxHisto] = max(weightingFunction.*N);

        % Plot Heaviside
        subplot(3,3,8)
        hold on
        plot(centers,weightingFunction.*N,'color',rgb('gray'))
        scatter(tSec(localMinTimeIdcs),localMinAmps(localMinTimeIdcs)...
            ,6,[0 0 0],'filled')
        scatter(tSec(absoluteMinTimeIdcs),abs(absoluteMinAmps)...
            ,4,rgb('lightseagreen'),'filled')
        hold off
        xlabel('time (s)')
        ylabel('local minimum prominence')
        xlim([0 1])
        axis square
        legend('weighted histogram','local minima','absolute minima')
        title('minima detections: scatter')
        xlabel('time (s)')
        ylabel('peak prominence')
        zlabel('local minima per bin')
        axis xy
        title('minima detections: histogram')
        sgtitle(sprintf('Burst # %d', pk));

        % stats summary
        subplot(3,3,3)
        text(0,0.6,'stats summary:','fontweight','bold')
        text(0,0.4,sprintf('LFP p-value = %.2f',lfpWaves(pk).p))
        text(0,0.2,sprintf('unit p-value = %.2f',unitWaves(pk).p))
        axis off

        % Figure 6A: Compass Plot of Unit and LFP Wave direction and speed during a burst
        subplot(3,3,5)
        compass(lfpV(1),lfpV(2),'r')
        hold on
        compass(unitV(1),unitV(2),'b')
        hold off
        legend('LFP wave','unit wave')

        % Figure 8A: Unit Wave propagating across UMA Visualization
        subplot(3,3,6)
        if (unitsclim(2)-unitsclim(1))==0
            imagesc(respMatunit);
        else
            imagesc(respMatunit,unitsclim);
        end
        colormap(copper);
        axis xy square off;
        xlabel('Microelectrodes');
        ylabel('Microelectrodes');
        h = colorbar;
        ylabel(h, 'Median Unit timmings (s)');

        % Figure 8A: LFP Wave propagating across UMA Visualization
        subplot(3,3,7)
        imagesc(respMatLFP,[min(min(respMatLFP)) max(max(respMatLFP))]); % [min(min(respMatLFP)) max(max(respMatLFP))]); %   [min(tSec(LFPtimeIdcs)) max(tSec(LFPtimeIdcs))]
        colormap(copper);
        axis xy square off;
        xlabel('Microelectrodes');
        ylabel('Microelectrodes');
        h = colorbar;
        ylabel(h, 'Maximal Descent LFP timings (s)');

        % saving figures.
        orient(gcf,'portrait')
        halfMaximize(gcf,'page')

        saveas(gcf,[whichLossFunction 'regularized_burst', num2str(pk), '_patient', num2str(pt),'.pdf'])
        close (gcf)
    end

end

% Extract angles (directions), velocities, and p-values from structures
lfpDir = [lfpWaves.angle];
unitDir = [unitWaves.angle];
lfpV = [lfpWaves.V];
unitV = [unitWaves.V];
lfp_p_values = [lfpWaves.p];
unit_p_values = [unitWaves.p];
p_threshold = 0.05;

% Filter for significant LFP and unit waves based on p-values and remove outliers
significant_lfp_indices = lfp_p_values < p_threshold & ~isoutlier(lfpV);
significant_unit_indices = unit_p_values < p_threshold & ~isoutlier(unitV) & unitV;

% Filter peak times, LFP, and unit speeds
peakTimes = tSecMUA(locs);
filtered_peakTimes = peakTimes(significant_unit_indices);
filtered_unitV = unitV(significant_unit_indices);

% Further filter based on time (e.g., > 150 seconds)
time_filtered_indices = filtered_peakTimes > 0;
TWtimes = filtered_peakTimes(time_filtered_indices);
filtered_unitV = filtered_unitV(time_filtered_indices);

% Apply filtering to LFPs
filtered_lfp_peakTimes = peakTimes(significant_lfp_indices);
filtered_lfpV = lfpV(significant_lfp_indices);
filtered_lfpDir = lfpDir(significant_lfp_indices);

% Apply filtering to Units
filtered_unitDir = unitDir(significant_unit_indices);

% Ensure LFP and unit directions are in the range [0, 2*pi]
%filtered_unitDir = wrapTo2Pi(unitDir(significant_unit_indices));
%filtered_unitDir = filtered_unitDir(time_filtered_indices);

% Create figure with tiled layout
t = tiledlayout(3, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

% Plot LFP p-value histogram
nexttile(1);
histogram(lfp_p_values, 20, 'FaceColor', 'k');
line([p_threshold, p_threshold], [0, 50], 'LineStyle', '--', 'Color', 'r');
axis square;
xlabel('p-values');
ylabel('Count');
title(sprintf('%d of %d LFP waves significant (%.1f%%)', sum(lfp_p_values <= p_threshold), length(lfp_p_values), 100 * mean(lfp_p_values <= p_threshold)));

% Plot unit p-value histogram
nexttile(2);
histogram(unit_p_values, 20, 'FaceColor', 'k');
line([p_threshold, p_threshold], [0, 30], 'LineStyle', '--', 'Color', 'r');
axis square;
xlabel('p-values');
ylabel('Count');
title(sprintf('%d of %d unit waves significant (%.1f%%)', sum(unit_p_values <= p_threshold), length(unit_p_values), 100 * mean(unit_p_values <= p_threshold)));

% LFP direction polar histogram
nexttile(3);
polarhistogram(filtered_lfpDir, 18, 'DisplayStyle', 'stairs', 'EdgeColor', 'r');
ax = gca;
ax.ThetaTick = [0 180];
ax.ThetaTickLabel = {'0', '\pi'};
title('LFP Waves');

%  Unit direction polar histogram
nexttile(4);
polarhistogram(filtered_unitDir, 18, 'DisplayStyle', 'stairs', 'EdgeColor', 'b');
ax = gca;
ax.ThetaTick = [0 180];
ax.ThetaTickLabel = {'0', '\pi'};
title('Unit Waves');

% Figure 6B:Combined LFP and unit direction histogram
nexttile(5);
polarhistogram(filtered_lfpDir, 18, 'DisplayStyle', 'stairs', 'EdgeColor', 'r');
hold on;
polarhistogram(filtered_unitDir, 18, 'DisplayStyle', 'stairs', 'EdgeColor', 'b');
ax.ThetaTick = [0 180];
ax.ThetaTickLabel = {'0', '\pi'};
title('Combined LFP and Unit Waves');
hold off;

% Figure 6C:Boxplot for LFP and unit speeds
nexttile(6);
hold on;
boxchart(ones(size(filtered_lfpV)) - 0.2, filtered_lfpV, 'BoxWidth', 0.1, 'BoxFaceColor', 'r');
boxchart(ones(size(filtered_unitV)) + 0.2, filtered_unitV, 'BoxWidth', 0.1, 'BoxFaceColor', 'b');
axis square;
xticks([1]);
xticklabels({'Speeds'});
ylabel('Speed (cm/s)');
title('LFP and Unit Speed Distribution');

% Add global title
sgtitle([ptID, ', file', num2str(nFiles), ', ', whichLossFunction, ', ', whichLFPfeat, ', Burst Summary Plots']);

% Save figure
filename = [ptID, '_', num2str(nFiles), '_', whichLossFunction, '_', whichLFPfeat, '_', 'allBursts.pdf'];
print(gcf, filename, '-dpdf', '-painters', '-fillpage');

%%%%%%%%%% Figure 2 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting a data summary
Sm = figure;
figure(Sm)

% Figure 2A: plot LFP
ax(1) = subplot(3,1,1); %The subplot handle is stored in the array ax, index 1.
plot(tSec,mean(LFP),'k') % tSec is a vector containing the time points and LFP is a matrix where each row represents LFP data at a different time point
xlabel('time (s)')
ylabel('LFP (uV)')

% Figure 2B:plot MUA
ax(4) = subplot(3,1,2);
plotSpikeRaster(MUAtimes','PlotType','VertLine');
ylabel('microelectrodes')
title('multiunit event times')

% Figure 2C:plot mean firing rate
ax(3) = subplot(3,1,3);
hold on
plot(tSecAP,rSmoAP,'color','r')
plot(tSecMUA,rSmoMUA,'color','b')
scatter(tSecMUA(locs),peaks,12,rgb('fuchsia'))
legend('sorted APs','sorted MUAs','detected peaks')
hold off
axis tight
ylabel('z-scored firing rate')

nFiles = num2str(nFiles);
sgtitle([ptID,', file', nFiles, ', Raw Signal Plots']);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperType', 'A4');
orient(gcf, 'landscape');
saveas(gcf,[whichLossFunction 'regularized_regressionplotUNITS_burst' num2str(pk) '_patient' num2str(pt)  '.pdf']);
close (gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Figure 7-Bimodality Analysis %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize results containers for patient data
circStatsCluster1 = struct();
circStatsCluster2 = struct();

sigDirections = lfpDir;
% use unitDir for unit analyzing unit waves
% sigDirections =unitDir

% Reduce precision to work with fitmvmdist.
sigDirections = (round(sigDirections * 10^2) / 10^2)';

% Figure 7A: Evaluate normality for overall direction distribution 
[p_valAll, k_statAll, k_permsAll] = evaluateNormality(sigDirections, 10000, 100, 'starmain_overall_IEDdistFit.pdf');

% Step: Fit von Mises Mixture Model (vMM) to direction data (Figure 7B: Clustered vMM distributions)
threshold = [0.3 0.7];
GMM = fitmvmdist(sigDirections, 2);
[clusters, nLogL, PostP] = cluster(GMM, sigDirections);

% Identify shared regions for clusters
idxBothSharedDiag = PostP(:, 1) >= threshold(1) & PostP(:, 1) <= threshold(2);
numInBoth = sum(idxBothSharedDiag);

% Figure 7B: Fit von Mises distribution for each cluster 
circStatsCluster1 = circ_stats(sigDirections(clusters == 1 | idxBothSharedDiag));
circStatsCluster2 = circ_stats(sigDirections(clusters == 2 | idxBothSharedDiag));

% Figure 7C: Evaluate normality for each cluster
[p_val_1, k_stat_1, k_perm_1] = evaluateNormality(sigDirections(clusters == 1 | idxBothSharedDiag), 10000, 100, 'starmain_cluster1_IEDdistFit.pdf');
[p_val_2, k_stat_2, k_perm_2] = evaluateNormality(sigDirections(clusters == 2 | idxBothSharedDiag), 10000, 100, 'starmain_cluster2_IEDdistFit.pdf');

% running simulations
a1 = - pi; a2 = pi;

% number of simulations to run.
nSims = 100;

% size of permutation samples. irrelevant if fixed to n + m in permutation_kuipertest.
permutationSizes = [50 100 200 500];

% looping over simulations and permutation sizes
pVals = ones(length(permutationSizes),nSims);
for ms = 1:nSims
    fprintf('\n doing simulation %d of 100...',ms)
    for sm = 1:length(permutationSizes)
        angles = (a2-a1).*rand(500,1) + a1;

        [p_val,k_stat,k_perm] = evaluateNormality(angles,10000,permutationSizes(sm));

        pVals(sm,ms) = p_val;

        if ms==nSims
            fprintf('\nFalse positive rate for permutation test with %d angles = %.2f',permutationSizes(sm),sum(pVals(sm,:)<0.05)./nSims)
        end
    end
end

% saving the data.
save('STARsimulations_100.mat')

% subsequent analysis of the simulation data
fprintf('\nOverall false positive rate = %.2f',sum(sum(pVals<0.05))./numel(pVals))

imagesc(pVals)


%%%%%%%%%% Function Section %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Beta, V, pValue] = SpatialLinearRegression(Input,P,varargin)
% [Beta, V, pValue] = SpatialLinearRegression(T,P,varargin)
%
% Estimating traveling wave velocity by multivariate linear regression.
%
% This function is used for the following methods
%
% - Multiunit estimator
% - Negative peak estimator
% - Maximal descent estimator
% - High gamma envelope estimator
%
% Please see the following description to cater the data format for each
% estimator
%
% Input: 1) Input:
%           ** Important **
%           Input type can be cell or numeric matrix, depends on which
%           estimator you are going to use.
%
%           * [Data type: cell] - used for
%                               * multiunit estimator
%                               * negative peak estimator
%                               * maximal descent estimator
%
%           In this case, Input is 'T' (information of event timing)
%           Input needs to be n_channel by 1 cell, which each cell
%           corresponds to the information from one recording channel.
%           Within a cell, it contains a numeric vector which encodes
%           timing of events (multiunit spike, negative peak, or maximal
%           descent) happen.
%
%           For example: channel 1 detects 3 multiunit spikes at time 5,
%                        7, 10 during an ictal discharge.
%                        Input{1} = [5, 7, 10];
%
%           For example, channel 1 detects its negative peak at time 12
%                        Input{1} = 12;
%
%           * [Data type: numeric matrix] - used for
%                                         * high gamma envelope estimator
%
%           In this case, Input encodes the instantaneous power.
%           Input needs to be a n_channel by n_time matrix.  Every entry
%           of Input needs to be strictly positive.  In other words, Input
%           here is treated as 'W' (weight) for linear regression.
%
%           For example, channel 1 has instantaneous power ramped up and
%           down within 10 time steps of an ictal discharge
%
%           Input(1,:) = [1 2 3 4 5 5 4 3 2 1]
%
%           The sampling rate (n_time, the 2nd dimension of Input matrix)
%           is set to be 1000 Hz by default.  You can change it by giving
%           optional input. For example
%
%           SpatialLinearRegression(Input,P,'Fs',500);
%
%        2) P: The position information. Format: an n_channel by 2 numetric
%           matrix.  For example, a 2 by 2 grid will be presented as
%
%           P = [1 1;
%                1 2;
%                2 1;
%                2 2];
%
%        3) [Optional inputs]
%           SpatialLinearRegression(Input,P,'Parameter_Name',parameter);
%
%           switch_plot: (1/0) - plot the result or not
%           n_shuffle: n of shuffle test if shuffle test is called
%           Fs: sampling frequency for time series data - default:1000 Hz
%           Lossfun: Define the loss function in regression
%                    Options:
%                    'L2' - least squares
%                    'L1' - least absolute deviation
%                    Default: 'L2'
%
% Output: Beta: The estimated regression coefficients
%         V: The estimated traveling wave velocity
%         pValue: p Value of F-test or shuffle test
%
% Author: Jyun-you Liou
% Final update: 2016/11/30, version 1

p = inputParser;
addParameter(p,'alpha',0.05); % Significance level
addParameter(p,'switch_plot',0);
% In point process, it will plot data average at its physical position.
% In continuous process, it will return shuffle test results.
addParameter(p,'n_shuffle',100); % # of shuffle test in continuous process case
addParameter(p,'Fs',1000); % Only effective in time series data (Hz)
addParameter(p,'Lossfun','L2'); % Definition of loss function
parse(p,varargin{:});

if rank([ones(size(P,1),1), P]) <= size(P,2)
    error('The spatial information contains less dimensions than it shows.')
end

if iscell(Input)
    %% Point process case (multiunit, negative peak, maximal descent)
    n_channel = numel(Input);
    %P = mat2cell(P,ones(n_channel,1),size(P,2));
    P = mat2cell(P, ones(size(P, 1), 1), size(P, 2));
    for i = 1:n_channel
        Input{i} = Input{i}(:); % make sure the alignment is correct
        P{i} = repmat(P{i}(:)',[numel(Input{i}),1]);
    end
    P = cell2mat(P);
    P1 = P(:,1);
    P2 = P(:,2);
    X = P;
    Y = cell2mat(Input(:));
    % Check sample size
    [Xm,Xn] = size(X);
    if Xm<=Xn
        warning('Not enough event sample (under-determined system), not enough data to perform estimation.')
        Beta = nan(1,Xn);
        V = nan(1,Xn);
        pValue = nan;
        return;
    end

    % Linear regression - by loss function
    switch p.Results.Lossfun
        case 'L2'
            XC = [X,ones(Xm,1)];
            Beta = lscov(XC,Y);
            Fstat = ((var(Y,0)-var(Y-XC*Beta,0))/(Xn)) ...
                /(var(Y-XC*Beta,0)/(Xm-Xn-1));
            F_CDF = @(f,df1,df2) betainc(f/(f+df2/df1),df1/2,df2/2);
            try
                pValue = 1- F_CDF(Fstat,Xn,Xm-Xn-1); % F-test
            catch
                pValue = 1;
            end
        case 'L1'
            [Beta, E, E_s, pValue] = L1LinearRegression(X,Y); % Nested function
    end
elseif ismatrix(Input) && isnumeric(Input)
    %% Time series case
    W = Input;
    [n_channel, n_t] = size(W);
    T_vector = (1:n_t)/p.Results.Fs;
    P1 = P(:,1);
    P2 = P(:,2);
    P1 = repmat(P1,[1,n_t]);
    P2 = repmat(P2,[1,n_t]);
    P1 = P1(:);
    P2 = P2(:);
    % Prepare design matrix
    X = [P1, P2];
    % Prepare response matrix % weight
    Y = repmat(T_vector,[n_channel,1]);
    Y = Y(:);
    W = W(:);
    % Check sample size
    [Xm,Xn] = size(X);
    if Xm<Xn
        warning('Not enough event sample (under-determined system), not enough data to perform estimation.')
        V = nan(1,Xn);
        return;
    end
    switch p.Results.Lossfun
        case 'L2'
            [Beta, E, E_s, pValue] = L2LinearRegression(X, Y, W);
        case 'L1'
            [Beta, E, E_s, pValue] = L1LinearRegression(X, Y, W);
    end
end
V = pinv(Beta(1:2));

%% Plotting the results
if p.Results.switch_plot
    % Plot the regression plane
    figure(gcf)
    subplot(2,2,1)
    if iscell(Input)
        scatter3(P1,P2,Y,'filled');hold on;
    elseif ismatrix(Input) && isnumeric(Input)
        scatter3(P1,P2,Y,100*W/max(W(:)),'filled');hold on;
    end
    [P1U,P2U] = meshgrid(sort(unique(P1)),sort(unique(P2)));
    f = scatteredInterpolant(P1,P2,X*Beta(1:2) + Beta(end));
    Z = f(P1U,P2U);
    mesh(P1U,P2U,Z) %interpolated

    %   legend('Data','Regression fit');
    xlabel('cm');ylabel('cm');zlabel('Second');

    % Plot the projection along the velocity axis
    subplot(2,2,2)
    P_v_axis = X*Beta(1:2)/norm(Beta(1:2));
    if iscell(Input)
        plot(P_v_axis,Y,'.');
    elseif ismatrix(Input) && isnumeric(Input)
        f = scatteredInterpolant(P_v_axis,Y,W);
        P_v_axis_lin = linspace(min(P_v_axis),max(P_v_axis),n_channel);
        [P1U,P2U] = meshgrid(P_v_axis_lin,Y);
        Z = f(P1U,P2U);
        imagesc(P_v_axis_lin,Y,Z);
    end
    hold on;
    plot(P_v_axis,[P1,P2]*Beta(1:2)+Beta(end));
    title('Projection along the velocity vector');
    xlabel('cm');
    ylabel('Second');colormap('hot')

    % Plot shuffled residual error distribution if L1-regression is used
    if strcmpi(p.Results.Lossfun,'L1')
        subplot(2,2,3)
        hist(E_s);hold on;
        plot([E E],get(gca,'YLim'),'r');hold off
        title(['Mean residual error, p = ' num2str(pValue)]);
    end
end

%% Nested function - L2 regression

    function [B, E, E_s, L2_p_shuffle] = L2LinearRegression(X,Y,Weight)
        if nargin < 3
            Weight = 1;
        end
        % Determine size of predictor data
        [n,~] = size(X);
        % Add constant term
        X = [X,ones(n,1)];
        [B,~,E] = lscov(X,Y,Weight);
        % pValue
        if p.Results.n_shuffle
            % Shuffle the channel's spatial information
            E_s = nan(p.Results.n_shuffle,1);
            for i_shuffle = 1:p.Results.n_shuffle
                [~,~,E_s(i_shuffle)] = lscov(Shuffle_Position(X),Y,Weight);
            end
            L2_p_shuffle = sum(E_s<=E)/p.Results.n_shuffle;
        else
            E_s = [];
            L2_p_shuffle = nan;
        end
    end

%% Nested function - L1 regression

    function [B, E, E_s, L1_p_shuffle] = L1LinearRegression(X,Y,Weight)
        if nargin < 3
            Weight = 1;
        end
        % Determine size of predictor data
        [n,m] = size(X);
        % Add constant term
        X = [X,ones(n,1)];
        [B,E] = L1Optimize(X,Y,Weight);

        if p.Results.n_shuffle
            % Shuffle the channel's spatial information
            E_s = nan(p.Results.n_shuffle,1);
            for i_shuffle = 1:p.Results.n_shuffle
                [~,E_s(i_shuffle)] = L1Optimize(Shuffle_Position(X),Y,Weight);
            end
            L1_p_shuffle = sum(E_s<=E)/p.Results.n_shuffle;
        else
            E_s = [];
            L1_p_shuffle = nan;
        end
        function [Bout,Eout] = L1Optimize(Xin,Yin,Win)
            % Use least-squares fit as initial guess
            Bout = Xin\Yin;
            % Least squares regression
            BOld = Bout;
            BOld(end) = BOld(end) + 1e-5;
            iter = 1;
            % Force divergence % Repeat until convergence
            while (max(abs(Bout - BOld)) > 1e-6)
                % Move old coefficients
                BOld = Bout;
                % Calculate new observation weights (based on residuals from old coefficients)
                w = sqrt(Win ./ max(abs(Xin * BOld - Yin),1e-6));
                % Floor to avoid division by zero
                % Calculate new coefficients
                Bout = (repmat(w,[1 m+1]) .* Xin) \ (w .* Yin);
                iter = iter + 1;
                if iter > 30;break;end
            end
            % Report mean absolute deviation
            Eout = sum(w.* abs(Xin * BOld - Yin));
        end
    end

    function Xout = Shuffle_Position(Xin)
        [X_u,~,IB_X] = unique(Xin,'rows');
        Per = randperm(size(X_u,1));
        X_u_per = X_u(Per,:);
        Xout = X_u_per(IB_X,:);
    end

end

% UMA pinout
function map = electrodepinout()

map = [
    -1  96  93  92  90  88  85  83  81  -1
    95  63  94  91  89  87  86  84  80  79
    32  61  59  57  55  53  49  82  78  77
    30  64  60  58  56  51  47  45  76  75
    28  31  62  52  46  44  43  41  74  73
    26  29  21  54  50  42  40  39  72  71
    24  27  25  19  15  48  38  37  70  69
    22  20  23  13  17   5  36  35  68  67
    18  16  12  11   9   7  34  33  66  65
    -1  14  10   8   6   4   3   1   2  -1
    ];

end


function [p_val, k_stat, k_perm] = evaluateNormality(angles, nPerms, permSize, fName)
% This function will fit a VonMises distribution to the input angles
% and return a measure of how well the empirical distribution fits,
% based on a Kuiper test that evaluates significance using [nPerms] permutations
% of angles in groups of [permSize].
%
% Adding a fourth input argument will save an analysis figure with that file name
%

% Fitting the VonMises distribution
[thetahat, kappa] = circ_vmpar(angles);

% Now sampling from that distribution
nAngles = length(angles);
alpha = circ_vmrnd(thetahat, kappa, nAngles);

% Ensure permutationSize is valid
totalSampleSize = length(angles) + length(alpha);
if permSize > totalSampleSize
    permSize = totalSampleSize;
    fprintf('Adjusted permSize to %d to match total sample size.\n', permSize);
end

% Running the Kuiper test between the original and ideal distributions
fprintf('\nRunning permutation Kuiper test. This may take a moment.');
[p_val, k_stat, k_perm] = permutation_kuipertest(angles, alpha, nPerms, permSize, 90);

% Plotting the results if fName is provided
if nargin > 3
    figure(1);

    subplot(3, 2, 1);
    polarhistogram(angles, 18, 'DisplayStyle', 'stairs');
    title('Histogram of input angles');

    subplot(3, 2, 2);
    polarhistogram(alpha, 18, 'DisplayStyle', 'stairs');
    title('Histogram of generated angles');

    subplot(3, 2, [3 4]);
    H = histfit(k_perm);
    xlabel('Permutation statistics');
    title(sprintf('True data test stat = %.2f', k_stat));

    saveas(1, fName);
    close(1);
end
end


function [peaks] = find_inflections (timeseries, peaktype, interval, slope_threshold)
% function to find maxima and minima of time series

% inflections: finds all inflection points of the time series,
% ie both maxima and minima
% defined by slope changing sign
% returns index of points found
% optional argument "interval" specifies window of interest

peaks = [];
if nargin < 1
    warning ('Usage: find_inflections <time series vector> ["maxima"|"minima"|"both"],[<time interval of interest ([pstart end])>], [slope threshold]');
    return;
end

if nargin < 2 || isempty(peaktype)
    peaktype = 'both';
end

if nargin < 3 || isempty(interval)
    interval = [1 length(timeseries)];
end

if nargin < 4 || isempty(slope_threshold)
    slope_threshold = 0;
end

%num2str(timeseries)
firstderiv = diff(timeseries(interval(1):interval(2)));
positive_slope = find (firstderiv >= slope_threshold);
%num2str(positive_slope)
sign_slope = zeros (1,length(firstderiv));
sign_slope(positive_slope) = 1;
%num2str(sign_slope)
sign_changes = diff(sign_slope);

switch peaktype
    case 'both'
        % if finding both positive and negative inflections,
        % look for any nonzero sign changes
        sign_changes = abs(sign_changes);
        peaks = find(sign_changes > 0) + 1;
    case 'maxima'
        peaks = find(sign_changes < 0) + 1;
    case 'minima'
        peaks = find(sign_changes > 0) + 1;
end

end


function [] = halfMaximize(figureHandle,side)
% HALFMAXIMIZE maximizes figure in figureHandle to half the horizontal area
% of the currently used monitor.
%
% halfMaximize(gcf,'right')
%
% author: Elliot Smith
% Version Date: 20170530

% default side
if ~exist(side,'var')
    side = 'left';
end

switch side
    case {'left'}
        % maximizing
        screenDimensions = get(0,'Screensize');
        screenDimensions(3) = ceil(screenDimensions(3)./2);
        set(figureHandle, 'Position', screenDimensions); % Maximize figure.
    case {'right'}
        % maximizing
        screenDimensions = get(0,'Screensize');
        screenDimensions(3) = ceil(screenDimensions(3)./2);
        screenDimensions(1) = floor(screenDimensions(3)./2);
        set(figureHandle, 'Position', screenDimensions); % Maximize figure.
end
end

function [] = maximize(figureHandle)

% MAXIMIZE maximizes figure in figureHandle
%
% author: Elliot Smith
% Version Date: 20141020

% maximizing
screenDimensions = get(0,'Screensize');
set(figureHandle, 'Position', screenDimensions); % Maximize figure.
end



function [p_val, k_data, K_dist] = permutation_kuipertest(alpha1, alpha2, nPerms, permutationSize, vis_on)
% [pval, k, K] = permutation_kuipertest(alpha1, alpha2, nPerms, permuationSize, vis_on)
%
% The Kuiper two-sample permutation test tests whether the two samples
% differ significantly based on probabilities derived from a permutation
% distribution. The difference can be in any property, such as mean
% location and dispersion. It is a circular analogue of the
% Kolmogorov-Smirnov test.
%
% H0: The two distributions are identical.
% HA: The two distributions are different.
%
% Input:
%   alpha1    first sample (in radians)
%   alpha2    second sample (in radians)
%   res       resolution at which the cdf is evaluated
%   vis_on    display graph
%
% Output:
%   p_val       permutation p-value
%   k_data      test statistic for real data
%   K_dist      permutation distribution
%
% Modified from the Circular Statistics Toolbox for Matlab
% EHS20160521

% Update 2012
% By Marc J. Velasco and Philipp Berens, 2009
% velasco@ccs.fau.edu

if nargin < 5
    vis_on = 0;
end

res = 108;

n = length(alpha1(:));
m = length(alpha2(:));

% create cdfs of both samples
[phis1, cdf1, phiplot1, cdfplot1] = circ_samplecdf(datasample(alpha1, floor(permutationSize / 2)), res);
[foo, cdf2, phiplot2, cdfplot2] = circ_samplecdf(datasample(alpha2, floor(permutationSize / 2)), res); %#ok<ASGLU>

% maximal difference between sample cdfs
[dplus, gdpi] = max([0, cdf1 - cdf2]);
[dminus, gdmi] = max([0, cdf2 - cdf1]);

% calculate k-statistic
k_data = sqrt(n + m) * (dplus + dminus);

% find p-value
% parfor_progress(nPerms);
k_perm = NaN(1, nPerms);
for ps = 1:nPerms
    perm_samp = datasample([alpha1(:); alpha2(:)], permutationSize, 'Replace', false);

    % create cdfs of both samples
    [phis1, cdf1, phiplot1, cdfplot1] = circ_samplecdf(perm_samp(1:floor(permutationSize / 2)), res);
    [foo, cdf2, phiplot2, cdfplot2] = circ_samplecdf(perm_samp(floor(permutationSize / 2):end), res); %#ok<ASGLU>

    % maximal difference between sample cdfs
    [dplus, gdpi] = max([0, cdf1 - cdf2]);
    [dminus, gdmi] = max([0, cdf2 - cdf1]);

    % calculate k-statistic
    k_perm(ps) = sqrt(permutationSize) * (dplus + dminus);
    % parfor_progress
end
% parfor_progress(0);

% output values
n_stat_bins = 100;
[K_dist, perm_edges] = histcounts(k_perm, n_stat_bins);
p_val = sum(k_perm > k_data) / nPerms;

% visualize
if vis_on
    figure;
    hold on;
    histfit(k_perm, n_stat_bins);
    % line([k_data k_data], [0 20], 'linewidth', 2, 'color', 'k')
    hold off;
    xlabel('Test statistic');
    ylabel('Count');
    title(sprintf('True data test statistic = %d', k_data));
end
end


function [xPoints, yPoints] = plotSpikeRaster(spikes,varargin)
% PLOTSPIKERASTER Create raster plot from binary spike data or spike times
%   Efficiently creates raster plots with formatting support. Faster than
%   common implementations. Multiple plot types and parameters available!
%   Look at Parameters section below.
%
%   Inputs:
%       M x N logical array (binary spike data):
%           where M is the number of trials and N is the number of time
%           bins with maximum of 1 spike per bin. Assumes time starts at 0.
%       M x 1 cell of spike times:
%           M is the number of trials and each cell contains a 1 x N vector
%           of spike times. Units should be in seconds.
%
%   Output:
%       xPoints - vector of x points used for the plot.
%       yPoints - vector of y points used for the plot.
%
%   Parameters:
%       PlotType - default 'horzline'. Several types of plots available:
%           1. 'horzline' -     plots spikes as gray horizontal lines.
%           2. 'vertline' -     plots spikes as vertical lines, centered
%               vertically on the trial number.
%           3. 'scatter' -      plots spikes as gray dots.
%
%           ONLY FOR BINARY SPIKE DATA:
%           4. 'imagesc' -      plots using imagesc. Flips colormap so
%               black indicates a spike. Not affected by SpikeDuration,
%               RelSpikeStartTime, and similar timing parameters.
%           5. 'horzline2' -    more efficient plotting than horzline if
%               you have many timebins, few trials, and high spike density.
%               Note: SpikeDuration parameter DOES NOT WORK IF LESS THAN
%               TIME PER BIN.
%           6. 'vertline2' -    more efficient plotting than vertline if
%               you have few timebins, many trials, and high spike density.
%           Note: Horzline and vertline should be fine for most tasks.
%
%       FigHandle - default gcf (get current figure).
%           Specify a specific figure or subplot to plot in. If no figure
%           is specified, plotting will occur on the current figure. If no
%           figure is available, a new figure will be created.
%
%       LineFormat - default line is gray. Used for 'horzline' and
%           'vertline' plots only. Usage example:
%               LineFormat = struct()
%               LineFormat.Color = [0.3 0.3 0.3];
%               LineFormat.LineWidth = 0.35;
%               LineFormat.LineStyle = ':';
%               plotSpikeRaster(spikes,'LineFormat',LineFormat)
%
%       MarkerFormat - default marker is a gray dot with size 1. Used for
%           scatter type plots only. Usage is the same as LineFormat.
%
%       AutoLabel - default 0.
%           Automatically labels x-axis as 'Time (ms)' or 'Time (s)' and
%           y-axis as 'Trial'.
%
%       XLimForCell - default [NaN NaN].
%           Sets x-axis window limits if using cell spike time data. If
%           unchanged, the default limits will be 0.05% of the range. For
%           better performance, this parameter should be set.
%
%       TimePerBin - default 0.001 (1 millisecond).
%           Sets the duration of each timebin for binary spike train data.
%
%       SpikeDuration - default 0.001 (1 millisecond).
%           Sets the horizontal spike length for cell spike time data.
%
%       RelSpikeStartTime - default 0 seconds.
%           Determines the starting point of the spike relative to the time
%           indicated by spike times or time bins. For example, a relative
%           spike start time of -0.0005 would center 1ms spikes for a
%           horzline plot of binary spike data.
%
%       rasterWindowOffset - default NaN
%           Exactly the same as relSpikeStartTime, but unlike
%           relSpikeStartTime, the name implies that it can be used to make
%           x-axis start at a certain time. If set, takes precedence over
%           relSpikeStartTime.
%
%       VertSpikePosition - default 0 (centered on trial).
%           Determines where the spike position is relative to the trial. A
%           value of 0 is centered on the trial number - so a spike on
%           trial 3 would have its y-center on 3. Example: A common type of
%           spike raster plots vertical spikes from previous trial to
%           current trial. Set VertSpikePosition to -0.5 to center the
%           spike between trials.
%
%       VertSpikeHeight - default 1 (spans 1 trial).
%           Determines height of spike for 'vertline' plots. Decrease to
%           separate trials with a gap.
%
%   Examples:
%       plotSpikeRaster(spikeTimes);
%               Plots raster plot with horizontal lines.
%
%       plotSpikeRaster(spikeTimes,'PlotType','vertline');
%               Plots raster plot with vertical lines.
%
%       plotSpikeRaster(spikeTimes,'FigHandle',h,'AutoLabel',1,...
%           'XLimForCell',[0 10],'HorzSpikeLength',0.002,);
%               Plots raster plot on figure with handle h using horizontal
%               lines of length 0.002, with a window from 0 to 10 seconds,
%               and automatic labeling.
%
%       plotSpikeRaster(spikeTimes,'PlotType','scatter',...
%           'MarkerFormat',MarkerFormat);
%               Plots raster plot using dots with a format specified by
%               MarkerFormat.


%% AUTHOR    : Jeffrey Chiou
%% $DATE     : 07-Feb-2014 12:15:47 $
%% $Revision : 1.2 $
%% DEVELOPED : 8.1.0.604 (R2013a)
%% FILENAME  : plotSpikeRaster.m

%% Set Defaults and Load optional arguments
LineFormat.Color = [0.2 0.2 0.2];
MarkerFormat.MarkerSize = 1;
MarkerFormat.Color = [0.2 0.2 0.2];
MarkerFormat.LineStyle = 'none';

p = inputParser;
p.addRequired('spikes',@(x) islogical(x) || iscell(x));
p.addParamValue('FigHandle',gcf,@isinteger);
p.addParamValue('PlotType','horzLine',@ischar);
p.addParamValue('LineFormat',LineFormat,@isstruct)
p.addParamValue('MarkerFormat',MarkerFormat,@isstruct);
p.addParamValue('AutoLabel',0, @islogical);
p.addParamValue('XLimForCell',[NaN NaN],@(x) isnumeric(x) && isvector(x));
p.addParamValue('TimePerBin',0.001,@(x) isnumeric(x) && isscalar(x));
p.addParamValue('SpikeDuration',0.001,@(x) isnumeric(x) && isscalar(x));
p.addParamValue('RelSpikeStartTime',0,@(x) isnumeric(x) && isscalar(x));
p.addParamValue('RasterWindowOffset',NaN,@(x) isnumeric(x) && isscalar(x));
p.addParamValue('VertSpikePosition',0,@(x) isnumeric(x) && isscalar(x));
p.addParamValue('VertSpikeHeight',1,@(x) isnumeric(x) && isscalar(x));
p.parse(spikes,varargin{:});

spikes = p.Results.spikes;
figH = p.Results.FigHandle;
plotType = lower(p.Results.PlotType);
lineFormat = struct2opt(p.Results.LineFormat);
markerFormat = struct2opt(p.Results.MarkerFormat);
autoLabel = p.Results.AutoLabel;
xLimForCell = p.Results.XLimForCell;
timePerBin = p.Results.TimePerBin;
spikeDuration = p.Results.SpikeDuration;
relSpikeStartTime = p.Results.RelSpikeStartTime;
rasterWindowOffset = p.Results.RasterWindowOffset;
vertSpikePosition = p.Results.VertSpikePosition;
vertSpikeHeight = p.Results.VertSpikeHeight;

if ~isnan(rasterWindowOffset) && relSpikeStartTime==0
    relSpikeStartTime = rasterWindowOffset;
elseif ~isnan(rasterWindowOffset) && relSpikeStartTime~=0
    disp(['Warning: RasterWindoWOffset and RelSpikeStartTime perform the same function. '...
        'The value set in RasterWindowOffset will be used over RelSpikesStartTime']);
    relSpikeStartTime = rasterWindowOffset;
end

%% Initialize figure and begin plotting logic
figure(figH);
hold on;

if islogical(spikes)
    %% Binary spike train matrix case. Initialize variables and set axes.
    nTrials = size(spikes,1);
    nTimes = size(spikes,2);
    % Convert Parameters to correct units using TimePerBin. Default is 1 ms
    % per bin (0.001s)
    spikeDuration = spikeDuration/timePerBin;
    relSpikeStartTime = relSpikeStartTime/timePerBin;

    % Note: xlim and ylim are much, much faster than axis or set(gca,...).
    xlim([0+relSpikeStartTime nTimes+1+relSpikeStartTime]);
    ylim([0 nTrials+1]);

    switch plotType
        case 'horzline'
            %% Horizontal Lines
            % Find the trial (yPoints) and timebin (xPoints) of each spike
            [trials,timebins] = find(spikes);
            trials = trials';
            timebins = timebins';

            xPoints = [ timebins + relSpikeStartTime;
                timebins + relSpikeStartTime + spikeDuration;
                NaN(size(timebins)) ];
            yPoints = [ trials + vertSpikePosition;
                trials + vertSpikePosition;
                NaN(size(trials)) ];

            xPoints = xPoints(:);
            yPoints = yPoints(:);
            plot(xPoints,yPoints,'k',lineFormat{:});
        case 'vertline'
            %% Vertical Lines
            % Find the trial (yPoints) and timebin (xPoints) of each spike
            [trials,timebins] = find(spikes);
            trials = trials';
            timebins = timebins';
            halfSpikeHeight = vertSpikeHeight/2;

            xPoints = [ timebins + relSpikeStartTime;
                timebins + relSpikeStartTime;
                NaN(size(timebins)) ];
            yPoints = [ trials - halfSpikeHeight + vertSpikePosition;
                trials + halfSpikeHeight + vertSpikePosition;
                NaN(size(trials)) ];

            xPoints = xPoints(:);
            yPoints = yPoints(:);
            plot(xPoints,yPoints,'k',lineFormat{:});
        case 'horzline2'
            %% Horizontal lines, for many timebins
            % Plots a horizontal line the width of a time bin for each
            % spike. Efficient for fewer trials and many timebins.
            xPoints = [];
            yPoints = [];

            for trials = 1:nTrials
                % If there are spikes, plot a line. Otherwise, do nothing.
                if sum(spikes(trials,:)) > 0

                    % Find the difference in spike times. Padding at the
                    % front and back with a zero accounts for spikes in
                    % the first and last indices, and keeps startY and endY
                    % the same size
                    spikeDiff = diff([0 spikes(trials,:) 0]);
                    % Ex. [0 1 1] -> [0 0 1 1 0]. diff(...) -> [0 1 0 -1]

                    % Find line segments to plot (line segments longer than
                    % one trial are for spikes on consecutive trials)
                    startX = find(spikeDiff > 0);
                    % Ex. cont. from above: find(...) -> 2
                    endX = find(spikeDiff < 0);
                    % Ex. cont. from above: find(...) -> 4. Thus a line
                    % segment will be defined from 2 to 4 (x-axis)

                    % Combine x points and adjust x points according to
                    % parameters. Add Add NaNs to break up line segments.
                    trialXPoints = [startX + relSpikeStartTime;...
                        endX + relSpikeStartTime + (spikeDuration - 1);...
                        NaN(size(startX)) ];
                    % Explanation for 'spikeDuration - 1': spikeDuration
                    % has been converted already using timePerBin, so the
                    % offset at the end is simply duration - one timeBin,
                    % since 1 timebin's worth of spikes has passed. Only
                    % affects last spike if less than time per bin.

                    % Convert x points array to vector
                    trialXPoints = trialXPoints(:)';

                    % Add y points centered on the trial by default
                    % (adjustable with vertSpikePosition parameter)
                    trialYPoints = trials*ones(size(trialXPoints)) + vertSpikePosition;

                    % Store this trial's x and y points. Unfortunately,
                    % preallocating and trimming may actually be slower,
                    % depending on data.
                    xPoints = [xPoints trialXPoints];
                    yPoints = [yPoints trialYPoints];

                end
            end

            plot(xPoints, yPoints,'k', lineFormat{:});

        case 'vertline2'
            %% Vertical lines, for many trials
            % Plot one long line for each timebin. Reduces the number of
            % objects to plot. Efficient for many trials and fewer timebins
            xPoints = [];
            yPoints = [];

            for time = 1:nTimes
                if sum(spikes(:,time)) > 0

                    % Find the difference in spike times. Same principle as
                    % horzline2. See comments above for explanation.
                    spikeDiff = diff([0 spikes(:,time)' 0]);

                    % Find line segments to plot (line segments longer than
                    % one trial are for spikes on consecutive trials)
                    startY = find(spikeDiff > 0);
                    endY = find(spikeDiff < 0);

                    % Add NaNs to break up line segments
                    timeBinYPoints = [startY + vertSpikePosition;...
                        endY + vertSpikePosition; NaN(size(startY))];
                    % Convert to vector
                    timeBinYPoints = timeBinYPoints(:)';
                    timeBinXPoints = time*ones(size(timeBinYPoints));

                    % Store this timebin's x and y points.
                    xPoints = [xPoints timeBinXPoints];
                    % Subtract 0.5 from each y point so spikes are centered
                    % by default (adjustable with vertSpikePosition)
                    yPoints = [yPoints timeBinYPoints-0.5];

                end
            end

            plot(xPoints, yPoints, 'k', lineFormat{:});

        case 'scatter'
            %% Dots or other markers (scatterplot style)
            % Find the trial (yPoints) and timebin (xPoints) of each
            % spike
            [yPoints,xPoints] = find(spikes==1);
            xPoints = xPoints + relSpikeStartTime;
            plot(xPoints,yPoints,'.k',markerFormat{:});

        case 'imagesc'
            %% Imagesc
            imagesc(spikes);
            % Flip the colormap since the default is white for 1, black for
            % 0.
            colormap(flipud(colormap('gray')));

        otherwise
            error('Invalid plot type. Must be horzline, vertline, horzline2, vertline2, scatter, or imagesc');
    end % switch
    set(gca,'YDir','reverse');

    %% Label
    if autoLabel
        xlabel('Time (ms)');
        ylabel('Trial');
    end

else % Equivalent to elseif iscell(spikes).
    %% Cell case

    % Validation: First check to see if cell array is a vector, and each
    % trial within is a vector.
    if ~isvector(spikes)
        error('Spike cell array must be an M x 1 vector.')
    end
    trialIsVector = cellfun(@isvector,spikes);
    if sum(trialIsVector) < length(spikes)
        error('Cells must contain 1 x N vectors of spike times.');
    end

    % Now make sure cell array is M x 1 and not 1 x M.
    if size(spikes,2) > 1 && size(spikes,1) == 1
        spikes = spikes';
    end

    % Make sure each trial is 1 x N and not N x 1
    nRowsInTrial = cellfun(@(x) size(x,1),spikes);
    % If there is more than 1 row in any trial, add a warning, and
    % transpose those trials. Allows for empty trials/cells (nRows > 1
    % instead of > 0).
    if sum(nRowsInTrial > 1) > 0
        trialsToReformat = find(nRowsInTrial > 1);
        disp('Warning - some cells (trials) have more than 1 row. Those trials will be transposed.');
        for t = trialsToReformat
            spikes{trialsToReformat} = spikes{trialsToReformat}';
        end
    end

    nTrials = length(spikes);

    % Find x-axis limits that aren't manually set (default [NaN NaN]), and
    % automatically set them. This is because we don't assume spikes start
    % at 0 - we can have negative spike times.
    limitsToSet = isnan(xLimForCell);
    try
        if sum(limitsToSet) > 0
            % First find range of spike times
            minTimes = cellfun(@min,spikes,'UniformOutput',false);
            minTime = min( [ minTimes{:} ] );
            maxTimes = cellfun(@max,spikes,'UniformOutput',false);
            maxTime = max( [ maxTimes{:} ] );
            timeRange = maxTime - minTime;

            % Find 0.05% of the range.
            xStartOffset = relSpikeStartTime - 0.0005*timeRange;
            xEndOffset = relSpikeStartTime + 0.0005*timeRange + spikeDuration;
            newLim = [ minTime+xStartOffset, maxTime+xEndOffset ];
            xLimForCell(limitsToSet) = newLim(limitsToSet);
            % End result, if both limits are automatically set, is that the x
            % axis is expanded 0.1%, so you can see initial and final spikes.
        end
    end
    xlim(xLimForCell);
    ylim([0 nTrials+1]);

    if strcmpi(plotType,'vertline') || strcmpi(plotType,'horzline')
        %% Vertical or horizontal line logic
        nTotalSpikes = sum(cellfun(@length,spikes));

        % Preallocation is possible since we know how many points to
        % plot, unlike discrete case. 3 points per spike - the top pt,
        % bottom pt, and NaN.
        xPoints = NaN(nTotalSpikes*3,1);
        yPoints = xPoints;
        currentInd = 1;

        if strcmpi(plotType,'vertline')
            %% Vertical Lines
            halfSpikeHeight = vertSpikeHeight/2;
            for trials = 1:nTrials
                nSpikes = length(spikes{trials});
                nanSeparator = NaN(1,nSpikes);

                trialXPoints = [ spikes{trials} + relSpikeStartTime;...
                    spikes{trials} + relSpikeStartTime; nanSeparator ];
                trialXPoints = trialXPoints(:);

                trialYPoints = [ (trials-halfSpikeHeight)*ones(1,nSpikes);...
                    (trials+halfSpikeHeight)*ones(1,nSpikes); nanSeparator ];
                trialYPoints = trialYPoints(:);

                % Save points and update current index
                xPoints(currentInd:currentInd+nSpikes*3-1) = trialXPoints;
                yPoints(currentInd:currentInd+nSpikes*3-1) = trialYPoints;
                currentInd = currentInd + nSpikes*3;
            end

        else % (if plotType is 'horzline')
            %% Horizontal Lines
            for trials = 1:nTrials
                nSpikes = length(spikes{trials});
                nanSeparator = NaN(1,nSpikes);

                trialXPoints = [ spikes{trials} + relSpikeStartTime;...
                    spikes{trials} + relSpikeStartTime + spikeDuration;...
                    nanSeparator ];
                trialXPoints = trialXPoints(:);

                trialYPoints = [ trials*ones(1,nSpikes);...
                    trials*ones(1,nSpikes); nanSeparator ];
                trialYPoints = trialYPoints(:);

                % Save points and update current index
                xPoints(currentInd:currentInd+nSpikes*3-1) = trialXPoints;
                yPoints(currentInd:currentInd+nSpikes*3-1) = trialYPoints;
                currentInd = currentInd + nSpikes*3;
            end

        end

        % Plot everything at once! We will reverse y-axis direction later.
        plot(xPoints, yPoints, 'k', lineFormat{:});

    elseif strcmpi(plotType,'scatter')
        %% Dots or other markers (scatterplot style)
        % Combine all spike times into one vector
        xPoints = [ spikes{:} ];

        % Getting the trials is trickier. 3 steps:
        % 1. First convert all the spike times into ones.
        trials = cellfun( @(x) {ones(size(x))}, spikes );
        % 2. Then multiply by trial number.
        for trialNum = 1:length(spikes)
            trials{trialNum} = trialNum*trials{trialNum};
        end
        % 3. Finally convert into a vector
        yPoints = [ trials{:} ];

        % Now we can plot! We will reverse y-axis direction later.
        plot(xPoints,yPoints,'.k',markerFormat{:});

    elseif strcmpi(plotType,'imagesc') || strcmpi(plotType,'vertline2') || strcmpi(plotType,'horzline2')
        error('Can''t use imagesc/horzline2/vertline2 with cell array. Use with logical array of binary spike train data.');
    else
        error('Invalid plot type. With cell array of spike times, plot type must be horzline, vertline, or scatter.');
    end % plot type switching

    %% Reverse y-axis direction and label
    set(gca,'YDir','reverse');
    if autoLabel
        xlabel('Time (s)');
        ylabel('Trial');
    end

end % logical vs cell switching

%% Figure formatting
% Draw the tick marks on the outside
set(gca,'TickDir','out')

% Use special formatting if there is only a single trial.
% Source - http://labrigger.com/blog/2011/12/05/raster-plots-and-matlab/
if size(spikes,1) == 1
    set(gca,'YTick', [])                        % don't draw y-axis ticks
    set(gca,'PlotBoxAspectRatio',[1 0.05 1])    % short and wide
    set(gca,'YColor',get(gcf,'Color'))          % hide the y axis
    ylim([0.5 1.5])
end

hold off;

end % main function


function paramCell = struct2opt(paramStruct)
% Converts structure to parameter-value pairs
%   Example usage:
%       formatting = struct()
%       formatting.color = 'black';
%       formatting.fontweight = 'bold';
%       formatting.fontsize = 24;
%       formatting = struct2opt(formatting);
%       xlabel('Distance', formatting{:});
% Adapted from:
% http://stackoverflow.com/questions/15013026/how-can-i-unpack-a-matlab-structure-into-function-arguments
% by user 'yuk'

fname = fieldnames(paramStruct);
fval = struct2cell(paramStruct);
paramCell = [fname, fval]';
paramCell = paramCell(:);

end % struct2opt


function  [r,rSmo,t] = psthBins(times, binsize, fs, ntrials, triallen)
% PSTHBINS Computes the peri-stimulus time histogram from spike times.
%
% H = PSTH(TIMES, BINSIZE, FS,NTRIALS,TRIALLEN)
% TIMES - spike times (samples)
% BINSIZE - binwidth (ms)
% FS - sampling rate (hz)
% NTRIALS - number of trials
% TRIALLEN - length of a trial (samples)
%
% R - binned firing rate.
% rSmo - smoothed firing rate (look inside function for smoothing window
% T - peri stimulus time histogram edges



% Author: Rajiv Narayan
% askrajiv@gmail.com
% Boston University, Boston, MA
%
% edited: [EHS::20211108]

%Compute PSTH
lastBin = binsize * ceil((triallen)*(1000/(fs*binsize)));
edges = 0 : binsize : lastBin;
x = (mod(times,triallen))*(1000/fs);
r = (histc(x,edges)*1000) / (ntrials*binsize);
t = edges;

% and smoothing
rSmo = smoothdata(r,'gaussian',10);

end

function rgb = rgb(s)

% RGB  Rgb triple for given CSS color name
%
%   RGB = RGB('COLORNAME') returns the red-green-blue triple corresponding
%     to the color named COLORNAME by the CSS3 proposed standard [1], which
%     contains 139 different colors (an rgb triple is a 1x3 vector of
%     numbers between 0 and 1). COLORNAME is case insensitive, and for gray
%     colors both spellings (gray and grey) are allowed.
%
%   RGB CHART creates a figure window showing all the available colors with
%     their names.
%
%   EXAMPLES
%     c = rgb('DarkRed')               gives c = [0.5430 0 0]
%     c = rgb('Green')                 gives c = [0 0.5 0]
%     plot(x,y,'color',rgb('orange'))  plots an orange line through x and y
%     rgb chart                        shows all the colors
%
%   BACKGROUND
%     The color names of [1] have already been ratified in [2], and
%     according to [3] they are accepted by almost all web browsers and are
%     used in Microsoft's .net framework. All but four colors agree with
%     the X11 colornames, as detailed in [4]. Of these the most important
%     clash is green, defined as [0 0.5 0] by CSS and [0 1 0] by X11. The
%     definition of green in Matlab matches the X11 definition and gives a
%     very light green, called lime by CSS (many users of Matlab have
%     discovered this when trying to color graphs with 'g-'). Note that
%     cyan and aqua are synonyms as well as magenta and fuchsia.
%
%   ABOUT RGB
%     This program is public domain and may be distributed freely.
%     Author: Kristjn Jnasson, Dept. of Computer Science, University of
%     Iceland (jonasson@hi.is). June 2009.
%
%   REFERENCES
%     [1] "CSS Color module level 3", W3C (World Wide Web Consortium)
%         working draft 21 July 2008, http://www.w3.org/TR/css3-color
%
%     [2] "Scalable Vector Graphics (SVG) 1.1 specification", W3C
%         recommendation 14 January 2003, edited in place 30 April 2009,
%         http://www.w3.org/TR/SVG
%
%     [3] "Web colors", http://en.wikipedia.org/wiki/Web_colors
%
%     [4] "X11 color names" http://en.wikipedia.org/wiki/X11_color_names

persistent num name
if isempty(num) % First time rgb is called
    [num,name] = getcolors();
    name = lower(name);
    num = reshape(hex2dec(num), [], 3);
    % Divide most numbers by 256 for "aesthetic" reasons (green=[0 0.5 0])
    I = num < 240;  % (interpolate F0--FF linearly from 240/256 to 1.0)
    num(I) = num(I)/256;
    num(~I) = ((num(~I) - 240)/15 + 15)/16; + 240;
end
if strcmpi(s,'chart')
    showcolors()
else
    k = find(strcmpi(s, name));
    if isempty(k)
        error(['Unknown color: ' s]);
    else
        rgb = num(k(1), :);
    end
end
end

function showcolors()
[num,name] = getcolors();
grp = {'White', 'Gray', 'Red', 'Pink', 'Orange', 'Yellow', 'Brown'...
    , 'Green', 'Blue', 'Purple', 'Grey'};
J = [1,3,6,8,9,10,11];
fl = lower(grp);
nl = lower(name);
for i=1:length(grp)
    n(i) = strmatch(fl{i}, nl, 'exact');
end
clf
p = get(0,'screensize');
wh = 0.6*p(3:4);
xy0 = p(1:2)+0.5*p(3:4) - wh/2;
set(gcf,'position', [xy0 wh]);
axes('position', [0 0 1 1], 'visible', 'off');
hold on
x = 0;
N = 0;
for i=1:length(J)-1
    N = max(N, n(J(i+1)) - n(J(i)) + (J(i+1) - J(i))*1.3);
end
h = 1/N;
w = 1/(length(J)-1);
d = w/30;
for col = 1:length(J)-1;
    y = 1 - h;
    for i=J(col):J(col+1)-1
        t = text(x+w/2, y+h/10 , [grp{i} ' colors']);
        set(t, 'fontw', 'bold', 'vert','bot', 'horiz','cent', 'fontsize',10);
        y = y - h;
        for k = n(i):n(i+1)-1
            c = rgb(name{k});
            bright = (c(1)+2*c(2)+c(3))/4;
            if bright < 0.5, txtcolor = 'w'; else txtcolor = 'k'; end
            rectangle('position',[x+d,y,w-2*d,h],'facecolor',c);
            t = text(x+w/2, y+h/2, name{k}, 'color', txtcolor);
            set(t, 'vert', 'mid', 'horiz', 'cent', 'fontsize', 9);
            y = y - h;
        end
        y = y - 0.3*h;
    end
    x = x + w;
end
end

function [hex,name] = getcolors()
css = {
    %White colors
    'FF','FF','FF', 'White'
    'FF','FA','FA', 'Snow'
    'F0','FF','F0', 'Honeydew'
    'F5','FF','FA', 'MintCream'
    'F0','FF','FF', 'Azure'
    'F0','F8','FF', 'AliceBlue'
    'F8','F8','FF', 'GhostWhite'
    'F5','F5','F5', 'WhiteSmoke'
    'FF','F5','EE', 'Seashell'
    'F5','F5','DC', 'Beige'
    'FD','F5','E6', 'OldLace'
    'FF','FA','F0', 'FloralWhite'
    'FF','FF','F0', 'Ivory'
    'FA','EB','D7', 'AntiqueWhite'
    'FA','F0','E6', 'Linen'
    'FF','F0','F5', 'LavenderBlush'
    'FF','E4','E1', 'MistyRose'
    %Grey colors'
    '80','80','80', 'Gray'
    'DC','DC','DC', 'Gainsboro'
    'D3','D3','D3', 'LightGray'
    'C0','C0','C0', 'Silver'
    'A9','A9','A9', 'DarkGray'
    '69','69','69', 'DimGray'
    '77','88','99', 'LightSlateGray'
    '70','80','90', 'SlateGray'
    '2F','4F','4F', 'DarkSlateGray'
    '00','00','00', 'Black'
    %Red colors
    'FF','00','00', 'Red'
    'FF','A0','7A', 'LightSalmon'
    'FA','80','72', 'Salmon'
    'E9','96','7A', 'DarkSalmon'
    'F0','80','80', 'LightCoral'
    'CD','5C','5C', 'IndianRed'
    'DC','14','3C', 'Crimson'
    'B2','22','22', 'FireBrick'
    '8B','00','00', 'DarkRed'
    %Pink colors
    'FF','C0','CB', 'Pink'
    'FF','B6','C1', 'LightPink'
    'FF','69','B4', 'HotPink'
    'FF','14','93', 'DeepPink'
    'DB','70','93', 'PaleVioletRed'
    'C7','15','85', 'MediumVioletRed'
    %Orange colors
    'FF','A5','00', 'Orange'
    'FF','8C','00', 'DarkOrange'
    'FF','7F','50', 'Coral'
    'FF','63','47', 'Tomato'
    'FF','45','00', 'OrangeRed'
    %Yellow colors
    'FF','FF','00', 'Yellow'
    'FF','FF','E0', 'LightYellow'
    'FF','FA','CD', 'LemonChiffon'
    'FA','FA','D2', 'LightGoldenrodYellow'
    'FF','EF','D5', 'PapayaWhip'
    'FF','E4','B5', 'Moccasin'
    'FF','DA','B9', 'PeachPuff'
    'EE','E8','AA', 'PaleGoldenrod'
    'F0','E6','8C', 'Khaki'
    'BD','B7','6B', 'DarkKhaki'
    'FF','D7','00', 'Gold'
    %Brown colors
    'A5','2A','2A', 'Brown'
    'FF','F8','DC', 'Cornsilk'
    'FF','EB','CD', 'BlanchedAlmond'
    'FF','E4','C4', 'Bisque'
    'FF','DE','AD', 'NavajoWhite'
    'F5','DE','B3', 'Wheat'
    'DE','B8','87', 'BurlyWood'
    'D2','B4','8C', 'Tan'
    'BC','8F','8F', 'RosyBrown'
    'F4','A4','60', 'SandyBrown'
    'DA','A5','20', 'Goldenrod'
    'B8','86','0B', 'DarkGoldenrod'
    'CD','85','3F', 'Peru'
    'D2','69','1E', 'Chocolate'
    '8B','45','13', 'SaddleBrown'
    'A0','52','2D', 'Sienna'
    '80','00','00', 'Maroon'
    %Green colors
    '00','80','00', 'Green'
    '98','FB','98', 'PaleGreen'
    '90','EE','90', 'LightGreen'
    '9A','CD','32', 'YellowGreen'
    'AD','FF','2F', 'GreenYellow'
    '7F','FF','00', 'Chartreuse'
    '7C','FC','00', 'LawnGreen'
    '00','FF','00', 'Lime'
    '32','CD','32', 'LimeGreen'
    '00','FA','9A', 'MediumSpringGreen'
    '00','FF','7F', 'SpringGreen'
    '66','CD','AA', 'MediumAquamarine'
    '7F','FF','D4', 'Aquamarine'
    '20','B2','AA', 'LightSeaGreen'
    '3C','B3','71', 'MediumSeaGreen'
    '2E','8B','57', 'SeaGreen'
    '8F','BC','8F', 'DarkSeaGreen'
    '22','8B','22', 'ForestGreen'
    '00','64','00', 'DarkGreen'
    '6B','8E','23', 'OliveDrab'
    '80','80','00', 'Olive'
    '55','6B','2F', 'DarkOliveGreen'
    '00','80','80', 'Teal'
    %Blue colors
    '00','00','FF', 'Blue'
    'AD','D8','E6', 'LightBlue'
    'B0','E0','E6', 'PowderBlue'
    'AF','EE','EE', 'PaleTurquoise'
    '40','E0','D0', 'Turquoise'
    '48','D1','CC', 'MediumTurquoise'
    '00','CE','D1', 'DarkTurquoise'
    'E0','FF','FF', 'LightCyan'
    '00','FF','FF', 'Cyan'
    '00','FF','FF', 'Aqua'
    '00','8B','8B', 'DarkCyan'
    '5F','9E','A0', 'CadetBlue'
    'B0','C4','DE', 'LightSteelBlue'
    '46','82','B4', 'SteelBlue'
    '87','CE','FA', 'LightSkyBlue'
    '87','CE','EB', 'SkyBlue'
    '00','BF','FF', 'DeepSkyBlue'
    '1E','90','FF', 'DodgerBlue'
    '64','95','ED', 'CornflowerBlue'
    '41','69','E1', 'RoyalBlue'
    '00','00','CD', 'MediumBlue'
    '00','00','8B', 'DarkBlue'
    '00','00','80', 'Navy'
    '19','19','70', 'MidnightBlue'
    %Purple colors
    '80','00','80', 'Purple'
    'E6','E6','FA', 'Lavender'
    'D8','BF','D8', 'Thistle'
    'DD','A0','DD', 'Plum'
    'EE','82','EE', 'Violet'
    'DA','70','D6', 'Orchid'
    'FF','00','FF', 'Fuchsia'
    'FF','00','FF', 'Magenta'
    'BA','55','D3', 'MediumOrchid'
    '93','70','DB', 'MediumPurple'
    '99','66','CC', 'Amethyst'
    '8A','2B','E2', 'BlueViolet'
    '94','00','D3', 'DarkViolet'
    '99','32','CC', 'DarkOrchid'
    '8B','00','8B', 'DarkMagenta'
    '6A','5A','CD', 'SlateBlue'
    '48','3D','8B', 'DarkSlateBlue'
    '7B','68','EE', 'MediumSlateBlue'
    '4B','00','82', 'Indigo'
    %Gray repeated with spelling grey
    '80','80','80', 'Grey'
    'D3','D3','D3', 'LightGrey'
    'A9','A9','A9', 'DarkGrey'
    '69','69','69', 'DimGrey'
    '77','88','99', 'LightSlateGrey'
    '70','80','90', 'SlateGrey'
    '2F','4F','4F', 'DarkSlateGrey'
    };
hex = css(:,1:3);
name = css(:,4);
end


function [] = updateUser(operationMessage,loopVariable,period,NloopVariables)
% UPDATEUSER displays progress through a for loop in the command window.
%
%   [updateMessage] = updateUser(operationMessage,loopVariable,period,NloopVariables)
%
% operationMessage: a message about what kind of operation is being
%   performed (e.g. calculating spectrograms)
% loopVariable: the loop variable
% period: number of loop variables after which to update the user
% NloopVariables: max number of loop variables

% author EHS20170609
% updated 20181010: more intuitive input args and more efficient

if isequal(mod(loopVariable,period),0)
    fprintf('\n %s %d out of %d',operationMessage,loopVariable,NloopVariables)
end
end


function t = circ_confmean(alpha, xi, w, d, dim)
%
% t = circ_mean(alpha, xi, w, d, dim)
%   Computes the confidence limits on the mean for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [xi   (1-xi)-confidence limits are computed, default 0.05]
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%   Output:
%     t     mean +- d yields upper/lower (1-xi)% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 5
    dim = 1;
end

if nargin < 4 || isempty(d)
    % per default do not apply correct for binned data
    d = 0;
end

if nargin < 3 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% set confidence limit size to default
if nargin < 2 || isempty(xi)
    xi = 0.05;
end

% compute ingredients for conf. lim.
r = circ_r(alpha,w,d,dim);
n = sum(w,dim);
R = n.*r;
c2 = chi2inv((1-xi),1);

% check for resultant vector length and select appropriate formula
t = zeros(size(r));

for i = 1:numel(r)
    if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
        t(i) = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));  % equ. 26.24
    elseif r(i) >= .9
        t(i) = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      % equ. 26.25
    else
        t(i) = NaN;
        warning('Requirements for confidence levels not met.');
    end
end

% apply final transform
t = acos(t./R);

end


function kappa = circ_kappa(alpha,w)
%
% kappa = circ_kappa(alpha,[w])
%   Computes an approximation to the ML estimate of the concentration
%   parameter kappa of the von Mises distribution.
%
%   Input:
%     alpha   angles in radians OR alpha is length resultant
%     [w      number of incidences in case of binned angle data]
%
%   Output:
%     kappa   estimated value of kappa
%
%   References:
%     Statistical analysis of circular data, Fisher, equation p. 88
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


alpha = alpha(:);

if nargin<2
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) > size(w,1)
        w = w';
    end
end

N = length(alpha);

if N>1
    R = circ_r(alpha,w);
else
    R = alpha;
end

if R < 0.53
    kappa = 2*R + R^3 + 5*R^5/6;
elseif R>=0.53 && R<0.85
    kappa = -.4 + 1.39*R + 0.43/(1-R);
else
    kappa = 1/(R^3 - 4*R^2 + 3*R);
end

if N<15 && N>1
    if kappa < 2
        kappa = max(kappa-2*(N*kappa)^-1,0);
    else
        kappa = (N-1)^3*kappa/(N^3+N);
    end
end
end


function [mu ul ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
    dim = 1;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
    t = circ_confmean(alpha,0.05,w,[],dim);
    ul = mu + t;
    ll = mu - t;
end
end


function r = circ_r(alpha, w, d, dim)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
%
%   Output:
%     r		mean resultant length
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
    dim = 1;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

if nargin < 3 || isempty(d)
    % per default do not apply correct for binned data
    d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
    c = d/2/sin(d/2);
    r = c*r;
end
end


function [phis, cdf, phiplot, cdfplot] = circ_samplecdf(thetas, resolution)

% [phis, cdf, phiplot, cdfplot] = circ_samplecdf(thetas, resolution)
%
%   Helper function for circ_kuipertest.
%   Evaluates CDF of sample in thetas.
%
% Input:
%   thetas      sample (in radians)
%   resolution  resolution at which the cdf is evaluated
%
% Output:
%   phis        angles at which CDF is evaluated
%   cdf         CDF values at these angles
%   phiplot     as phi, for plotting
%   cdfplot     as cdf, for plotting
%
%
% Circular Statistics Toolbox for Matlab

% By Marc J. Velasco, 2009
% velasco@ccs.fau.edu

if nargin < 2
    resolution = 100;
end

phis = 0;
cdf = zeros(1, length(phis));

phis = linspace(0,2*pi,resolution+1);
phis = phis(1:end-1);

% ensure all points in thetas are on interval [0, 2pi)
x = thetas(thetas<0);
thetas(thetas<0) = (2*pi-abs(x));

% compute cdf
thetas = sort(thetas);
dprob = 1/length(thetas); %incremental change in probability
cumprob = 0; %cumultive probability so far

% for a little bit, we'll add on 2pi to the end of phis
phis = [phis 2*pi];

for j=1:resolution
    minang = phis(j);
    maxang = phis(j+1);
    currcount = sum(thetas >= minang & thetas < maxang);
    cdf(j) = cumprob + dprob*currcount;
    cumprob = cdf(j);
end

phis = phis(1:end-1);

% for each point in x, duplicate it with the preceding value in y
phis2 = phis;
cdf2 = [0 cdf(1:end-1)];

cdfplottable = [];
phisplottable = [];

for j=1:length(phis);
    phisplottable = [phisplottable phis(j) phis2(j)]; %#ok<AGROW>
    cdfplottable = [cdfplottable cdf2(j) cdf(j)]; %#ok<AGROW>
end

phiplot = [phisplottable 2*pi];
cdfplot = [cdfplottable 1];

end


function [thetahat kappa] = circ_vmpar(alpha,w,d)

% r = circ_vmpar(alpha, w, d)
%   Estimate the parameters of a von Mises distribution.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%
%   Output:
%     thetahat		preferred direction
%     kappa       concentration parameter
%
% PHB 3/23/2009
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de

alpha = alpha(:);
if nargin < 2
    w = ones(size(alpha));
end
if nargin < 3
    d = 0;
end

r = circ_r(alpha,w,d);
kappa = circ_kappa(r);

thetahat = circ_mean(alpha,w);

end


function alpha = circ_vmrnd(theta, kappa, n)

% alpha = circ_vmrnd(theta, kappa, n)
%   Simulates n random angles from a von Mises distribution, with preferred
%   direction thetahat and concentration parameter kappa.
%
%   Input:
%     [theta    preferred direction, default is 0]
%     [kappa    width, default is 1]
%     [n        number of samples, default is 10]
%
%     If n is a vector with two entries (e.g. [2 10]), the function creates
%     a matrix output with the respective dimensionality.
%
%   Output:
%     alpha     samples from von Mises distribution
%
%
%   References:
%     Statistical analysis of circular data, Fisher, sec. 3.3.6, p. 49
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens and Marc J. Velasco, 2009
% velasco@ccs.fau.edu


% default parameter
if nargin < 3
    n = 10;
end

if nargin < 2
    kappa = 1;
end

if nargin < 1
    theta = 0;
end

if numel(n) > 2
    error('n must be a scalar or two-entry vector!')
elseif numel(n) == 2
    m = n;
    n = n(1) * n(2);
end

% if kappa is small, treat as uniform distribution
if kappa < 1e-6
    alpha = 2*pi*rand(n,1);
    return
end

% other cases
a = 1 + sqrt((1+4*kappa.^2));
b = (a - sqrt(2*a))/(2*kappa);
r = (1 + b^2)/(2*b);

alpha = zeros(n,1);
for j = 1:n
    while true
        u = rand(3,1);

        z = cos(pi*u(1));
        f = (1+r*z)/(r+z);
        c = kappa*(r-f);

        if u(2) < c * (2-c) || ~(log(c)-log(u(2)) + 1 -c < 0)
            break
        end


    end

    alpha(j) = theta +  sign(u(3) - 0.5) * acos(f);
    alpha(j) = angle(exp(i*alpha(j)));
end

if exist('m','var')
    alpha = reshape(alpha,m(1),m(2));
end
end


function [mua,timestamps,thresh] = IED_MUA(data)
% does quick multiunit action potential detection with the 'RQQ' threshold.
%
%   input is a data matrix [channels x samples] of neural data sampled at
%   30 kilosamples per second.

Fs = 30000;
% number of channels
numChans=size(data,1);

% Band pass between 500 and 3000 Hz
band = [300 3000];
[b,a] = fir1(96,[band(1)/(Fs/2) band(2)/(Fs/2)]);
for ch3=1:numChans

    mua(ch3,:) = filtfilt(b,a,double(data(ch3,:)));
    display(['filtering channel ' num2str(ch3)])

    % removing the mean from each channel's signal.
    mua(ch3,:) = mua(ch3,:)-mean(mua(ch3,:));

    % finding a threshold based on the period defined in ThreshPeriod
    % thresh = 4*median(abs(mua(ch3,:)./0.6745));
    thresh = -4*rms(mua(ch3,:));
    display(['Thresholding channel ' num2str(ch3)])

    % finding the minima of the filtered signal.
    mua_peaks = find_inflections(mua(ch3,:),'minima');
    fprintf('Threshold value (uV) for channel %d: %.2f\n\n',ch3,thresh)

    % saving the times of the minima that are greater than the threshold.
    spiketimes = mua_peaks(mua(ch3,mua_peaks)<thresh);
    timestamps{ch3} = spiketimes./Fs;
end

end
