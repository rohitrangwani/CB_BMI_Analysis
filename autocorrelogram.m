%% Rohit 2025

function autocorrelogram(spikeTimestamps, binSize, maxLag)
    % Generates and plots an autocorrelogram from spike timestamps.
    %
    % Parameters:
    % spikeTimestamps: Array of spike timestamps (in milliseconds or seconds).
    % binSize: Bin size for the histogram (same units as spike_times).
    % maxLag: Maximum lag for the autocorrelogram (same units as spike_times).


    % Compute histogram edges and initialize counts
edges = -maxLag:binSize:maxLag; % Bin edges
autocorrCounts = zeros(1, length(edges) - 1); % Preallocate histogram counts

% Process each spike incrementally
nSpikes = length(spikeTimestamps);
for i = 1:nSpikes
    % Compute time differences relative to the current spike
    diffs = spikeTimestamps - spikeTimestamps(i);
    
    % Exclude self-comparison (time lag = 0)
    diffs(i) = NaN;
    
    % Restrict to maximum lag range
    validDiffs = diffs(abs(diffs) <= maxLag);
    
    % Accumulate histogram counts
    autocorrCounts = autocorrCounts + histcounts(validDiffs, edges);
end

% Normalize by the number of spikes and bin size
autocorrCounts = autocorrCounts / (nSpikes * binSize);

% Plot the autocorrelogram
binCenters = edges(1:end-1) + binSize / 2; % Compute bin centers
% figure;
bar(binCenters, autocorrCounts, 'histc');
xlabel('Lag (s)');
ylabel('Normalized Count');
title('Spike Autocorrelogram');
xlim([-maxLag, maxLag]);

end
