%% Rohit - 2025

function isolation_quality = compute_unit_isolation(waveform_cells, label, num_pca_components)
% COMPUTE_UNIT_ISOLATION - Optimized unit isolation calculation with Z-score normalization & PCA

fprintf('Starting unit isolation computation...\n');
start_total = tic;

[num_clusters_rows, num_clusters_cols] = size(waveform_cells);
num_clusters = num_clusters_rows * num_clusters_cols;

% Collect all waveforms into a single matrix
fprintf('Collecting all waveforms...\n');
start_collect = tic;
all_waveforms = vertcat(waveform_cells{:});
if isempty(all_waveforms)
    error('No waveforms found in input data.');
end
time_collect = toc(start_collect);
fprintf('Time taken to collect waveforms: %.2f seconds\n', time_collect);

% **Z-score normalization**
% fprintf('Applying Z-score normalization...\n');
% start_zscore = tic;
% all_waveforms = zscore(all_waveforms); % Normalize each feature to zero mean, unit variance
% time_zscore = toc(start_zscore);
% fprintf('Z-score normalization time: %.2f seconds\n', time_zscore);

% Reduce dimensionality using PCA
fprintf('Performing PCA...\n');
start_pca = tic;

% [~, score, ~, ~, ~] = pca(all_waveforms(:,11:50)); % Compute PCA

[~, score, ~, ~, explained] = pca(all_waveforms(:,11:50));
% cumulative_variance = cumsum(explained) / sum(explained); % fixed or variable number of pca components to use
% num_pca_components = find(cumulative_variance >= 0.95, 1);
all_waveforms_pca = score(:, 1:num_pca_components);

% cumulative_variance = cumsum(explained) / sum(explained);
% num_pca_components = min(5, find(cumulative_variance >= 0.95, 1)); % Retain 95% variance
% all_waveforms_pca = score(:, 1:num_pca_components); % Project data onto top components

time_pca = toc(start_pca);
fprintf('PCA computation time: %.2f seconds\n', time_pca);

% Split PCA-transformed waveforms back into clusters
fprintf('Splitting PCA-transformed waveforms back into clusters...\n');
start_split = tic;

% Compute sizes of each cluster
cluster_sizes = cellfun(@(x) size(x, 1), waveform_cells);

% Remove zero-sized clusters
valid_clusters = cluster_sizes > 0;  
filtered_cluster_sizes = cluster_sizes(valid_clusters);

% Check if all clusters are empty after filtering
if isempty(filtered_cluster_sizes)
    error('No valid clusters found. Check waveform_cells input.');
end

% Ensure filtered_cluster_sizes is a row vector and use mat2cell
split_indices = mat2cell(1:sum(filtered_cluster_sizes), 1, filtered_cluster_sizes(:)');

% cluster_sizes = cellfun(@(x) size(x, 1), waveform_cells);
% split_indices = mat2cell(1:sum(cluster_sizes), 1, cluster_sizes(:)');  % Ensure cluster_sizes is a row vector

% split_indices = mat2cell(1:sum(cluster_sizes), 1, cluster_sizes);
waveform_cells_pca = cellfun(@(idx) all_waveforms_pca(idx, :), split_indices, 'UniformOutput', false);
time_split = toc(start_split);
fprintf('Splitting time: %.2f seconds\n', time_split);

% Initialize results structure
isolation_quality(num_clusters) = struct('cluster_id', [], 'label', '', 'isolation_distance', [], 'L_ratio', []);

% valid_clusters = cluster_sizes > 0;  % Identify valid clusters
num_clusters = sum(valid_clusters);  % Ensure num_clusters reflects the count of non-empty clusters

disp(num_clusters);
num_clusters = sum(num_clusters);
disp(num_clusters);

% Compute isolation distance in parallel
fprintf('Starting parallel computation of isolation metrics...\n');
parfor idx = 1:num_clusters
    cluster_start = tic;
    
    if isempty(waveform_cells{idx})
        continue;
    end
    if ~strcmp(label{idx}, 'good')
        continue;
    end

    % Extract in-cluster waveforms in PCA space
    in_cluster_pca = waveform_cells_pca{idx};

    % Collect all out-cluster waveforms
    out_cluster_pca = vertcat(waveform_cells_pca{setdiff(1:num_clusters, idx)});

    if isempty(out_cluster_pca)
        warning('Skipping cluster %d due to empty out-cluster.', idx);
        continue;
    end

    % Compute Mahalanobis distance
    mahal_start = tic;
    mahal_dist = mahal(in_cluster_pca, out_cluster_pca);

    time_mahal = toc(mahal_start);
    
    % Compute Isolation Distance (ID)
    ID_start = tic;
    sorted_mahal = sort(mahal_dist);
    ID = sorted_mahal(min(size(in_cluster_pca, 1), numel(sorted_mahal)));

    
    time_ID = toc(ID_start);
    
    % Compute L-ratio
    L_ratio_start = tic;
    dof = size(in_cluster_pca, 2);
    L_ratio = sum(1 - chi2cdf(mahal_dist, dof)) / size(in_cluster_pca, 1);
    time_L_ratio = toc(L_ratio_start);

    % Store results
    isolation_quality(idx).cluster_id = idx;
    isolation_quality(idx).label = label{idx};
    isolation_quality(idx).isolation_distance = ID;
    isolation_quality(idx).L_ratio = L_ratio;

    cluster_time = toc(cluster_start);
    fprintf('Cluster %d: Mahal Time: %.2fs, ID Time: %.2fs, L-ratio Time: %.2fs, Total: %.2fs\n', ...
        idx, time_mahal, time_ID, time_L_ratio, cluster_time);
end

total_time = toc(start_total);
fprintf('Total execution time: %.2f seconds\n', total_time);

end
