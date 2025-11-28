load('no_controls.mat');
all_test_acc1=all_test_acc;
all_test_auc1=all_test_auc;
load('with_controls.mat')

[p_acc, stats] = perm_test_paired(all_test_acc1, all_test_acc, 10000);

fprintf('Permutation test (accuracy): p = %.4g\n', p_acc);

[p, stats] = perm_test_paired(all_test_auc1, all_test_auc, 10000);

fprintf('Permutation test (accuracy): p = %.4g\n', p);

function [p, stats] = perm_test_paired(metric1, metric2, nPerm)
%PERM_TEST_PAIRED  Paired permutation test for two models on same splits.
%
%   [p, stats] = perm_test_paired(metric1, metric2, nPerm)
%
%   INPUT:
%       metric1 : [N x 1] or [1 x N] vector (e.g., accuracies of model 1)
%       metric2 : [N x 1] or [1 x N] vector (e.g., accuracies of model 2)
%       nPerm   : number of permutations (e.g., 10000)
%
%   OUTPUT:
%       p       : two-sided permutation p-value
%       stats   : struct with fields
%                   .obs_diff      - observed mean(metric1 - metric2)
%                   .null_diffs    - [nPerm x 1] permuted mean differences
%                   .null_mean     - mean(null_diffs)
%                   .null_std      - std(null_diffs)
%                   .N             - number of valid paired observations
%
%   NOTE:
%       - This is a *paired* test: assumes metric1(i) and metric2(i)
%         come from the same split / subject / fold.
%       - Under the null of no difference, the label "model1 vs model2"
%         is exchangeable within each pair; we simulate this by randomly
%         swapping the two values within each pair.
%
%   EXAMPLE:
%       [p, stats] = perm_test_paired(all_test_acc1, all_test_acc2, 10000);
%

    if nargin < 3 || isempty(nPerm)
        nPerm = 10000;
    end

    metric1 = metric1(:);
    metric2 = metric2(:);

    if numel(metric1) ~= numel(metric2)
        error('metric1 and metric2 must have the same length.');
    end

    % Remove pairs with NaNs in either metric
    valid = ~isnan(metric1) & ~isnan(metric2);
    m1 = metric1(valid);
    m2 = metric2(valid);

    N = numel(m1);
    if N < 2
        error('Not enough valid paired observations for a test.');
    end

    % Observed mean difference
    diff_obs = mean(m1 - m2);

    % Allocate null distribution
    null_diffs = zeros(nPerm,1);

    % For speed, pre-allocate temp arrays
    m1_perm = m1;
    m2_perm = m2;

    % Permutation loop
    for b = 1:nPerm
        % Randomly choose which pairs to swap
        swapMask = rand(N,1) < 0.5;

        % Swap within those pairs
        m1_perm(swapMask) = m2(swapMask);
        m2_perm(swapMask) = m1(swapMask);

        % Compute permuted mean difference
        null_diffs(b) = mean(m1_perm - m2_perm);

        % Restore original (cheap, no need to copy every time)
        m1_perm = m1;
        m2_perm = m2;
    end

    % Two-sided p-value
    p = (sum(abs(null_diffs) >= abs(diff_obs)) + 1) / (nPerm + 1);

    % Package stats
    stats = struct();
    stats.obs_diff  = diff_obs;
    stats.null_diffs = null_diffs;
    stats.null_mean = mean(null_diffs);
    stats.null_std  = std(null_diffs);
    stats.N         = N;
end
