function [table4, caption, note] = make_table4(matFile, outputFile)
%MAKE_TABLE4 Build Table 4 using saved summary values.


    if nargin < 1 || isempty(matFile)
        matFile = fullfile(fileparts(mfilename('fullpath')), 'table4_data.mat');
    end
    loaded = load(matFile, 'data');
    data = loaded.data;
    caption = data.caption;
    note = data.note;

    nFeatures = numel(data.featureSets);
    nMetrics = numel(data.metrics);
    globalText = strings(nMetrics, nFeatures);
    localText = strings(nMetrics, nFeatures);
    differences = strings(nMetrics, nFeatures);
    intervals = strings(nMetrics, nFeatures);
    significance = strings(nMetrics, nFeatures);

    for feature = 1:nFeatures
        for metric = 1:nMetrics
            globalStats = data.globalSummary(metric, :, feature);
            localStats = data.localSummary(metric, :, feature);
            paired = data.pairedSummary(metric, :, feature);
            globalText(metric, feature) = format_summary(globalStats);
            localText(metric, feature) = format_summary(localStats);
            differences(metric, feature) = display_minus(sprintf('%.2f', paired(1)));
            intervals(metric, feature) = display_minus( ...
                sprintf('%.2f to %.2f', paired(2), paired(3)));
            significance(metric, feature) = significance_label(paired(4));
        end
    end

    features = repelem(data.featureSets(:), nMetrics);
    metrics = repmat(data.metrics(:), nFeatures, 1);
    table4 = table(features, metrics, globalText(:), localText(:), ...
        differences(:), intervals(:), significance(:), ...
        'VariableNames', {'Feature set', 'Metric', 'Global SVM (%)', ...
        'Age-Adaptive Local SVM (%)', 'Difference (%)', ...
        '95% CI of difference (%)', 'Significance'});

    if nargin >= 2 && ~isempty(outputFile)
        writetable(table4, outputFile);
    end
    if nargout == 0
        disp(caption);
        disp(table4);
        disp(note);
    end
end

function text = format_summary(stats)
    text = string(sprintf('%.2f ± %.2f(%.2f–%.2f)', ...
        stats(1), stats(2), stats(3), stats(4)));
end

function text = display_minus(text)
    text = replace(string(text), "-", "−");
end

function label = significance_label(p)
    if p < 0.001
        label = "***";
    elseif p < 0.01
        label = "**";
    elseif p < 0.05
        label = "*";
    else
        label = "n.s.";
    end
end

