%% Frequency and final-FDR analysis of 500-run  selections
clear
script_dir = fileparts(mfilename('fullpath'));
result_file = fullfile(script_dir,'S03_results_rep500.mat');
load(result_file,'out_PLS','nfeat_pls');
load(fullfile(script_dir,'S02_data.mat'),'mdls_fit');
roi_names_file = 'name.mat';
load(roi_names_file,'name');

q_threshold = 0.05;
model_names = {'linear','quadratic','cubic'};
results = struct;
json_models = cell(1,3);

for md = 1:3
    model_name = model_names{md};
    feature_to_roi_idx = find(mdls_fit.map1.pVal_FDR_all.(model_name));
    pfinal = numel(feature_to_roi_idx);
    n_total = size(out_PLS,1);
    selected_sets = cell(n_total,1);
    n_selected = nan(n_total,1);
    valid = false(n_total,1);

    for r = 1:n_total
        if isempty(out_PLS{r,md}), continue; end
        x = out_PLS{r,md};
        if ~isfield(x,'RFE') || ~isfield(x.RFE,'results') || ...
                ~isfield(x.RFE.results,'geneTop'), continue; end
        idx = unique(x.RFE.results.geneTop(:).');
        idx = idx(idx>=1 & idx<=pfinal);
        selected_sets{r} = idx;
        n_selected(r) = numel(idx);
        valid(r) = true;
    end

    selected_sets = selected_sets(valid);
    n_selected = n_selected(valid);
    n_rep = numel(selected_sets);
    counts = zeros(pfinal,1);
    for r = 1:n_rep
        counts(selected_sets{r}) = counts(selected_sets{r}) + 1;
    end
    frequency = counts/n_rep;
    avg_n = mean(n_selected);
    median_n = median(n_selected);
    sd_n = std(n_selected);
    pr = avg_n/pfinal;

    fdr_by_k = nan(n_rep,1);
    expected_by_k = nan(n_rep,1);
    observed_by_k = nan(n_rep,1);
    for k = 1:n_rep
        observed = sum(counts>=k);
        expected = pfinal*(1-binocdf(k-1,n_rep,pr));
        observed_by_k(k) = observed;
        expected_by_k(k) = expected;
        if observed>0, fdr_by_k(k) = expected/observed; end
    end
    feature_fdr = nan(pfinal,1);
    for j = 1:pfinal
        if counts(j)>0, feature_fdr(j)=fdr_by_k(counts(j)); end
    end
    significant = feature_fdr<q_threshold;

    [~,order] = sortrows([-counts feature_fdr],[1 2]);
    reduced_idx = (1:pfinal)';
    original_idx = feature_to_roi_idx(:);
    roi_name = string(name(original_idx));
    frequency_table = table(reduced_idx,original_idx,roi_name,counts,frequency, ...
        feature_fdr,significant,'VariableNames',{'ReducedFeatureIndex', ...
        'OriginalROIIndex','ROIName','SelectionCount','SelectionFrequency', ...
        'FinalFDR','FDRSignificant'});
    frequency_table = frequency_table(order,:);
    significant_table = frequency_table(frequency_table.FDRSignificant,:);

    first_k = find(fdr_by_k<q_threshold & observed_by_k>0,1,'first');
    if isempty(first_k), first_k=NaN; cutoff_frequency=NaN;
    else, cutoff_frequency=first_k/n_rep; end

    result = struct;
    result.model = model_name;
    result.nRep = n_rep;
    result.pfinal = pfinal;
    result.average_selected = avg_n;
    result.median_selected = median_n;
    result.sd_selected = sd_n;
    result.null_selection_probability = pr;
    result.q_threshold = q_threshold;
    result.minimum_significant_count = first_k;
    result.minimum_significant_frequency = cutoff_frequency;
    result.n_FDR_significant = height(significant_table);
    result.frequency_table = frequency_table;
    result.significant_table = significant_table;
    result.fdr_by_k = fdr_by_k;
    result.expected_by_k = expected_by_k;
    result.observed_by_k = observed_by_k;
    result.n_selected_each_rep = n_selected;
    results.(model_name) = result;

    json_models{md} = struct('model',model_name,'nRep',n_rep,'pfinal',pfinal, ...
        'averageSelected',avg_n,'medianSelected',median_n,'sdSelected',sd_n, ...
        'nullSelectionProbability',pr,'qThreshold',q_threshold, ...
        'minimumSignificantCount',first_k, ...
        'minimumSignificantFrequency',cutoff_frequency, ...
        'nFDRSignificant',height(significant_table), ...
        'frequencyRows',table2struct(frequency_table), ...
        'significantRows',table2struct(significant_table));

    fprintf('%s: mean %.3f, median %.3f, SD %.3f, FDR significant %d/%d\n', ...
        model_name,avg_n,median_n,sd_n,height(significant_table),pfinal);
end

save(fullfile(script_dir,'final_FDR_PLS_selected_ROIs_rep500.mat'), ...
    'results','q_threshold','model_names','result_file','roi_names_file','-v7.3');
json_text = jsonencode(struct('sourceResult',result_file, ...
    'roiNamesSource',roi_names_file,'models',{json_models}));
fid = fopen(fullfile(script_dir,'PLS_feature_frequency_FDR_rep500.json'),'w');
cleanup = onCleanup(@() fclose(fid));
fwrite(fid,json_text,'char');
