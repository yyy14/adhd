
clear
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir,'..','PLS_improved'));
load(fullfile(script_dir,'S01_data.mat'),'T');
load(fullfile(script_dir,'S02_data.mat'),'CT','mdls_fit');

rep_total = 500;             % diagnostic first; set to 500 after inspection
binranges = 7:2.5:22;     % fixed ADHD analysis age range
K = 50;                     % fixed ADHD analysis subsample size
m = 100000;                 % Nunes original
num_subs = 70;              % Nunes original
str_md = {'linear','quadratic','cubic'};

accu = nan(rep_total,3); sensi = accu; spesi = accu;
accu_pls = accu; sensi_pls = accu; spesi_pls = accu;
nfeat_pls = nan(rep_total,3);
out_PLS = cell(rep_total,3);
split_subject_indices = cell(rep_total,4);
rng('shuffle');

for rep = 1:rep_total
    idx_a = find(T.group>0); idx_c = find(T.group==0);
    age_a = T.age(idx_a); age_c = T.age(idx_c);
    [~,~,ba] = histcounts(age_a,binranges);
    [~,~,bc] = histcounts(age_c,binranges);
    [a1,a2] = split_by_bins(ba); [c1,c2] = split_by_bins(bc);
    group = {idx_a(a1),idx_a(a2),idx_c(c1),idx_c(c2)};
    age = {age_a(a1),age_a(a2),age_c(c1),age_c(c2)};
    split_subject_indices(rep,:) = group;

    for g = 1:4
        if numel(group{g}) < K
            error('Group %d has %d subjects, less than Nunes K=%d.', ...
                g,numel(group{g}),K);
        end
    end

    % Generate candidates and retain only the 70 maximum-entropy samples.
    Idx_max = zeros(num_subs,K,4);
    for g = 1:4
        entropy = zeros(m,1); candidates = zeros(m,K,'uint16');
        for kk = 1:m
            idx = randperm(numel(age{g}),K);
            p = histcounts(age{g}(idx),[binranges inf]); p = p/sum(p);
            q = p(p>0);
            entropy(kk) = -sum(q.*log2(q));
            candidates(kk,:) = uint16(idx);
        end
        [~,order] = maxk(entropy,num_subs);
        Idx_max(:,:,g) = double(candidates(order,:));
    end

    % Dimension order follows Nunes exactly:
    % subsample x ROI x group x polynomial coefficient.
    n_roi = size(CT.map1.raw,2);
    mdls_coeff.linear = zeros(num_subs,n_roi,4,2);
    mdls_coeff.quadratic = zeros(num_subs,n_roi,4,3);
    mdls_coeff.cubic = zeros(num_subs,n_roi,4,4);
    for g = 1:4
        for s = 1:num_subs
            local = Idx_max(s,:,g);
            x = age{g}(local);
            x = x(:);
            design.linear = [x ones(K,1)];
            design.quadratic = [x.^2 x ones(K,1)];
            design.cubic = [x.^3 x.^2 x ones(K,1)];
            for md = 1:3
                y = CT.map1.(str_md{md})(group{g}(local),:);
                beta = design.(str_md{md})\y; % coefficient x ROI
                n_beta = size(beta,1);
                mdls_coeff.(str_md{md})(s,:,g,:) = ...
                    reshape(beta.',1,n_roi,1,n_beta);
            end
        end
    end

    for md = 1:3
        keep = logical(mdls_fit.map1.pVal_FDR_all.(str_md{md})(:).');
        if numel(keep) ~= n_roi
            error('%s FDR mask has %d elements, but CT has %d ROIs.', ...
                str_md{md},numel(keep),n_roi);
        end
        if ~any(keep)
            error('%s has no FDR-significant ROIs for SVM classification.', ...
                str_md{md});
        end
        coeff = mdls_coeff.(str_md{md})(:,keep,:,1); % highest order
        train = [coeff(:,:,1);coeff(:,:,3)];
        test  = [coeff(:,:,2);coeff(:,:,4)];
        labels_train = [-ones(num_subs,1);ones(num_subs,1)];
        labels_test  = labels_train;

        % Exact Nunes classifier family: MATLAB default linear SVM and no
        % Standardize. This is intentionally not LIBSVM's default RBF model.
        svm_model = fitcsvm(train,labels_train,'Standardize',false, ...
            'KernelFunction','linear');
        predicted = predict(svm_model,test);
        sensi(rep,md) = mean(predicted(1:num_subs)==-1);
        spesi(rep,md) = mean(predicted(num_subs+1:end)==1);
        accu(rep,md) = (sensi(rep,md)+spesi(rep,md))/2;

        % Fair comparison with ALL: PLS receives exactly the same unscaled
        % coefficient matrix. The only difference is the selected columns.
        train_pls = train;
        test_pls = test;
        labels_pls = [ones(num_subs,1);zeros(num_subs,1)];
        % evalc captures verbose output from PLS_FS_main and its subfunctions.
        % Runtime errors still propagate normally.
        evalc('out = PLS_FS_main(train_pls,labels_pls,0,1);');
        out_PLS{rep,md} = out;
        list_pls = unique(out.RFE.results.geneTop(:).','stable');
        list_pls = list_pls(list_pls>=1 & list_pls<=size(train_pls,2));
        nfeat_pls(rep,md) = numel(list_pls);
        if isempty(list_pls)
            warning('No PLS features selected at rep=%d, model=%s.', ...
                rep,str_md{md});
            continue
        end
        pls_model = fitcsvm(train_pls(:,list_pls),labels_train, ...
            'Standardize',false,'KernelFunction','linear');
        predicted_pls = predict(pls_model,test_pls(:,list_pls));
        sensi_pls(rep,md) = mean(predicted_pls(1:num_subs)==-1);
        spesi_pls(rep,md) = mean(predicted_pls(num_subs+1:end)==1);
        accu_pls(rep,md) = (sensi_pls(rep,md)+spesi_pls(rep,md))/2;
    end
    fprintf('rep %d/%d all: %.3f %.3f %.3f; PLS: %.3f %.3f %.3f\n', ...
        rep,rep_total,accu(rep,:),accu_pls(rep,:));
end

result_file = fullfile(script_dir,sprintf( ...
    'S03_results__rep%d.mat',rep_total));
save(result_file,'accu','sensi','spesi','accu_pls','sensi_pls', ...
    'spesi_pls','nfeat_pls','out_PLS','rep_total','binranges','K','m', ...
    'num_subs','split_subject_indices','-v7.3');
fprintf('Mean accuracy: %.4f %.4f %.4f\n',mean(accu,1,'omitnan'));
fprintf('PLS mean accuracy: %.4f %.4f %.4f\n', ...
    mean(accu_pls,1,'omitnan'));

function [half1,half2] = split_by_bins(bin_id)
half1 = []; half2 = [];
for b = 1:max(bin_id)
    members = find(bin_id==b);
    order = randperm(numel(members)); cut = round(numel(members)/2);
    half1 = [half1;members(order(1:cut))]; %#ok<AGROW>
    half2 = [half2;members(order(cut+1:end))]; %#ok<AGROW>
end
end
