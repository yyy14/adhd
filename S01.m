
clear
script_dir = fileparts(mfilename('fullpath'));
load(fullfile(script_dir, '..', 'T.mat'), 'T');
rng(9, 'twister');


T.siteID = T.Site;
T.sex    = T.Gender;
T.age    = T.Age;
T.group  = T.DX;


keep =  T.age >= 7 & T.age <= 22 & ...
       ~isnan(T.siteID) & ~isnan(T.group);
T = T(keep, :);

thickness = nan(height(T), 68);
for i = 1:height(T)
    thickness(i,:) = reshape(T.data_68{i}(4,:), 1, 68);
end
valid_ct = all(isfinite(thickness), 2);
T = T(valid_ct,:);
thickness = thickness(valid_ct,:);


centers = unique(T.siteID);
counts = zeros(numel(centers),2);
for c = 1:numel(centers)
    at_site = T.siteID == centers(c);
    counts(c,:) = [nnz(at_site & T.group > 0), nnz(at_site & T.group == 0)];
end
drop_centers = centers(any(counts < 10,2));
drop = ismember(T.siteID, drop_centers);
T(drop,:) = [];
thickness(drop,:) = [];


out_count = zeros(height(T),1);
centers = unique(T.siteID);
for c = 1:numel(centers)
    at_site = T.siteID == centers(c);
    z = zscore(thickness(at_site,:));
    z(~isfinite(z)) = 0;
    out_count(at_site) = sum(abs(z) > 3.5, 2);
end
drop = out_count > 1;
T(drop,:) = [];
thickness(drop,:) = [];


centers = unique(T.siteID);
counts = zeros(numel(centers),2);
for c = 1:numel(centers)
    at_site = T.siteID == centers(c);
    counts(c,:) = [nnz(at_site & T.group > 0), nnz(at_site & T.group == 0)];
end
drop_centers = centers(any(counts < 10,2));
drop = ismember(T.siteID, drop_centers);
T(drop,:) = [];
thickness(drop,:) = [];


centers = unique(T.siteID);
remove_idx = [];
min_tol = 0.13;
for c = 1:numel(centers)
    at_site = T.siteID == centers(c);
    adhd = find(at_site & T.group > 0);
    tdc  = find(at_site & T.group == 0);
    ratio = numel(adhd) / numel(tdc);
    if ratio > 1 + min_tol
        n_remove = round(numel(adhd) - numel(tdc)*(1 + min_tol));
        order = randperm(numel(adhd));
        remove_idx = [remove_idx; adhd(order(1:n_remove))]; %#ok<AGROW>
    elseif ratio < 1 - min_tol
        n_remove = ceil(numel(tdc) - numel(adhd)*(1 + min_tol/2));
        order = randperm(numel(tdc));
        remove_idx = [remove_idx; tdc(order(1:n_remove))]; %#ok<AGROW>
    end
end
remove_idx = unique(remove_idx);
T(remove_idx,:) = [];
thickness(remove_idx,:) = [];


centers = unique(T.siteID);
age_remove = [];
age_edges = 7:0.5:22;
for c = 1:numel(centers)
    at_site = T.siteID == centers(c);
    ia = find(at_site & T.group > 0);
    ic = find(at_site & T.group == 0);
    [h,~] = ttest2(T.age(ia), T.age(ic));
    if ~h, continue; end
    ha = histcounts(T.age(ia), [age_edges inf]);
    hc = histcounts(T.age(ic), [age_edges inf]);
    d = ha/sum(ha) - hc/sum(hc);
    for b = find(abs(d) > 0.049)
        lo = age_edges(b); hi = lo + 0.5;
        if d(b) > 0
            candidates = ia(T.age(ia) >= lo & T.age(ia) < hi);
            n_remove = round(d(b)*numel(ia)*(1-min_tol));
        else
            candidates = ic(T.age(ic) >= lo & T.age(ic) < hi);
            n_remove = round(abs(d(b))*numel(ic)*(1-min_tol));
        end
        n_remove = min(n_remove, numel(candidates));
        age_remove = [age_remove; candidates(1:n_remove)]; %#ok<AGROW>
    end
end
age_remove = unique(age_remove);
T(age_remove,:) = [];
thickness(age_remove,:) = [];

CT.map1.raw = thickness;
save(fullfile(script_dir,'S01_data.mat'),'T','CT','-v7.3');


