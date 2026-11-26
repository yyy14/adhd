
clear
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir,'functions'));
load(fullfile(script_dir,'S01_data.mat'),'T','CT');


str_md = {'linear','quadratic','cubic'};
str_at = {'map1'};

% One site is the reference category, as in an intercept-containing model.
centers = unique(T.siteID);
reference_center = centers(1);
dummy_centers = centers(2:end);
dummy_site = zeros(height(T),numel(dummy_centers));
for c = 1:numel(dummy_centers)
    dummy_site(:,c) = T.siteID == dummy_centers(c);
end

a0 = ones(size(T.age)); a1 = T.age; a2 = a1.^2; a3 = a1.^3;
g1 = T.group > 0;
mdl_site.linear    = [a1 a1.*g1 a0 g1 dummy_site];
mdl_site.quadratic = [a2 a2.*g1 a1 a1.*g1 a0 g1 dummy_site];
mdl_site.cubic     = [a3 a3.*g1 a2 a2.*g1 a1 a1.*g1 a0 g1 dummy_site];

for at = 1:numel(str_at)
    raw = CT.(str_at{at}).raw;
    for md = 1:numel(str_md)

            design = mdl_site.(str_md{md});

        for s = 1:size(raw,2)
            coeff = regress(raw(:,s),design);
            site_coeff = coeff(end-size(dummy_site,2)+1:end);
            CT.(str_at{at}).(str_md{md})(:,s) = raw(:,s)-dummy_site*site_coeff;
        end
    end
end

str_grp = {'adhd','tdc'};
group = {find(T.group>0),find(T.group==0)};
ages = {T.age(group{1}),T.age(group{2})};
forms = {'y ~ x1','y ~ x1^2','y ~ x1^3'};
mdls_fit = struct;
for at = 1:numel(str_at)
    for md = 1:numel(str_md)
        ct = CT.(str_at{at}).(str_md{md});
        for s = 1:size(ct,2)
            for g = 1:2
                gm = fitglm(ages{g},ct(group{g},s),forms{md});
                mdls_fit.(str_at{at}).(str_grp{g}).pVal(s,md) = ...
                    gm.devianceTest.pValue(2);
                mdls_fit.(str_at{at}).(str_grp{g}).AIC(s,md) = ...
                    gm.ModelCriterion.AIC;
            end
        end
    end
end

for md = 1:numel(str_md)
    for g = 1:2
        pvals = mdls_fit.map1.(str_grp{g}).pVal(:,md);
        mdls_fit.map1.(str_grp{g}).pVal_FDR.(str_md{md}) = ...
            fdr_bh(pvals,0.05,'dep','yes');
    end
    mdls_fit.map1.pVal_FDR_all.(str_md{md}) = ...
        mdls_fit.map1.adhd.pVal_FDR.(str_md{md}) & ...
        mdls_fit.map1.tdc.pVal_FDR.(str_md{md});
end

save(fullfile(script_dir,'S02_data.mat'), ...
    'mdls_fit','CT','center_mode','reference_center','-v7.3');

