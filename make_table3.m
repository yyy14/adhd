function [table3, caption, rankings] = make_table3(matFile, outputFile)
%MAKE_TABLE3 Build Table 3 from  table3_data.mat.
%   T = MAKE_TABLE3() uses the MAT file next to this function.
%   T = MAKE_TABLE3([], 'Table3.xlsx') also exports the significant-ROI table.
%   [T, CAPTION,  RANKINGS] = MAKE_TABLE3() returns the full lists too.

    if nargin < 1 || isempty(matFile)
        matFile = fullfile(fileparts(mfilename('fullpath')), 'table3_data.mat');
    end
    loaded = load(matFile, 'data');
    data = loaded.data;
    caption = data.caption;
    
    rankings = data.rankings;
    significant = cellfun(@(t) t(t.FDR < data.fdrThreshold, :), ...
        rankings, 'UniformOutput', false);

    rowCount = max(cellfun(@height, significant));
    columns = strings(rowCount, 2 * numel(significant));
    headers = strings(1, size(columns, 2));
    for k = 1:numel(significant)
        selected = significant{k};
        n = height(selected);
        formatted = compose("%.3f", selected.FDR);
        formatted(selected.FDR < 0.001) = "<0.001";
        columns(1:n, 2*k-1) = formatted;
        columns(1:n, 2*k) = selected.ROI;
        headers(2*k-1:2*k) = data.methods(k) + [" FDR", " ROI"];
    end

    table3 = [table((1:rowCount)', 'VariableNames', {'Rank'}), ...
        array2table(columns, 'VariableNames', cellstr(headers))];
    if nargin >= 2 && ~isempty(outputFile)
        writetable(table3, outputFile);
    end
    if nargout == 0
        disp(caption);
        disp(table3);
        
    end
end
