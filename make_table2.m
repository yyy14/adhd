function [table2, caption, note] = make_table2(matFile, outputFile)
%MAKE_TABLE2 Build Table 2 using saved values.
%   T = MAKE_TABLE2() loads the MAT file beside this function.


    if nargin < 1 || isempty(matFile)
        matFile = fullfile(fileparts(mfilename('fullpath')), 'table2_data.mat');
    end
    loaded = load(matFile, 'data');
    data = loaded.data;
    caption = data.caption;
    note = data.note;
    values = data.summaryValues;

    formatted = reshape(compose("%.2f", values(:)), size(values));
    formatted = regexprep(formatted, "0+$", "");
    formatted = regexprep(formatted, "\.$", "");
    formatted(isnan(values)) = "-";

    methods = repelem(data.methods(:), numel(data.models));
    models = repmat(data.models(:), numel(data.methods), 1);
    table2 = [table(methods, models, 'VariableNames', {'Method', 'Model'}), ...
        array2table(formatted, 'VariableNames', cellstr(data.columnNames))];

    if nargin >= 2 && ~isempty(outputFile)
        writetable(table2, outputFile);
    end
    if nargout == 0
        disp(caption);
        disp(table2);
        disp(note);
    end
end


