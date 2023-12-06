function unwrappedCell = LK_unwrapNestedCell(nestedCell)
%
% LK_unwrapNestedCell unwraps a nested cell array into an unnested cell.
%
% Lukas Kunz, 2021

unwrappedCell   = cell(sum(cellfun(@(x) numel(x), nestedCell)), 1);
idx             = 1;
for i = 1:size(nestedCell, 1) % loop through cells
    for j = 1:numel(nestedCell{i}) % loop through cells nested in cell i
        unwrappedCell{idx, 1}   = nestedCell{i, 1}{j};
        idx                     = idx + 1;
    end
end