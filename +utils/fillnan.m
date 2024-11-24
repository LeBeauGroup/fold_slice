
function [arr] = fixnan(arr, skipdims)
%FIXNAN Replaces NaN values by averaging adjacent values in the array
% Values are replaced in the first 2 dimensions, after skipping `skipdims`.

if anynan(arr)
    nanind = find(any(isnan(arr), 1:skipdims));
    sz = size(arr);
    sz = sz(skipdims+1:end);

    [row, col, rest] = ind2sub(sz, nanind);

    inds = cell(1, skipdims + 3);
    for i = 1:skipdims
        inds{i} = ':';
    end

    % this is probably very slow. don't look.
    for i = 1:size(row)
        sum = 0;
        count = 0;

        inds{skipdims + 1} = row(i);
        inds{skipdims + 2} = col(i);
        inds{skipdims + 3} = rest(i);

        for rowoffset = [-1, 1]
            if row(i) + rowoffset >= 1 && row(i) + rowoffset <= sz(1)
                inds{skipdims + 1} = row(i) + rowoffset;
                val = arr(inds{:});
                if ~isnan(val)
                    sum = sum + val;
                    count = count + 1;
                end
            end
        end
        inds{skipdims + 1} = row(i);

        for coloffset = [-1, 1]
            if col(i) + coloffset >= 1 && col(i) + coloffset <= sz(2)
                inds{skipdims + 2} = col(i) + coloffset;
                val = arr(inds{:});
                if ~isnan(val)
                    sum = sum + val;
                    count = count + 1;
                end
            end
        end
        inds{skipdims + 2} = col(i);

        arr(inds{:}) = sum / count;
    end

    %fprintf("Replaced %d NaN frame(s) in data\n", size(row, 1));
    utils.verbose(2, "Replaced %d NaN frame(s) in data", size(row, 1));
end
end
