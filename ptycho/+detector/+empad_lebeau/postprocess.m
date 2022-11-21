%POSTPROCESS postprocess functions
%   applied to the 4D raw data
%
% ** p      p structure
% returns:
% ++ p      updated p structure
%
% see also: detector.prep_data.matlab_ps.matlab_ps


function [ p ] = postprocess( p )

detStorage = p.detectors(p.scanID).detStorage;

%% apply some custom correction such as background subtraction on detStorage.data

if isfield(p.detector, 'bg_sub') && p.detector.bg_sub
    if islogical(p.detector.bg_sub)
        % subtract anything <0.2% of peak intensity
        factor = 0.2e-2;
    else
        factor = p.detector.bg_sub;
    end
    bg = factor * quantile(detStorage.data, 0.99, 'all');
    utils.verbose(2, "Subtracting background intensity (%d)...", bg);

    detStorage.data = max(detStorage.data - bg, 0.);
end

