%BINNING postprocess functions
% apply binning or upsampling on the measured data 
%
% ** p      p structure
% returns:
% ++ p      updated p structure
%
% see also: detector.prep_data.matlab_ps.matlab_aps


function [ p ] = binning( p )


data = p.detectors(p.scanID).detStorage.data;
fmask = p.detectors(p.scanID).detStorage.fmask;

%% pad measured data

data_size = size(data);
data_size = data_size(1:2);
if any(p.asize > data_size)
    utils.verbose(2, "Padding reciprocal space to %dx%d", p.asize(1), p.asize(2));
    %shift = (p.asize - data_size) / 2;
    data = utils.imshift_fast(data, 0, 0, p.asize, 'nearest');
    fmask = utils.imshift_fast(fmask, 0, 0, p.asize, 'nearest', 1);
end

%% apply binning or upsampling on the measured data 

if isfield(p.detector,'binning')&& p.detector.binning
    binning = 2^p.detector.binning;
    data = utils.binning_2D(data, binning) * binning^2; 
    fmask =utils.binning_2D(fmask,binning) == 1;  % remove all binned pixel where at least one was masked
end

        
if isfield(p.detector,'upsampling')&& p.detector.upsampling

    upsample = 2^p.detector.upsampling;
    utils.verbose(2, 'Upsample diffraction patterns by %d',upsample)

    if isfield(p.detector,'upsampling_method') && (strcmp(p.detector.upsampling_method,'bilinear') || strcmp(p.detector.upsampling_method,'bicubic'))
        %added by YJ
        %disp(p.detector.upsampling_method)
        utils.verbose(2, 'Upsample method: imresize with %s interpolation',p.detector.upsampling_method)
        data = imresize(data, upsample, p.detector.upsampling_method) / upsample^2;
    else
        utils.verbose(2, 'Upsample method: utils.unbinning_2D')
        data = utils.unbinning_2D(data, upsample) / upsample^2; %PSI method
    end
    fmask = utils.unbinning_2D(fmask,upsample) == 1;  % remove all binned pixel where at least one was masked       
end

p.detectors(p.scanID).detStorage.data = data;
p.detectors(p.scanID).detStorage.fmask = fmask;

end

