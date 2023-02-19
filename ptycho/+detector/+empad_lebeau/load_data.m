

function [ p ] = load_data( p )
det = p.detectors(p.scanID);
detStorage = det.detStorage;
if isempty(detStorage.files)
    [p] = det.params.get_filename(p);
end

files = detStorage.files;
file = files{p.scanID};

[~, filename, ~] = fileparts(file);
raw_parts = strsplit(filename, '_');
x_part = raw_parts(startsWith(raw_parts, "x"));
y_part = raw_parts(startsWith(raw_parts, "y"));
nx = str2double(x_part{1}(2:end));
ny = str2double(y_part{1}(2:end));

utils.verbose(2, "Opening '%s'", file);
if exist(file, 'file') ~= 2
    error("Raw data file '%s' not found.", file);
end
data = single(fread(fopen(file,'r'), 128*130*nx*ny, 'float32'));

% resize as 4D
data = reshape(data, 128, 130, ny, nx);
% crop junk rows, flip reciprocal y axis
data = data(:, 128:-1:1, :, :);
% transpose reciprocal space
data = permute(data, [2, 1, 3, 4]);

if isfield(p.detector, 'sim') && p.detector.sim
    % scale simulated data by beam current
    if isfield(p.detector, 'beam_current')
        current = p.detector.beam_current;
    else
        current = 30.;
    end
    data = data .* (current * 6241.51); % I*1e-12 C/s / (1.602e-19 C/elec) * 1 ms
else
    % scale experimental data by single-electron intensity
    data = data ./ 375;
end

if isfield(p.detector, 'crop') && ~isempty(p.detector.crop)
    crop = num2cell(p.detector.crop);
    [min_x, max_x, min_y, max_y] = crop{:};
    utils.verbose(2, "Cropping data to %d:%d x %d:%d", min_x, max_x, min_y, max_y)
    data = data(:, :, colon(min_x, max_x), colon(min_y, max_y));
    ny = size(data, 3);
    nx = size(data, 4);
end

% flip scan
data = data(:, :, ny:-1:1, nx:-1:1);

if isfield(p, 'src_positions') && p.src_positions == "matlab_pos"
    % check raster scan dimensions
    if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster" && (ny ~= p.scan.ny || nx ~= p.scan.nx)
        utils.verbose(1, "Warning: Specified raster scan size %dx%d does not match actual scan size %dx%d", p.scan.nx, p.scan.ny, nx, ny);
    end
end

if isfield(p.detector, 'poisson') && p.detector.poisson
    % TODO do this after bg subtract?
    utils.verbose(2, "Applying poisson noise.");
    data = single(poissrnd(max(0., data)));
end

if isfield(p.detector, 'psf_sigma') && p.detector.psf_sigma
    utils.verbose(2, "Applying gaussian PSF, sigma %f", p.detector.psf_sigma);
    sigma = p.detector.psf_sigma;
    ny = size(data, 1);
    nx = size(data, 2);
    ys = -floor(ny/2):ceil(ny/2)-1;
    xs = -floor(nx/2):ceil(nx/2)-1;

    gaussian = exp(-(ys'.^2 + xs.^2)/(2*sigma^2)) / (sigma*sqrt(2*pi));

    data = abs(fftshift(ifft2(fft2(fftshift(data)) .* fft2(fftshift(gaussian)))));
end

% if isfield(p, 'upsampling') && p.upsampling > 0
%     factor = 2^p.upsampling;
%     utils.verbose(2, "Upsampling data by factor of %d", factor);
%     data = utils.unbinning_2D(data, factor) / factor^2;
%     p.asize = p.asize .* factor;
% end

utils.verbose(1, "Loaded data from: '%s'", file);

detStorage.data = data;
end

