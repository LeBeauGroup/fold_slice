

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
% NOTE: 'nx' and 'ny' are in fold_slice orientation,
% which is transposed from output orientation
ny = str2double(x_part{1}(2:end));
nx = str2double(y_part{1}(2:end));

utils.verbose(2, "Opening '%s'", file);
if exist(file, 'file') ~= 2
    error("Raw data file '%s' not found.", file);
end
data = single(fread(fopen(file,'r'), 128*130*nx*ny, 'float32'));

% resize as 4D
data = reshape(data, 128, 130, nx, ny);
% crop junk rows, flip reciprocal y axis
data = data(:, 128:-1:1, :, :);
% transpose reciprocal space
data = permute(data, [2, 1, 3, 4]);
size(data)

if isfield(p.detector, 'crop') && ~isempty(p.detector.crop)
    crop = num2cell(p.detector.crop);
    [min_x, max_x, min_y, max_y] = crop{:};
    utils.verbose(2, "Cropping data to %d:%d x %d:%d", min_x, max_x, min_y, max_y)
    % crop in output orientation, not fold_slice orientation
    data = data(:, :, colon(min_y, max_y), colon(min_x, max_x));
    nx = size(data, 3);
    ny = size(data, 4);

    if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster"
        positions = reshape(p.positions, p.scan.nx, p.scan.ny, 2);
        positions = positions(min_y:max_y, min_x:max_x, :);
        p.positions = reshape(positions, nx*ny, 2);
        p.numpts = nx*ny;
        p.scanidxs{1, 1} = 1:nx*ny;

        p.scan.nx = nx;
        p.scan.ny = ny;
    end
end

if isfield(p.detector, 'step') && ~isempty(p.detector.step)
    % should broadcast
    step = zeros(1, 2);
    step(:) = p.detector.step;
    utils.verbose(1, "Using every %dx%d scan position", step(1), step(2)); 

    data = data(:, :, 1:step(1):size(data, 3), 1:step(2):size(data, 4));
    nx = size(data, 3);
    ny = size(data, 4);
    utils.verbose(1, "New size: %dx%d", nx, ny); 

    utils.verbose(1, "Data shape: %dx%dx%dx%d", size(data, 1), size(data, 2), size(data, 3), size(data, 4));

    if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster"
        positions = reshape(p.positions, p.scan.nx, p.scan.ny, 2);
        positions = positions(1:step(1):p.scan.nx, 1:step(2):p.scan.ny, :);
        p.positions = reshape(positions, nx*ny, 2);
        p.numpts = nx*ny;
        p.scanidxs{1, 1} = 1:nx*ny;

        p.scan.nx = nx;
        p.scan.ny = ny;
        p.scan.step_size_x = p.scan.step_size_x * step(1);
        p.scan.step_size_y = p.scan.step_size_y * step(2);
    end
end

if isfield(p.detector, 'fill_nan') && p.detector.fill_nan
    data = utils.fillnan(data, 2);
end

if isfield(p.detector, 'sim') && p.detector.sim
    if isfield(p.detector, 'beam_dose') && ~isempty(p.detector.beam_dose)
        utils.verbose(2, "Scaling data to dose: %.1f e/A^2", p.detector.beam_dose);
        if ~(isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster")
            error("Option 'beam_dose' requires a raster scan");
        end
        area = p.scan.step_size_x * p.scan.step_size_y; % A^2
        scale = p.detector.beam_dose * area;
    else
        if isfield(p.detector, 'beam_current')
            current = p.detector.beam_current;
        else
            current = 30.;
        end
        utils.verbose(2, "Scaling data to beam current: %.1f pA", current);
        scale = current * 6241.51; % I*1e-12 C/s / (1.602e-19 C/elec) * 1 ms
    end
    % scale simulated data by beam current or dose
    utils.verbose(3, "Scaling by: %.1f", scale);
    data = data .* scale;
else
    % scale experimental data by single-electron intensity
    % TODO grab from metadata
    data = data ./ 375;
end

% flip scan
data = data(:, :, nx:-1:1, ny:-1:1);

if isfield(p, 'src_positions') && p.src_positions == "matlab_pos"
    % check raster scan dimensions
    if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster" && (ny ~= p.scan.ny || nx ~= p.scan.nx)
        utils.verbose(1, "Warning: Specified raster scan size %dx%d does not match actual scan size %dx%d", p.scan.nx, p.scan.ny, nx, ny);
    end
end

if isfield(p, 'tile') && ~isempty(p.tile) && ~all(p.tile == 1)
    % tile patterns (useful for simulation)
    tile = num2cell(p.tile);
    [tile_x, tile_y] = tile{:};
    utils.verbose(1, "Tiling data '%dx%d' in realspace", tile_x, tile_y);

    % tile in output orientation
    data = repmat(data, 1, 1, tile_y, tile_x);

    if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster"
        p.scan.nx = p.scan.nx * tile_y;
        p.scan.ny = p.scan.ny * tile_x;
    end
    p.scanidxs{1, 1} = 1:nx*ny;
end

if isfield(p.detector, 'poisson') && ~isempty(p.detector.poisson) && p.detector.poisson
    % TODO do this after bg subtract?
    utils.verbose(2, "Applying poisson noise.");
    data = single(poissrnd(max(0., data)));
end

if isfield(p.detector, 'psf_sigma') && ~isempty(p.detector.psf_sigma) && p.detector.psf_sigma > 0
    utils.verbose(2, "Applying gaussian PSF, sigma %f", p.detector.psf_sigma);
    sigma = p.detector.psf_sigma;
    ny = size(data, 1);
    nx = size(data, 2);
    ys = -floor(ny/2):ceil(ny/2)-1;
    xs = -floor(nx/2):ceil(nx/2)-1;

    gaussian = exp(-(ys'.^2 + xs.^2)/(2*sigma^2));
    gaussian = gaussian / sum(gaussian, 'all');

    data = real(fftshift(ifft2(fft2(fftshift(data)) .* fft2(fftshift(gaussian)))));
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

