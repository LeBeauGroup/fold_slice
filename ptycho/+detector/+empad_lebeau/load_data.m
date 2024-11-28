

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
data = reshape(data, 128, 130, nx, ny);
% crop junk rows, flip reciprocal y axis
data = data(:, 128:-1:1, :, :);
% transpose reciprocal space (to column-major)
% we keep the scan row major because that's how the raster is created
% k_y, k_x, s_x, s_y
data = permute(data, [2, 1, 3, 4]);
utils.verbose(1, "Data shape: %dx%dx%dx%d", size(data, 1), size(data, 2), size(data, 3), size(data, 4));

% NOTE: nx and ny MUST be updated as data changes size

if isfield(p.detector, 'crop') && ~isempty(p.detector.crop)
    crop = num2cell(p.detector.crop);
    [min_x, max_x, min_y, max_y] = crop{:};
    utils.verbose(2, "Cropping data to %d:%d x %d:%d", min_x, max_x, min_y, max_y)
    data = data(:, :, min_x:max_x, min_y:max_y);

    % crop scan positions as well
    mask = false(nx, ny);
    mask(min_x:max_x, min_y:max_y) = true;
    mask = reshape(mask, [], 1);

    if ~all(mask)
        p.positions_real = p.positions_real(mask, :);
        p.positions_orig = p.positions_orig(mask, :);
        p.positions = p.positions(mask, :);

        p.numpts = size(p.positions, 1);
        p.scanidxs{1, 1} = 1:p.numpts;

        nx = size(data, 3);
        ny = size(data, 4);
        utils.verbose(1, "New size: %dx%d", nx, ny); 

        if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster"
            p.scan.nx = nx;
            p.scan.ny = ny;
        end
    end
end

if isfield(p.detector, 'step') && ~isempty(p.detector.step)
    % should broadcast
    step = zeros(1, 2);
    step(:) = p.detector.step; % x, y
    utils.verbose(1, "Using every %dx%d scan position", step(1), step(2)); 

    if any(step ~= 1)
        data = data(:, :, 1:step(1):nx, 1:step(2):ny);

        % crop scan positions as well
        mask = false(nx, ny);
        mask(1:step(1):nx, 1:step(2):ny) = true;
        mask = reshape(mask, [], 1);

        p.positions_real = p.positions_real(mask, :);
        p.positions_orig = p.positions_orig(mask, :);
        p.positions = p.positions(mask, :);

        p.numpts = size(p.positions, 1);
        p.scanidxs{1, 1} = 1:p.numpts;

        nx = size(data, 3);
        ny = size(data, 4);
        utils.verbose(1, "New size: %dx%d", nx, ny); 

        if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster"
            p.scan.nx = nx;
            p.scan.ny = ny;
            p.scan.step_size_x = p.scan.step_size_x * nx;
            p.scan.step_size_y = p.scan.step_size_y * ny;
        end
    end
end

if isfield(p.detector, 'fill_nan') && p.detector.fill_nan
    data = utils.fillnan(data, 2);
end

if isfield(p.detector, 'sim') && p.detector.sim
    if isfield(p.detector, 'beam_dose') && ~isempty(p.detector.beam_dose)
        utils.verbose(2, "Scaling data to dose: %.1f e/A^2", p.detector.beam_dose);
        if ~(isfield(p, 'src_positions') && p.src_positions == 'matlab_pos' && isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster")
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

if isfield(p.detector, 'tile') && ~isempty(p.detector.tile) && ~all(p.detector.tile == 1)
    % tile patterns (useful for simulation)
    tile = num2cell(p.detector.tile);
    [tile_x, tile_y] = tile{:};
    utils.verbose(1, "Tiling data %dx%d in realspace", tile_x, tile_y);

    if ~(isfield(p, 'src_positions') && p.src_positions == "matlab_pos" && isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster")
        error("Option 'tile' requires a raster scan");
    end

    % tile in output orientation
    data = repmat(data, 1, 1, tile_x, tile_y);
    nx = nx * tile_x;
    ny = ny * tile_y;

    % and re-generate scan
    p.scan.nx = p.scan.nx * tile_x;
    p.scan.ny = p.scan.ny * tile_y;
    utils.verbose(1, "Re-generating tiled scan...");
    p = scans.read_positions(p);
    p = core.ptycho_adjust_positions(p);

    utils.verbose(3, "Old object size: %dx%d", p.object_size(1), p.object_size(2));
    % and update object size
    p.object_size = p.asize + max(ceil(p.positions),[],1) + p.positions_pad;
    utils.verbose(2, "New object size: %dx%d", p.object_size(1), p.object_size(2));
end

if isfield(p, 'src_positions') && p.src_positions == "matlab_pos"
    % check raster scan dimensions
    if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster" && (ny ~= p.scan.ny || nx ~= p.scan.nx)
        utils.verbose(1, "Warning: Specified raster scan size %dx%d does not match actual scan size %dx%d", p.scan.nx, p.scan.ny, nx, ny);
    end
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

