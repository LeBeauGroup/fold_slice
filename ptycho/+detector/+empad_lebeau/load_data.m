

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
ny = str2double(x_part{1}(2:end));
nx = str2double(y_part{1}(2:end));

% TODO better errors here
data = single(fread(fopen(file,'r'), 128*130*nx*ny,'float32'));
% resize as 4D
data = reshape(data, 128, 130, nx, ny);
% crop junk rows
data = data(1:128, 1:128, :, :);

data = permute(data, [1, 2, 4, 3]); % in ptychoshelves y is fast scan dir
data = flip(data, 2); % reciprocal-space y dimension is flipped

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
    data = data ./ 580;
end

if isfield(p.detector, 'crop') && ~isempty(p.detector.crop)
    crop = num2cell(p.detector.crop);
    [min_x, max_x, min_y, max_y] = crop{:};
    utils.verbose(2, "Cropping data to %d:%d x %d:%d", min_x, max_x, min_y, max_y)
    data = data(:, :, colon(min_y, max_y), colon(min_x, max_x));
    ny = size(data, 3);
    nx = size(data, 4);
end

if isfield(p, 'src_positions') && p.src_positions == "matlab_pos"
    % check raster scan dimensions
    if isfield(p, 'scan') && isfield(p.scan, 'type') && p.scan.type == "raster" && (ny ~= p.scan.ny || nx ~= p.scan.nx)
        utils.verbose(1, "Warning: Specified raster scan size %dx%d does not match actual scan size %dx%d", p.scan.nx, p.scan.ny, nx, ny);
    end
end

utils.verbose(1, strcat('Loaded data from: ', file));

detStorage.data = data;
end

