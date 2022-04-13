

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
ny = str2double(raw_parts{end}(2:end));
nx = str2double(raw_parts{end-1}(2:end));

if ny ~= p.scan.ny || nx ~= p.scan.nx
    verbose(1, sprintf("Warning: Specified scan size %dx%d does not match filename scan size %dx%d", p.scan.nx, p.scan.ny, nx, ny));
end

% TODO better errors here
data = single(fread(fopen(file,'r'), 128*130*nx*ny,'float32'));
% resize as 4D
data = reshape(data, 128, 130, nx, ny);
% crop junk rows
data = data(1:128, 1:128, :, :);

data = permute(data, [2, 1, 4, 3]);
data = flip(data, 3);
data = flip(data, 4);
data = flip(data, 2);

utils.verbose(1, strcat('Loaded data from: ', file));

detStorage.data = data;
end

