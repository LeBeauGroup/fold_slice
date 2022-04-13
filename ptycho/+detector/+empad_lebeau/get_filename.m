
function [p] = get_filename(p)
det = p.detectors(p.scanID).params;
detStorage = p.detectors(p.scanID).detStorage;

if isfield(p, 'raw_data_path') && ~isempty(p.raw_data_path)
    read_path = p.raw_data_path;
else
    read_path = p.base_path;
end

if isfield(p, 'raw_data_filename') && ~isempty(p.raw_data_filename)
    filename = strcat("/", p.raw_data_filename);
else
    filename = "/scan_x%d_y%d.raw";
end

for ii = 1:p.numscans
    detStorage.files{ii} = strcat(read_path, sprintf(filename, p.scan.nx, p.scan.ny));
end
end