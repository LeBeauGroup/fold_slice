
function [p] = get_filename(p)
det = p.detectors(p.scanID).params;

read_path = p.base_path;
detStorage = p.detectors(p.scanID).detStorage;

if isfield(p, 'prepare_data_filename') && ~isempty(p.prepare_data_filename)
    filename = strcat("/", p.prepare_data_filename);
else
    filename = "/scan_x%d_y%d.raw";
end

for ii = 1:p.numscans
    detStorage.files{ii} = strcat(p.raw_data_path, sprintf(filename, p.scan.nx, p.scan.ny));
end
end