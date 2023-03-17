
function [p] = get_filename(p)
det = p.detectors(p.scanID).params;
detStorage = p.detectors(p.scanID).detStorage;

if isfield(p, 'raw_data_path') && ~isempty(p.raw_data_path)
    read_path = p.raw_data_path;
else
    read_path = p.base_path;
end

if isfield(p, 'raw_data_filename') && ~isempty(p.raw_data_filename)
    filename = p.raw_data_filename;
else
    filename = "scan_x%d_y%d.raw";
end

if isfield(p, 'scan') && isfield(p.scan, 'nx') && isfield(p.scan, 'ny')
    filename = sprintf(filename, p.scan.nx, p.scan.ny);
elseif contains(filename, '%d')
    error("scan.nx and scan.ny required to format empad filename '%s'", filename);
end

for ii = 1:p.numscans
    detStorage.files{ii} = strcat(read_path, "/", filename);
end
end