function [ ] = save_tiff_image_stack(input, par)
% Save 3D arry as a (8-bit) tiff stack 
% Written by YJ
% ** input: 3D real array
% ** par: parameters

input = double(input);

if ~isfield(par,'axis'); par.axis = 3; end
if ~isfield(par,'range'); par.range = [0,0]; end
if ~isfield(par,'norm_single_slice'); par.norm_single_slice = false; end
if ~isfield(par,'name'); par.name = 'image_stack.tiff'; end

if par.norm_single_slice
    if any(par.range)
        range = par.range;
    else
        % range based on entire image (TODO: use quantile here?)
        range = [min(input(:)),max(input(:))];
    end
else
    % call mat2gray without limits, auto-scaling for each slice
    range = [0, 0];
end

switch par.axis
    case 3
        imwrite(mat2uint16(input(:, :, 1), range), par.name, 'tiff')
        for i=2:size(input,3)
            imwrite(mat2uint16(input(:, :, i), range), par.name, 'tiff', 'WriteMode', 'append')
        end
    case 2
        imwrite(mat2uint16(input(:, 1, :), range), par.name, 'tiff')
        for i=2:size(input,3)
            imwrite(mat2uint16(input(:, i, :), range), par.name, 'tiff', 'WriteMode', 'append')
        end
    case 1
        imwrite(mat2uint16(input(1, :, :), range), par.name, 'tiff')
        for i=2:size(input,3)
            imwrite(mat2uint16(input(i, :, :), range), par.name, 'tiff', 'WriteMode', 'append')
        end
    otherwise
        disp('Wrong axis! par.axis should be 1, 2, or 3.')
end

end

function [out] = mat2uint16(A, limits)
    if any(limits)
        scaled = mat2gray(A, limits);
    else
        scaled = mat2gray(A);
    end
    out = uint16(2^16 * scaled);
end
