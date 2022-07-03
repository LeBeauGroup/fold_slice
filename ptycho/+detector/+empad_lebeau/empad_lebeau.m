function [ det ] = empad_lebeau( p )

det.pixel_size = 1;

det.data_stored = true;

det.asize = [128 128];

det.geometry.sz = [256 256];
det.geometry.mask = [];

det.orientation = [0 0 0];

det.image_read_extraargs = {};

det.basepath_dir = '';
det.get_filename = @detector.empad_lebeau.get_filename;

det.read_path = @utils.compile_cu_dirname; 

end

