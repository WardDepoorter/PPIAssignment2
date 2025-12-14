function [TopoMx] = img2mx(filename,px)
fid = fopen(filename, 'r', 'ieee-le');
img = fread(fid, [360*px, 180*px+1], 'int16'); 
fclose(fid);

img = img';
%height = double(img) * 0.5;   % scaling factor(see https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/DATA/LOLA_GDR/CYLINDRICAL/IMG/LDEM_4.LBL)
TopoMx = img;

end

