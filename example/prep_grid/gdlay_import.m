function gd = gdlay_import(gdlay_file)

% This function import gdlay file to gd structure

fid=fopen(gdlay_file,'r'); %

%-- first line: how many interfaces
gd.num_of_interfaces = fscanf(fid,'%d', 1);

%-- second line: how many cells
gd.num_of_cell_per_layer = fscanf(fid,'%d', gd.num_of_interfaces-1);

%-- third line: is dz equal of each layer
gd.dz_is_equal_of_layer = fscanf(fid,'%d', gd.num_of_interfaces-1);

%-- 4th line: nx, ny
gd.nx = fscanf(fid,'%d', 1);

%-- others: z of interfaces
A = fscanf(fid,'%f', 2 * gd.nx * gd.num_of_interfaces);

%- reshape to 4D array
A = reshape(A,[3,gd.nx,gd.num_of_interfaces]);

%- to 
gd.x2d = squeeze(A(1,:,:));
gd.z2d = squeeze(A(3,:,:));

fclose(fid);

end % function
