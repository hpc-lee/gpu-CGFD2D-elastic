clc;
clear all;
close all;
%x direction 21 stations
nx = 41;
is_coord = 1;
is_depth = 1;
depth = 0.0;
sta_num = nx;
%origin is (0,0)
%space 1000
dh =100;
%%
%==============================================================================
%-- write .gdlay file
%==============================================================================
station_file = 'station.list';

fid=fopen(station_file,'w'); % Output file name 

%-- first line: how many stations
fprintf(fid,'%6d\n',sta_num);

%-- second line: station name and coords
% on topography surface, z set 9999 
for i = 1:nx
        sta_name = ['recv',num2str(i),'_x_z'];
        x = dh*(i-1);
        fprintf(fid,'%s %d %d  %12.2f  %12.2f\n',sta_name, is_coord, is_depth, x, depth);
    end
end
fclose(fid);
