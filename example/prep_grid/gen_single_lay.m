%------------------------------------------------------------------------------
%-- Example of generating the interface file 
%------------------------------------------------------------------------------

clear

%------------------------------------------------------------------------------
%-- generate x and y
%------------------------------------------------------------------------------

%-- should 6 more than FD points due to ghosts points
nghost = 3;
nx = 306;
x0 = 0.0;
dx = 100.0;

x1d = [0 : nx-1] * dx + x0 - nghost * dx;

%------------------------------------------------------------------------------
%-- generate or load topography
%------------------------------------------------------------------------------

for i=1:nx
  topo(i) = sin((i-1)*pi/nx) * 100;
end
free_topo = topo;
%------------------------------------------------------------------------------
%-- generate grid interfaces
%------------------------------------------------------------------------------

%-- note: from bottom to top

num_of_interfaces = 2; 

% 1 less than num_of_interfaces, total cell should be nghost more than FD points
num_of_cell_per_layer = [ 350 ];
dz_is_equal_of_layer  = [ 1 ]; % 1:The grid spacing is equal; 0: Is not.
avg_dz_of_layer       = [ 100 ];
smooth_length         = [ 40,1];

%-- use avg_dz and num_of_cell to esti z-axis of each interface
%-- last elem is the free surface
z_of_interfaces(num_of_interfaces) = 0;
for ilay = num_of_interfaces-1 : -1 : 1
  z_of_interfaces(ilay) = z_of_interfaces(ilay+1) ...
        - avg_dz_of_layer(ilay) * num_of_cell_per_layer(ilay);
end

%-- construct grid interfaces from free_topo and z_of_interfaces
x2d = zeros(nx,num_of_interfaces);
z2d = zeros(nx,num_of_interfaces);

%-- set x2d
for n = 1 : num_of_interfaces
  x2d(:,n) = x1d;
end

%- first same to free surface
z2d(:,num_of_interfaces) = free_topo;

for ilay = num_of_interfaces-1 : -1 : 1
  %-- smooth free_topo 
  z2d(:,ilay) = z_of_interfaces(ilay);
end

hid = figure;
set(hid,'BackingStore','on');

plot(x2d,z2d,'k-');
hold on
plot(x2d',z2d','k-');
  
xlabel(['X axis m']);
ylabel(['Z axis m']);
%==============================================================================
%-- write .gdlay file
%==============================================================================
gdlay_file = 'random_topo_single.gdlay'

gdlay_export(gdlay_file, ...
             num_of_interfaces, nx,  ...
             num_of_cell_per_layer, ...
             dz_is_equal_of_layer, ...
             x2d, z2d);

