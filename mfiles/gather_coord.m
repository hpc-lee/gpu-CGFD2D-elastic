function [x,z] = gather_coord(parfnm,output_dir,subs,subc,subt)

% load
fnm_coord=[output_dir,'/','coord.nc'];
    
if ~ exist(fnm_coord,'file')
   error([mfilename ': file ' fnm_coord 'does not exist']);
end

xzs = nc_attget(fnm_coord,nc_global,'local_index_of_first_physical_points');
xzs = double(xzs);
xs = subs(1) -1 + xzs(1); 
zs = subs(2) -1 + xzs(2); 

xzc = nc_attget(fnm_coord,nc_global,'count_of_physical_points');
xzc = double(xzc);

if(subc(1) == -1)
  xc = ceil((xzc(1)-subs(1)+1)/subt(1));
else
  xc = subc(1);
end
if(subc(2) == -1)
  zc = ceil((xzc(2)-subs(2)+1)/subt(2));
else
  zc = subc(2);
end
%stride
xt = subt(1);
zt = subt(2);

i1 = 1;
i2 = i1 + xc - 1;
k1 = 1;
k2 = k1 + zc - 1;

% check dimension size of x,y,z to detmine cart or curv
xvar_info = ncinfo(fnm_coord,'x');
num_dim_x = length(xvar_info.Dimensions);
if num_dim_x == 1 % cart grid
  x1d = nc_varget(fnm_coord,'x',xs-1,xc,xt);
  x3d = repmat(x1d,1,(k2-k1+1));
  x(k1:k2,i1:i2) = x3d';
else % curv grid
  x(k1:k2,i1:i2)=nc_varget(fnm_coord,'x',[zs,xs],[zc,xc],[zt,xt]);
end

%- z coord
zvar_info = ncinfo(fnm_coord,'z');
num_dim_z = length(zvar_info.Dimensions);
if num_dim_z == 1 % cart grid
  z1d = nc_varget(fnm_coord,'z',zs-1,zc,zt);
  z3d = repmat(z1d,1,(i2-i1+1));
  z(k1:k2,i1:i2) = z3d;
else % curv grid
  z(k1:k2,i1:i2)=nc_varget(fnm_coord,'z',[zs,xs],[zc,xc],[zt,xt]);
end
    

x=x';
z=z';

end
