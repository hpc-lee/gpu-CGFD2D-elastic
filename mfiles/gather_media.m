function [v] = gather_media(parfnm,output_dir,varnm,subs,subc,subt)

% load
fnm_media=[output_dir,'/','media.nc'];
    
if ~ exist(fnm_media,'file')
   error([mfilename ': file ' fnm_media 'does not exist']);
end

xzs = nc_attget(fnm_media,nc_global,'local_index_of_first_physical_points');
xzs = double(xzs);
xs = subs(1) -1 + xzs(1); 
zs = subs(2) -1 + xzs(2); 

xzc = nc_attget(fnm_media,nc_global,'count_of_physical_points');
xzc = double(xzc);

if(subc(1) == -1)
  xc = floor(xzc(1)/subt(1))-subs(1)+1;
else
  xc = subc(1);
end
if(subc(2) == -1)
  zc = floor(xzc(2)/subt(2))-subs(2)+1;
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

v(k1:k2,i1:i2)=nc_varget(fnm_media,varnm,[zs,xs],[zc,xc],[zt,xt]);

v=v';

end
