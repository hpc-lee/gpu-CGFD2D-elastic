function [v,t]=gather_snap(parfnm,output_dir,nlayer,varnm,subs,subc,subt)

% load
fnm_snap=[output_dir,'/','volume_vel.nc'];

% check dir exists
if ~ exist(fnm_snap,'file')
    error([mfilename ': file ' fnm_snap ' does not exist']);
end

tdim=nc_getdiminfo(fnm_snap,'time');
if tdim.Length==0 | (nlayer-1)-1>=tdim.Length
   error([num2str(nlayer) 'th layer is beyond current time dim (' ...
        num2str(tdim.Length) ') in ' fnm_snap]);
end

% read parameters file
par=loadjson(parfnm);
snap_subc=par.snapshot{1}.grid_index_count;
snap_subc = double(snap_subc);

xs = subs(1) -1; 
zs = subs(2) -1; 
if(subc(1) == -1)
  xc = ceil((snap_subc(1)-subs(1)+1)/subt(1));
else
  xc = subc(1);
end
if(subc(2) == -1)
  zc = ceil((snap_subc(2)-subs(2)+1)/subt(2));
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
% get data
v(k1:k2,i1:i2)=nc_varget(fnm_snap,varnm, ...
      [nlayer-1,zs,xs],[1,zc,xc],[1,zt,xt]);
t=nc_varget(fnm_snap,'time',[nlayer-1],[1]);

%v=v';

end
