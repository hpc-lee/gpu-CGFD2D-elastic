clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project1/test.json';
output_dir='../project1/output';

% which grid profile to plot
subs=[5,5];     % start from index '1'
subc=[-1,-1];   % '-1' to plot all points in this dimension
subt=[2,2];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_title  = 1;
scl_daspect = [1 1 1];
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

[x,z]=gather_coord(parfnm,output_dir,subs,subc,subt);

%- set coord unit
flag_km     = 1;
if flag_km
   x=x/1e3;
   z=z/1e3;
   str_unit='km';
else
   str_unit='m';
end

%-----------------------------------------------------------
%-- set figure
%-----------------------------------------------------------
hid = figure;
set(hid,'BackingStore','on');

plot(x,z,'k-');
hold on
plot(x',z','k-');
  
xlabel(['X axis (' str_unit ')']);
ylabel(['Z axis (' str_unit ')']);
  
set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis tight;

% title
if flag_title
        gridtitle='XOZ-Grid';
    title(gridtitle);
end

% save and print figure
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[gridtitle '.png'],'-dpng');
end
