clear all;
close all;
clc;
addmypath;
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project1/test.json';
output_dir='../project1/output';

% which variable to plot
varnm='Vx';
% which station to plot (start from index '1')
startid=1;
endid = 7;

% figure control parameters
flag_print=0;

% ---------------------------------------------------------------------- %

% read parameter file
par=loadjson(parfnm);

fileID = fopen(par.in_source_file);
recvprefix = fgetl(fileID);
while(recvprefix(1) == "#")
    recvprefix = fgetl(fileID);
    if(recvprefix(1) ~= "#")
        break;
    end
end
fclose(fileID);



fileID = fopen(par.in_station_file);
%first line is number recv or station
%must read to skip
for i=1:startid
    recvnum = fgetl(fileID);
    while(recvnum(1) == "#")
    recvnum = fgetl(fileID);
    end
end
% load data
for irec=startid:1:endid
    recvinfo = fgetl(fileID);
    while(recvinfo(1) == "#")
        recvinfo = fgetl(fileID);
    end
    recvinfo = strsplit(recvinfo);
    recvnm = char(recvinfo(1));
    sacnm=[output_dir,'/',recvprefix,'.',recvnm,'.',varnm,'.sac'];
    sacdata=rsac(sacnm);
    seismodata(:,irec-startid+1)=sacdata(:,2)/max(abs(sacdata(:,2)));
%     seismodata(:,irec-startid+1)=sacdata(:,2);
    seismot(:,irec-startid+1)=sacdata(:,1);
end
% plot receiver
for irec=startid:1:endid
    figure(irec-startid+1)
    plot(seismot(:,irec-startid+1),seismodata(:,irec-startid+1),'b','linewidth',1.0);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title([varnm, ' recv No.',num2str(irec),' interpreter ','yes']);
    set(gcf,'color','white','renderer','painters');
    % save and print figure
    if flag_print
        width= 800;
        height=400;
        set(gcf,'paperpositionmode','manual');
        set(gcf,'paperunits','points');
        set(gcf,'papersize',[width,height]);
        set(gcf,'paperposition',[0,0,width,height]);
        print(gcf,[varnm,'_rec_no',num2str(irec),'.png'],'-dpng');
    end
end

seismodata1 = seismodata;
seismot1 = seismot;
save('seismodata1.mat','seismodata1');
save('seismot1.mat','seismot1');

% seismodata2 = seismodata;
% seismot2 = seismot;
% save seismodata2.mat;
% save seismot2.mat;
