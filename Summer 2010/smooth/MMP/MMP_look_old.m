%Look at eq'l Atlantic MMP data (from pilot exp't)
%
%









%load('D:\eql_atlantic_MMP\equatorial_atlantic\grd0149.mat')
%load('D:\eql_atlantic_MMP\equatorial_atlantic\grd0150.mat')








if 1==2%read profiles

pmin=1;
pmax=5500;
dec=1;



pathstr='D:\eql_atlantic_MMP\equatorial_atlantic';
files=dir(pathstr);
files(1:2)=[];
%filestr=files(:).name;

for jj=1:length(files)
  filename=files(jj).name;
  filestr(jj,:)=([pathstr '\' filename]);
end


[time,ss,tt,pp,pt,sg,u,v,w,dpdt,tm] = get_mmp(filestr,dec,pmin,pmax);
% OUTPUTS:
%	  time: the start time of each cast
%         ss:  rectangular array of salinity values (dimensions pp x time) 
%         tt:  rectangular array of temperature values  
%         pp:  vector of pressure level
%         pt:  rectangular array of potential temperature values
%	  sg:  rectangular array of sigma-theta values
%	   u:  rectangular array of velocity values  (dimensions acmpp x time)
%	   v:  rectangular array of velocity values 
%	   w:  rectangular array of vertical velocity values
%         tm:  rectangular array of time data for each data bin



figure
%plot(time(2:end),diff(time))
plot(diff(time))
ylabel('Time between samples (days)')

figure
%plot(time(2:end),diff(time))
plot(time-time(1),'.')
ylabel('Sample time (days)')



figure
pcolor(time,pp,u)
shading flat
set(gca,'ydir','reverse')
title('U')


figure
pcolor(time,pp,v)
shading flat
set(gca,'ydir','reverse')
title('V')

comment_mmp=['MMP profile data merged in MMP_look.m on ' date]
save eql_Atl_MMP_raw time ss tt pp pt sg u v w dpdt tm comment_mmp

else
load eql_Atl_MMP_raw
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid data in time:

pgrid=pp(1,1):4:pp(end,1);
tgrid=time(1):2:time(end);
%[hi,W] = loess2d_semi_regular2(x,t,h,Sx,St,xi,ti,xdim)

Sx=6;%half-power points are ~0.6*S
St=max(diff(time)).*2./0.6;

[ugrid] = loess2d_semi_regular2(pp,tm,u,Sx,St,pgrid,tgrid,1);
[vgrid] = loess2d_semi_regular2(pp,tm,v,Sx,St,pgrid,tgrid,1);
[Sgrid] = loess2d_semi_regular2(pp,tm,ss,Sx,St,pgrid,tgrid,1);
[Tgrid] = loess2d_semi_regular2(pp,tm,tt,Sx,St,pgrid,tgrid,1);


figure
pcolor(tgrid,pgrid,ugrid)
shading flat
set(gca,'ydir','reverse')
title('U')
caxis([-30 30])

figure
pcolor(tgrid,pgrid,vgrid)
shading flat
set(gca,'ydir','reverse')
title('V')
caxis([-30 30])

figure
pcolor(tgrid,pgrid,Sgrid)
shading flat
set(gca,'ydir','reverse')
title('S')
caxis([34.6 35.1])

figure
pcolor(tgrid,pgrid,Tgrid)
%contour(tgrid,pgrid,Tgrid,[4.5 4.55],'k')
shading flat
set(gca,'ydir','reverse')
title('T')
caxis([2 5.5])


figure
pcolor(tgrid,pgrid,ugrid)
shading flat
set(gca,'ydir','reverse')
title('U')
caxis([-30 30])
hold on
contour(tgrid,pgrid,Tgrid,[4.5 4.58],'k')
contour(tgrid,pgrid,Sgrid,[1 1].*34.975,'m')



save eql_Atl_MMP_gridded tgrid pgrid Sgrid Tgrid ugrid vgrid Sx St



load eql_Atl_MMP_gridded





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid and incorporate ADCP data
% 
load('D:\eql_atlantic_MMP\kiel_adcp_ctd\data_KPO1001_all\u_lr_23w_2006_2008.mat')

yday_adcp=datenum(gregorian(jd));
zadcp=ones(length(yday_adcp),1)*z;
tadcp=yday_adcp*ones(1,length(z));

figure
imagesc(jd,z,u')


[ugrid2] = loess2d_semi_regular2(zadcp,tadcp,u,Sx,St,z,tgrid,1);
[vgrid2] = loess2d_semi_regular2(zadcp,tadcp,v,Sx,St,z,tgrid,1);



figure
%pcolor(tgrid,z,ugrid2.*100)
pcolor(tgrid,z,ugrid_adcp)
hold on
pcolor(tgrid,pgrid,ugrid)
shading flat
set(gca,'ydir','reverse')
title('U')
caxis([-1 1].*80)
axis tight


figure
%pcolor(tgrid,z,vgrid2.*100)
pcolor(tgrid,z,vgrid_adcp)
hold on
pcolor(tgrid,pgrid,vgrid)
shading flat
set(gca,'ydir','reverse')
title('V')
caxis([-1 1].*80)
axis tight


ugrid_adcp=100.*ugrid2;
vgrid_adcp=100.*vgrid2;
comment=['MMP/ADCP data smoothed in MMP_look.m on ' date]

save eql_Atl_MMP_ADCP_gridded comment tgrid pgrid Sgrid Tgrid ugrid vgrid ugrid_adcp vgrid_adcp z Sx St





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Load current meter data and smooth in same way as MMP data

load('D:\eql_atlantic_MMP\kiel_adcp_ctd\data_KPO1001_all\vel_687m_842m_998m_mid_depth.mat')









