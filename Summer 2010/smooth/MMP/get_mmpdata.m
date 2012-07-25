function [stas,time,ss,tt,pp,pt,sg,u,v,w,dpdt,tm] = get_mmpdata;
%function [stas,time,ss,tt,pp,pt,sg,u,v,w,dpdt,tm] = get_mmpdata;
%
% driver for get_MMP
%  routine to load multiple profiles of gridded MMP data from one directory

arch=computer;
if(strcmp(arch,'PCWIN'));
   slsh='\';
 else
   slsh='/';
end

   

stas=input('Enter the range of profiles that you want to load ');
directory=input(' Enter the name of the data directory ','s');
pmin=input('Enter the min pressure to access ');
pmax=input('Now the max pressure to load ');
decimation=input('Lastly, enter the decimation rate in dbars ');

for i=1:length(stas)
flnames(i,:)=[directory,slsh,'grd',sprintf('%4.4i',stas(i)),'.mat'];
end

[time,ss,tt,pp,pt,sg,u,v,w,dpdt,tm] = get_MMP(flnames,decimation,pmin,pmax);
