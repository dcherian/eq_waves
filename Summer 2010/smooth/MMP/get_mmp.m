function [time,ss,tt,pp,pt,sg,u,v,w,dpdt,tm] = get_mmp(flnames,decimation,pmin,pmax);
%function [time,ss,tt,pp,pt,sg,u,v,w,dpdt,tm] = get_mmp(flnames,decimation,pmin,pmax);
%
% loads in selected  stations of gridded MMP data from user specified 
%        list of files
%
% INPUTS:
%   flnames: list of file names to be loaded, specify full path
%   decimatation: (opt) decimate the data to this pressure interval
%                    default is 20
%   pmin and pmax are the min and max pressures of the output data 
%
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
%
% 4/23/97 Gwyneth Hufford
%3/1/01 JMT

ss = []; tt = [] ; pp = []; sg = []; pt = []; 
u = []; v = []; w = []; dpdt=[]; time = []; tm = [];

deci=decimation;

for st =1:size(flnames,1)
  disp(['Loading ',flnames(st,:)]);

if(exist(flnames(st,:))==2);

  eval(['load ' flnames(st,:)])

%first, extract the data between pmin and pmax

j=find(pgrid>=pmin & pgrid<=pmax);

pgrid=pgrid(j);
tave=tave(j);
thetave=thetave(j);
s_ave=s_ave(j);
sigthave=sigthave(j);
uave=uave(j);
vave=vave(j);
wave=wave(j);
dpdtave=dpdtave(j);
ctimave=ctimave(j);


%perform decimation
  l=length(pgrid);
  dp=abs(pgrid(1)-pgrid(2));	
  step=round(deci/dp);
     if(step<=1)
	step =1;
     end
  PRES=pgrid(1:step:l); PRES=PRES(:);
  TEMP=tave(1:step:l); TEMP=TEMP(:);
  PTMP=thetave(1:step:l); PTMP=PTMP(:);
  SALT=s_ave(1:step:l); SALT=SALT(:);
  SIGT=sigthave(1:step:l); SIGT=SIGT(:);
  U=uave(1:step:l);  U=U(:); 
  V=vave(1:step:l);  V=V(:);
  W=wave(1:step:l);  W=W(:);
  DPDT=dpdtave(1:step:l); DPDT=DPDT(:);
  TM=ctimave(1:step:l); TM=TM(:);

	else;  %if file doesn't exist
  PRES=nan;
  TEMP=nan;
  PTMP=nan;
  SALT=nan;
  SIGT=nan;
  U=nan; 
  V=nan;
  W=nan;
  DPDT=nan;
  TM=nan;
  startdaytime=nan;

	end;  %if exist loop


%attach these values and vectors to a matrix for the set
    pp = merge(pp,PRES);   
    tt = merge(tt,TEMP); 
    ss = merge(ss,SALT);
    pt = merge(pt,PTMP);
    sg = merge(sg,SIGT); 
     u = merge(u,U);
     v = merge(v,V);
     w = merge(w,W);
    dpdt=merge(dpdt,DPDT);
     tm = merge(tm,TM);

    time=merge(time,startdaytime);

end	%end looping on files

%[nz,np]=size(pp);
%longest=find(~isnan(pp(nz,:)));
%pp=pp(:,longest(1));
