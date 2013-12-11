G=20;%true number of modes
M=5;%number of modes fit
atrue=randn(G,1).*[G:-1:1]';%3.*[G:-1:1]';atrue=atrue(1:G);
NN=0:5000;
%depth_ind=1:length(NN);
depth_ind=[1:100:3000];


for n=1:G
  Ftrue(:,n) = sin(n.*pi*NN'./NN(end));
end


F1=Ftrue(depth_ind,1:M);


eta_true=Ftrue(depth_ind,:)*atrue;


eta1=eta_true+10.*randn(size(eta_true));



ahat=F1\eta1;

normfac=atrue(1)./ahat(1);
ahat.*normfac
atrue(1:M)

figure
hold on
for n=1:M
 plot(NN(depth_ind),normfac.*ahat(n).*F1(:,n),'r')
 plot(NN,atrue(n).*Ftrue(:,n),'b')
end
title('Blue: true contribution of first M modes; Red: estimated')
xlim([depth_ind(1) depth_ind(end)])

eta_recon=(F1*ahat).*normfac;