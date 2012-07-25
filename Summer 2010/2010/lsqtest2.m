G=20;%true number of modes
M=5;%number of modes fit


for ii=1:1
atrue=randn(G,1).*[G:-1:1]';%3.*[G:-1:1]';atrue=atrue(1:G);
NN=0:5000;
%depth_ind=1:length(NN);
depth_ind=[1:50:500];


for n=1:G
  Ftrue(:,n) = cos(n.*pi*NN'./NN(end));
end


F1=Ftrue(depth_ind,1:M);


eta_true=Ftrue(depth_ind,:)*atrue;


eta1=eta_true+10.*randn(size(eta_true));



ahat(:,ii)=F1\eta1;
%ahat1(:,ii) = F1\eta1;

normfac(ii)=atrue(1)./ahat(1,ii);
%normfac1(ii)=atrue(1)./ahat1(1,ii);

%atrue(1:M)

figure
hold on
for n=1:M
    %plot(NN(depth_ind),normfac1.*ahat1(n).*F1(:,n),'r')
    plot(NN(depth_ind),normfac.*ahat(n).*F1(:,n),'g')
    plot(NN,atrue(n).*Ftrue(:,n),'b')
end
title('Blue: true contribution of first M modes; Red: estimated')

end
%plot(ahat.*repmat(normfac,M,1) - ahat1.*repmat(normfac1,M,1));
coher(ahat(1,:),ahat(2,:),1,5,'1','2');
%coher(ahat1(1,:),ahat1(2,:),1,5,'1','2');