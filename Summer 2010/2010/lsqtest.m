function [] = lsqtest(n_modes,n_depths,drange)

	for i=1:n_modes
        F2(:,i) = sin(i*pi*[0:n_depths]'/(n_depths+1));
        FR(:,i) = sin((i)*pi*[0:n_depths]'/(n_depths+1));    
    end
    
	F3 = F2(2:n_depths+1,:);
	eta = zeros(n_depths,1);
	
	At = randn(1000,n_modes);
	
	eta = At*FR(2:n_depths+1,:)';
    
	At1 = eta(:,drange)*pinv(F3(drange,:)');
	
	error = (At1-At)./At*100
    plot(At1);  
    
    coher(At1(:,1),At1(:,2),10/24/60,4,'1','2',1);