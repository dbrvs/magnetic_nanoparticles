 function tau = getTAU(t,M1)
t=reshape(t,length(t),1);
[a] = polyfit(t,log(M1(:,1)),1);
tau = abs(1/a(1)); %us