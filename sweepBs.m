function sweepBs

T2=500; T1=300;

for ii=1:6; 
    Bs(ii)=ii/3;
    
    [hx1(ii,:),rc(ii)]=Brown(Bs(ii),T1); 
    [hx2(ii,:),rh(ii)]=Brown(Bs(ii),T2); 
end

for j=1:ii; 
    hh1(j)=norm(hx1(j,:));
    hh2(j)=norm(hx2(j,:)); 
end

subplot(1,2,1); plot(Bs,hh1,Bs,hh2,Bs*T1/T2,hh2,'o')
subplot(1,2,2); plot(Bs,rc,Bs,rh,Bs*T1/T2,rh,'o')
