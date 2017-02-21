function A = initrand(N)

%make matrix of normalized random vectors
 A = zeros(N,3); 
 
 %%
 for ii=1:N
 
     r = 2*rand(1,3)-1;
     %r = rand(1,3);
 a=r/norm(r);
 
 A(ii,:) = a;
 end
