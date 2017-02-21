function S = ScaleCorr(f,Hr,Hb,NumMonoSeg,flag)
% function S = ScaleCorr(f,Hr,Hb,NumMonoSeg,flag)
% f - frequencies at which functions are given
% Hr- first function
% Hb- seconf function
% numMonoSeg - number of monotonic segments
% flag - "C" (default) aligns common points - "P" uses polynomial fits to align

if nargin<4; NumMonoSeg=1; end; OutFlag=0;
if nargin<5; flag='C'; end; if length(flag)~=1; flag='C'; end; if flag~='P'; flag='C'; end;

if size(f ,1)<size(f ,2); f = f'; end;  % make them column vectors
if size(Hr,1)<size(Hr,2); Hr=Hr'; end;
if size(Hb,1)<size(Hb,2); Hb=Hb'; end;

if NumMonoSeg==1 & flag=='C';   if OutFlag==1; disp('One monotonic segment'); end;
  % Monotonic harmonics
  in=find(Hb<=max(Hr) & Hb>=min(Hr));
  if ~isempty(in);  % do the fit on the basis of the overlapping points if there are any
    fb=interp1(Hr,f,Hb,'spline');
    S=(f(in)'*fb(in))/(f(in)'*f(in));
  else;
    B=[f,ones(size(f))]; LC=lscov(B,Hr); fb=(Hb-LC(2))./LC(1); S=mean(fb./f);
    disp('exterior points only'); if OutFlag==1; plot(f,Hr,':',f.*S,Hb,'d'); end;
  end;
end;

if NumMonoSeg==2;   if OutFlag==1; disp('Two monotonic segments'); end;
  % Find the signs of the derivative at each point
  dfp=[diff(Hr);0]; dfn=[0;-flip(diff(flip(Hr)))];
  df=zeros(size(Hr));
  inp=find(dfp>=0 & dfn>=0); if ~isempty(inp); df(inp)=+1; end;
  inn=find(dfp<=0 & dfn<=0); if ~isempty(inn); df(inn)=-1; end;
   if ~isempty(inp); bf=Consecutive(df',+1)'; inHrp=find(bf==1 | (bf==0&df==0)); else; inHrp=[]; end;
   if ~isempty(inn); bf=Consecutive(df',-1)'; inHrn=find(bf==1 | (bf==0&df==0)); else; inHrn=[]; end;
  
  dfp=[diff(Hb);0]; dfn=[0;-flip(diff(flip(Hb)))];
  dfb=zeros(size(Hb));
  inp=find(dfp>=0 & dfn>=0); if ~isempty(inp); dfb(inp)=+1; end;
  inn=find(dfp<=0 & dfn<=0); if ~isempty(inn); dfb(inn)=-1; end;
   if ~isempty(inp); bf=Consecutive(dfb',+1)'; inHbp=find(bf==1 | (bf==0&df==0)); else; inHbp=[]; end;
   if ~isempty(inn); bf=Consecutive(dfb',-1)'; inHbn=find(bf==1 | (bf==0&df==0)); else; inHbn=[]; end;

  % Interpolate the frequencies for each monotonic segment
  Fb=[]; Fbint=[];
  if ~isempty(inHrp) & ~isempty(inHbp);
    Hrt=Hr(inHrp); fr=f(inHrp); Hbt=Hb(inHbp);
    in=find(Hbt<=max(Hrt) & Hbt>=min(Hrt));
    Fbint=interp1(Hrt,fr,Hbt(in)); Fb=f(inHbp(in));
  end;

  if ~isempty(inHrn) & ~isempty(inHbn);
    Hrt=Hr(inHrn); fr=f(inHrn); Hbt=Hb(inHbn); fb=f(inHbn);
    in=find(Hbt<=max(Hrt) & Hbt>=min(Hrt));
    bf=interp1(Hrt,fr,Hbt(in)); Fbint=[Fbint;bf]; Fb=[Fb;f(inHbn(in))];
  end;

  % Find the least squares scale
  S=(Fbint'*Fb)/(Fb'*Fb);

end;

if NumMonoSeg==1 & flag=='P';  if OutFlag==1; disp('Polynomial fit, one monotonic segment'); end;
  B=[f.^2,f,ones(size(f))]; H=[Hr;Hb];
  SL=0; SU=10; Tol=1e-4;
  N=0; NL=0; NU=0;
  while NL*NU==0;
  while SU-SL>Tol; if OutFlag==1; disp([SL,SU,N,NL,NU]); end;
   Sl=SL+0.382*(SU-SL); Sr=SL+0.618*(SU-SL);
   S=[[Sl^2,0,0];[0,Sl,0];[0,0,1]]; FM=[B;B*S]; Coef=lscov(FM,H); El=normf(H-FM*Coef);
   S=[[Sr^2,0,0];[0,Sr,0];[0,0,1]]; FM=[B;B*S]; Coef=lscov(FM,H); Er=normf(H-FM*Coef);
   if OutFlag==1; disp([El,Er]); end;
   if El<Er; SU=Sr; Sr=Sl; Sl=SL+0.382*(SU-SL); NU=NU+1;
       else; SL=Sl; Sl=Sr; Sr=SL+0.618*(SU-SL); NL=NL+1; end; N=N+1;
  end;   if OutFlag==1; disp('end first loop'); end;
  if OutFlag==1; disp([NL,NU]); end;
  if NL==0; SL=SL/2; end;
  if NU==0; SU=SU*2; end;
  if NL==0 | NU>1000; disp('Optimization Not Successful'); break; end;
  end;    % disp('end second loop');
  S=(SU+SL)/2;
end;


return;

C=lscov([f.^2',f',ones(size(f))'],mr(:,1));
plot(f,mr(:,1),'o',f,[f.^2',f',ones(size(f))']*C,'-');
normf(mr(:,1)-[f.^2',f',ones(size(f))']*C)
normf(mr(:,1)-[f.^2',f',ones(size(f))']*C)/normf(mr(:,1))




FM=[B;B*1]; Coef=lscov(FM,H); normf(H-FM*Coef)



if NumMonoSeg==1;
%  Monotonic harmonics
in=find(Hb<=max(Hr) & Hb>=min(Hr));

if ~isempty(in);  % do the fit on the basis of the overlapping points if there are any
    fb=interp1(Hr,f,Hb);
    S=mean(fb(in)./f(in)); disp('interior points');  disp(S); plot(f(in),Hr(in),':o',fb(in),Hb(in),'--d');
    Lo=min(fb(in)./f(in))*.5; Uo=max(fb(in)./f(in))*2; Tol=.001;
    L=Lo; U=Uo;
    N=0; NL=0; NU=0;
    while NL*NU==0;
    while U-L>Tol;    if OutFlag==1; disp(round([L,U,N,NL,NU])); end;
     l=L+0.382*(U-L); u=L+0.618*(U-L);
     El=normf(fb(in)-l*f(in));  % El=MError;
     Eu=normf(fb(in)-u*f(in));  % Er=MError;
     if OutFlag==1; disp('Error Bounds'); disp(round([El,Eu])); end;
     if El<Eu; U=u; u=l; l=L+0.382*(U-L); NU=NU+1;
         else; L=l; l=u; u=L+0.618*(U-L); NL=NL+1; end; N=N+1;
    end;    disp('end first loop');
     if OutFlag==1; disp([NL,NU]); end;
     if NL==0; L=L/2; end;
     if NU==0; U=U*2; end;
     if NL==0 | NU>100000; break; end;
    end;  if OutFlag==1; disp('end second loop'); end;
    S=(U+L)/2;
    if OutFlag==1; disp('Data is fit'); plot(f(in),Hr(in),':',f(in)*S,Hb(in),'o'); end;
else;
    B=[f,ones(size(f))]; LC=lscov(B,Hr); fb=(Hb-LC(2))./LC(1); S=mean(fb./f);
    disp('exterior points only'); plot(f,Hr,':',f,LC(1)*f+LC(2),'*',f,Hb,'--*',f.*S,Hb,'d');
end;
end;



return




if NumMonoSeg==2;
  [Cr,Hfr]=QuadFit(Hr,f);
  in=find(Hb<=max(Hr) & Hb>=min(Hr));
  [Cb,Hfb]=QuadFit(Hb(in),f(in));
  fc=f(in); Nf=length(fc);
  B=[fc.^2,fc,ones(Nf,1)]*Cr; A=[fc.^2*Cb(1),fc*Cb(2),ones(Nf,1)*Cb(3)];
    S=zeros(Nf,2);
  for k=1:Nf;
    a=A(k,1); b=A(k,2); c=A(k,1)-B(k);
    S(k,1)=(-b+sqrt(b^2-4*a*c))/(2*a);
    S(k,2)=(-b-sqrt(b^2-4*a*c))/(2*a);
  end;
  S=reshape(S,2*Nf,1); in=find(imag(S)~=0); S(in)=[];
  % minEr=normf(A*[S(1)^2;S(1);1]-B); inM=1;
  % for k=2:length(S);
  %   tEr=normf(A*[S(k)^2;S(k);1]-B);
  %   if tEr<minEr; minEr=tEr; inM=k; end;
  % end;
    Lo=min(S)*.5; Uo=max(S)*2; Tol=.001; S=zeros(1);
    L=Lo; U=Uo; keyboard;
    N=0; NL=0; NU=0;
    while NL*NU==0;
    while U-L>Tol;    if OutFlag==1; disp(round([L,U,N,NL,NU])); end;
     l=L+0.382*(U-L); u=L+0.618*(U-L);
     El=normf(A*[l^2;l;1]-B);  % El=MError;
     Eu=normf(A*[u^2;u;1]-B);  % Er=MError;
     if OutFlag==1; disp('Error Bounds'); disp(round([El,Eu])); pause; end;
     if El<Eu; U=u; u=l; l=L+0.382*(U-L); NU=NU+1;
         else; L=l; l=u; u=L+0.618*(U-L); NL=NL+1; end; N=N+1;
    end;    disp('end first loop');
     if OutFlag==1; disp([NL,NU]); end;
     if NL==0; L=L/2; end;
     if NU==0; U=U*2; end;
     if NL==0 | NU>100000; break; end;
    end;  if OutFlag==1; disp('end second loop'); end;
    S=(U+L)/2;
    if OutFlag==1; disp('Data is fit'); plot(f(in),Hr(in),':o',fb(in)*S,Hb(in),'--d'); pause; end;
end;

return


% Outer loop makes sure that both limits change; if they don't bounds were wrong
% Inner loop optimizes objective function within the bounds, Lo and Uo
L=Lo; U=Uo;
N=0; NL=0; NU=0;
while NL*NU==0;
while U-L>Tol;    if OutFlag==1; disp(round([L,U,N,NL,NU])); end;
 l=L+0.382*(U-L); u=L+0.618*(U-L);
 El=MError;
 Eu=MError;   if OutFlag==1; disp(round([El,Eu])); end;
 if El<Eu; U=u; u=l; l=L+0.382*(U-L); NU=NU+1;
     else; L=l; l=u; u=L+0.618*(U-L); NL=NL+1; end; N=N+1;
end;    disp('end first loop');
if OutFlag==1; disp([NL,NU]); pause; end;
if NL==0; L=L/2; end;
if NU==0; U=U*2; end;
if NL==0 | NU>100000; break; end;
end;  if OutFlag==1; disp('end second loop'); end;





  %S=[1:200]/100; for k=1:32; oEr(k)=normf(B-A*[S(k)^2;S(k);1]); end;
  %S=1; oEr=normf(B-A*[S^2;S;1]); dS=.01; Stop=0;
  %while Stop==0; S=S+dS; nEr=normf(B-A*[S^2;S;1]); if nEr>oEr & Stop==0; ; end; end;
  

% Monotonic harmonics
% in=find(Hb<=max(Hr) & Hb>=min(Hr)); Win=zeros(size(Hb)); Win(in)=1; [a,b,c]=LocalError(Hr,Hb,Win);
% s=find(c==max(c)); bf=shift1d(Hb,-s-1);
% delx=corr1d(Hb(in),Hr,nummax); bf=shift1d(,delx);

% if ~isempty(in);  % do the fit on the basis of the overlapping points if there are any
%     fb=interp1(Hr,f,Hb);
%     a=mean(fb(in)./f(in)); disp('interior points');  plot(f(in),Hr(in),':o',fb(in),Hb(in),'--d');
% else;
%     B=[f,ones(size(f))]; LC=lscov(B,Hr); fb=(Hb-LC(2))./LC(1); a=mean(fb./f);
%     disp('exterior points only'); plot(f,Hr,':o',f,LC(1)*f+LC(2),'-.*',f,Hb,'--*',f.*a,Hb,'d');
% end;
% end;



