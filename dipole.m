function V = dipole

m  = 9.42e-18;%4.38e-20;%[Am^2]
th = 0;
B0 = 10^-7*m;
d  =.00001;
A  = pi*.005^2;
f  = 1000;
mgW=.001;
N  = 6e11*mgW;



Bx=B0*((3*cos(th)^2-1)/d^3);
By=B0*(3*cos(th)*sin(th)/d^3);

B = sqrt(Bx.^2+By.^2);

flux = -B*A*cos(th);%*cos(2*pi*t*f); assume max amplitude

V = -N*flux*2*pi*f*10^6; %uV

disp(['Voltage is ' num2str(V) ' uV'])
