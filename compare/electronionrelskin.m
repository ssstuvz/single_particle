% electron movement
close all
clear all
dt=0.005;
dtdiv2=dt/2;

a1=10;

tau=4*2.67;
N=30000;

q=-1; % charge
Z=0;

x(1)=0.0;
y(1)=0;
z(1)=0;
px(1)=0;
py(1)=0;
pz(1)=0;
recgamma(1)=1;
tstart1=60;
dens=400;
ls=1/sqrt(dens);
kappa=sqrt(dens);%/sqrt(sqrt(1+4*a1^2/401));
alpha=atan(kappa);
frcoef=0;
phi0=pi/2;
sp=0;

for n=2:N
    
    disp(n)
   t(n)=(n-1)*dt;
   %ey=f(0, t(n), tstart1, x(n-1));
   %ey1=a1*exp(-(t(n)-x(n-1)-tstart1)^2/2/tau^2)*sin(t(n)-x(n-1)-tstart1);
   ey1=2*a1/sqrt(1+kappa^2)*exp(-(t(n)-x(n-1)-tstart1)^2/2/tau^2)*cos(t(n)-x(n-1)+alpha);
   bz1=2*kappa*a1/sqrt(1+kappa^2)*exp(-(t(n)-x(n-1)-tstart1)^2/2/tau^2)*sin(t(n)-x(n-1)+alpha);
   ez1=2*a1/sqrt(1+kappa^2)*exp(-(t(n)-x(n-1)-tstart1)^2/2/tau^2)*cos(t(n)-x(n-1)+alpha+phi0);
   by1=-2*kappa*a1/sqrt(1+kappa^2)*exp(-(t(n)-x(n-1)-tstart1)^2/2/tau^2)*sin(t(n)-x(n-1)+alpha+phi0);
   
   ey=ey1;
   
   eyn(n)=ey;
   ezn(n)=0;
   bz=bz1;
   by=sp*by1;
   ex=dens*x(n-1)+frcoef*(px(n-1));
   ez=sp*ez1;
     
   pxminus=px(n-1)+q*dtdiv2*ex;
   pyminus=py(n-1)+q*dtdiv2*ey;
   pzminus=pz(n-1)+q*dtdiv2*ez;
   
   reciprocalgamma=1/sqrt(1+pxminus^2+pyminus^2+pzminus^2);
   
   ty=by*q*dtdiv2*reciprocalgamma;
   tz=bz*q*dtdiv2*reciprocalgamma;
   
   sy=2*ty/(1+ty^2+tz^2);
   sz=2*tz/(1+tz^2+ty^2);
   
   pxplus=(1-tz*sz-ty*sy)*pxminus+sz*pyminus-sy*pzminus;
   pyplus=-sz*pxminus+(1-tz*sz)*pyminus+sz*ty*pzminus;
   pzplus=sy*pxminus+sy*tz*pyminus+(1-sy*ty)*pzminus;
   
   px(n)=pxplus+q*dtdiv2*ex;
   py(n)=pyplus+q*dtdiv2*ey;
   pz(n)=pzplus+q*dtdiv2*ez;
   
   recgamma(n)=1/sqrt(1+px(n)^2+py(n)^2+pz(n)^2);
   gamma(n)=sqrt(1+px(n)^2+py(n)^2+pz(n)^2);
   
   x(n)=x(n-1)+recgamma(n)*px(n)*dt;
   y(n)=y(n-1)+recgamma(n)*py(n)*dt;
   z(n)=z(n-1)+recgamma(n)*pz(n)*dt;
 
   %erefl(n)=1*exp(-(t(n)-x(n)-tstart1)^2/2/tau^2)*cos(t(n)-x(n)+2*alpha);
%    
%     subplot(2,2,1)
  %   plot3(x(n), y(n), z(n), '.')
   %  hold on
     %axis([-0.1 0.1 -0.2 0.2])
% %    
%     subplot(2,2,2)
%     plot(x)
% %    
%     subplot(2,2,3)
%     if px(n)<0
%         pxplot=px(n);
%     else 
%         pxplot=0;
%     end
%     gammaplot=sqrt(1+pxplot^2);
% %    
%     plot(t(n), gammaplot)
%     hold on
% %        
%     subplot(2,2,4)
%     plot(erefl)
 %   pause(0.0001);
end

figure
plot(x/2/pi, y/2/pi)
xlabel('x, \mum', 'FontSize', 16)
ylabel('y, \mum', 'FontSize', 16)
set(gca, 'FontSize', 16)

figure
plot(t/2/pi, x/2/pi)
xlabel('t, periods', 'FontSize', 16)
ylabel('x, \mum', 'FontSize', 16)
set(gca, 'FontSize', 16)

figure
plot(t/2/pi, y/2/pi)
xlabel('t, periods', 'FontSize', 16)
ylabel('y, \mum', 'FontSize', 16)
set(gca, 'FontSize', 16)

figure
subplot(3,4,[1,2,3,5,6,7])
plot(x/2/pi, y/2/pi, 'k', 'LineWidth', 2)
xlabel('x, \mum', 'FontSize', 16)
ylabel('y, \mum', 'FontSize', 16)
set(gca, 'FontSize', 16)
axis([-0.5e-3 4e-3 -0.06 0.06])

subplot(3,4,[9 10 11])
plot(x/2/pi, t/2/pi, 'k', 'LineWidth', 2)
ylabel('t, periods', 'FontSize', 16)
xlabel('x, \mum', 'FontSize', 16)
set(gca, 'FontSize', 16)
axis([-0.5e-3 4e-3 2 10])


subplot(3,4,[4 8])
plot(t/2/pi, y/2/pi, 'k', 'LineWidth', 2)
xlabel('t, periods', 'FontSize', 16)
ylabel('y, \mum', 'FontSize', 16)
set(gca, 'FontSize', 16)
axis([2 10 -0.06 0.06])

%print -f4 -dpsc 'fig.ps'
%break

%plot(t/2/pi, gamma)
datasave(1,:)=t/2/pi;
datasave(2,:)=x/2/pi;
save('x-t.dat', 'datasave', 'ASCII')
%break

figure
plot(t, x/2/pi)
hold on

%break

%figure
%a=[0.1 1 2 3 4 5];
%pxmax=[1.7e-4 0.015 0.05 0.105 0.175 0.3]
%close all
%plot(a, pxmax, '*-')
%hold on
%plot(a, 0.012*a.^2)


%figure 
%plot(erefl)

%sp=fft(erefl);

%figure
%semilogy(abs(sp))

%figure
%plot(x(1:4000)/2/pi, y(1:4000)/2/pi)


%clear all
%a=[0 1 2 3 4 5 7 10 15 20];
%x=[0 5e-4 18e-4 3.5e-3 5.5e-3 7.5e-3 0.01 0.0161 0.023 0.026];
%y=[0 0.03 0.06 0.08 0.1 0.104 0.12 0.131 0.15 0.151];

%subplot(1,3,1)
%plot(a,x, '*--')
%subplot(1,3,2)
%plot(a, y, '*--')

%subplot(1,3,3)
%plot(py)