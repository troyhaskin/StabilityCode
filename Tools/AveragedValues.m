clc;
clear('all');
clf();

f     = @(u) [(u-20).*(u-1).*(u+15);u.*sin(pi*u).^2];
df    = @(u) [-295 + 3*(-4 + u).*u, sin(pi*u).*(2*pi*u.*cos(pi*u)+sin(pi*u))];
duNew = @(u) df(u)*f(u);
duSec = @(unm0,unm1) -mean(f(unm0)./(f(unm0)-f(unm1)).*(unm0-unm1));

u = linspace(-5,5,1E3);

subplot(2,1,1);
plot(u,f(u))
hold('on');
subplot(2,1,2);
plot(u,diff(f(u)))
hold('on');

NotDone = true;
uNew = 10;
uSec = 10;
plot([uNew,uNew],[-20,20],'r--',[uSec,uSec],[-20,20],'k-.');
axis([-5,5,-1500,1500]);
uOld = 1.01*uSec;

while NotDone
    uNew  = uNew + duNew(uNew);

    uSecN = uSec + duSec(uSec,uOld);
    uOld  = uSec;
    uSec  = uSecN;

    subplot(2,1,1)
    plot([uNew,uNew],[-20,20],'r--',[uSec,uSec],[-20,20],'k-.');
    axis([-5,5,-10,10]);
    subplot(2,1,2)
    plot([uNew,uNew],[-20,20],'r--',[uSec,uSec],[-20,20],'k-.');
    axis([-5,5,-10,10]);
end