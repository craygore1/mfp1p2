clear
close all
clc

rng("default")

p1 = 1;
p2 = 0.75;
k = 8;
sz = 1024;

initmat = ones(sz);

pvec = 1;

for i = 1:k
    temp = [pvec,pvec];
    
    len = 2^i;
    left = 1:len/2;
    right = len/2+1:len;


    half = rand;
    if half < 0.5
        temp(left) = temp(left).*p1;
        temp(right) = temp(right).*p2;
    else
        temp(left) = temp(left).*p2;
        temp(right) = temp(right).*p1;
    end
    pvec = temp;
end

centerx = floor(sz/2); %setting ellipse axes
centery = floor(sz/2); 
radiusx = centerx;
radiusy = centery;
fullradius = radiusx+radiusy;
[imcols, imrows] = meshgrid(1:sz, 1:sz);
myellipse = (imrows - centery).^2 ./ radiusy^2 ...
    + (imcols - centerx).^2 ./ radiusx^2 <= 1;

theta = linspace(-pi,pi,2^k + 1);

xvec = -floor(sz/2):floor(sz/2);
yvec = -floor(sz/2):floor(sz/2);
[Xim,Yim] = meshgrid(xvec,yvec);

angmat = zeros(sz);

for x = 1:sz
    for y = 1:sz
        angmat(x,y) = atan2(Yim(x,y),Xim(x,y));
    end
end

angcell = cell(length(theta)-1,1);

for t = 1:length(theta)-1
    angcell{t} = angmat > theta(t) & angmat <= theta(t+1);
end

for j = 1:length(theta)-1
    angcell{j} = cell2mat(angcell(j)).*pvec(j);
end

pmat = sum(cat(3,angcell{:}),3);
pmat = pmat.*myellipse;

for u = 1:sz
    for v = 1:sz
        mydec = rand;
        if mydec <= pmat(u,v)
            initmat(u,v) = 0;
        end
    end
end

figure(1)
imshow(initmat)

h = 0.1;

q = -10:h:10;

Dqtheory = zeros(length(q),1);

a = (p1+p2)/2;
b = p1/p2;

% for currq = 1:length(Dqtheory)
%     if q(currq) == 1
%         Dqtheory(currq) = (log(2*(p1+p2)) - (p1*log(p1) + p2*log(p2)))/log(2);
%     else
%         Dqtheory(currq) = 1 + (log(p1^q(currq) + p2^q(currq)) - q(currq)*log(p1+p2))/((1-q(currq))*log(2));
%     end
% end

for currq = 1:length(Dqtheory)
    if q(currq) == 1
        Dqtheory(currq) = log2(b+1) - (b*log2(b))/(b+1);
    else
        Dqtheory(currq) = (log2(b^q(currq) + 1) - q(currq)*log2(b+1))/(1-q(currq));
    end
end

tauq = (q'-1).*Dqtheory;

alphatheory = zeros(length(Dqtheory),1);
alphatheory(1) = (tauq(2) - tauq(1))/h;
alphatheory(end) = (tauq(end) - tauq(end-1))/h;

for step = 2:length(alphatheory)-1
    alphatheory(step) = (tauq(step+1) - tauq(step-1))/(2*h);
end

ftheory = q'.*alphatheory - tauq;
%% Plots
[Dqrec,alpharec,falpharec] = mfrectanglebinarized(initmat,0,q,0);
[Dqpol,alphapol,falphapol] = mfradialellipse(initmat,0,q,0);
[Dqtheta,alphatheta,falphatheta] = mfthetacoordinate(initmat,0,q,0);

% alphapol = alphapol - 1;
% falphapol = falphapol - 1;
% 
% alpharec = alpharec - 1;
% falpharec = falpharec - 1;

figure(2)
plot(q,Dqpol,'--r',q,Dqtheory,'k')
box on
grid, grid minor
legend("Box Counting","Theoretical")
xlabel('$q$','Interpreter','latex')
ylabel('$D_q$','Interpreter','latex')
% ylim([0 2])
fontname(gcf,"Times")

figure(3)
hold on
plot(alphatheory,ftheory,'k',LineWidth=1.25)
plot(alpharec,falpharec,'--r',LineWidth=1.25)
plot(alphapol,falphapol,'--b',LineWidth=1.25)
plot(alphatheta,falphatheta,'--',LineWidth=1.25,Color="#EDB120")
box on
grid, grid minor
legend("Theoretical","Rectangular","Polar","Radial")
fontname(gcf,"Times")
xlabel('$\alpha$','Interpreter','latex')
ylabel('$f(\alpha)$','Interpreter','latex')
% xlim([0 2])
ylim([0 2.1])
hold off

figure(4)
hold on
plot(alphatheory,ftheory,'k',LineWidth=1.25)
box on
grid, grid minor
fontname(gcf,"Times")
xlabel('$\alpha$','Interpreter','latex')
ylabel('$f(\alpha)$','Interpreter','latex')
xlim([0.2 1.8])
ylim([0 1.1])
hold off

errorrec = rms(ftheory - falpharec);
errorpol = rms(ftheory - falphapol);
errortheta = rms(ftheory - falphatheta);

