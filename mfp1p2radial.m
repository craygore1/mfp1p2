clear
close all
clc

rng("default")

p1 = 1;
p2 = 0.75;
k = 10; % only use even k

sz = 2048;
sqsz = 2.^(1:k);
sqsz = repelem(sqsz,2);
sqsz(k+1:end) = [];

initmat = zeros(sz);
probmat = ones(sz);

% probmat(1:end/2,:) = p1;
% probmat((end/2)+1:end,:) = p2;

rowvec = sz;
colvec = sz;

counter = 0;

for x = sqsz
    counter = counter + 1;

    if mod(counter,2) == 0
        pxsz = floor(sz./x); 

        numboxC = floor(sz / pxsz); % vertical cut
        colvec = pxsz * ones(1, numboxC);

        probcell = mat2cell(probmat,rowvec,colvec);

        for i = 1:length(probcell(:,1))
            for j = 1:2:length(probcell(1,:))-1
                left = cell2mat(probcell(i,j));
                right = cell2mat(probcell(i,j+1));

                half = rand;
                if half < 0.5
                    left = left.*p1;
                    right = right.*p2;
                else
                    left = left.*p2;
                    right = right.*p1;
                end

                probcell{i,j} = left;
                probcell{i,j+1} = right;

                
            end
        end

        probmat = cell2mat(probcell);

    else
        pxsz = floor(sz./x);

        numboxR = floor(sz / pxsz); % horizontal cut
        rowvec = pxsz * ones(1, numboxR);

        probcell = mat2cell(probmat,rowvec,colvec);

        for i = 1:2:length(probcell(:,1))-1
            for j = 1:length(probcell(1,:))
                top = cell2mat(probcell(i,j));
                bot = cell2mat(probcell(i+1,j));

                half = rand;
                if half < 0.5
                    top = top.*p1;
                    bot = bot.*p2;
                else
                    top = top.*p2;
                    bot = bot.*p1;
                end

                probcell{i,j} = top;
                probcell{i+1,j} = bot;
            end
        end

        probmat = cell2mat(probcell);

    end
end

bandaid = zeros(length(probcell(:,1)),length(probcell(1,:)));

for ii = 1:length(probcell(:,1))
    for jj = 1:length(probcell(1,:))
        bandaid(ii,jj) = mean(mean(probcell{ii,jj}));
    end
end

xvec = round(linspace(-sz/2,sz/2,sz));
yvec = round(linspace(-sz/2,sz/2,sz));
[Xim,Yim] = meshgrid(xvec,yvec);

[theta,rho] = cart2pol(Xim,Yim); %polarcoordinates of all points

d1 = pi.*ones(2^(k/2),1);
d2 = -pi.*ones((2^(k/2))-1,1);
A = diag(d1) + diag(d2,-1);

areas = (pi/(2^(k/2))).*ones(2^(k/2),1);
rhorange = [0; sqrt(A\areas)];

rhorange = rescale(rhorange,0,sz/2);
thetarange = linspace(-pi,pi,(2^(k/2))+1);

for i = 1:length(rhorange)-1
    for j = 1:length(thetarange)-1
        temp1 = (rho < rhorange(i+1)) & (rho >= rhorange(i));
        temp2 = (theta < thetarange(j+1)) & (theta >= thetarange(j));
        temp3 = bandaid(i,j).*(temp1&temp2);
%         imagesc(temp3)
%         pause(0.1)
        initmat = temp3+initmat;
    end
end

mypic = ones(sz);

for len1 = 1:sz
    for len2 = 1:sz
        dec = rand;
        if initmat(len1,len2) < dec
            mypic(len1,len2) = 0;
        end
    end
end

mypic = ~mypic;
figure
imshow(mypic)
%% Multifractal analysis
h = 0.1;

q = -10:h:10;

Dqtheory = zeros(length(q),1);

a = (p1+p2)/2;
b = p1/p2;

for currq = 1:length(Dqtheory)
    if q(currq) == 1
        Dqtheory(currq) = 2*log2(b+1) - (2*b*log2(b))/(b+1);
    else
        Dqtheory(currq) = (2*log2(b^q(currq) + 1) - 2*q(currq)*log2(b+1))/(1-q(currq));
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

objcolor = 0;

[Dqrec,alpharec,falpharec] = mfrectanglebinarized(mypic,objcolor,q,0);
[Dqpol,alphapol,falphapol] = mfradialellipse(mypic,objcolor,q,0);
[Dqtheta,alphatheta,falphatheta] = mfthetacoordinate(mypic,objcolor,q,0);

% falphatheta = falphatheta + 1;
% alphatheta = alphatheta + 1;

% falphatheta = falphatheta*2;
% alphatheta = alphatheta*2;

%% Plots

figure(2)
plot(q,Dqtheory,'k',q,Dqrec,'--r',q,Dqpol,'--b')
box on
grid, grid minor
legend("Theoretical","Rectangular","Radial")
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
plot(alphatheta,falphatheta,'--r',LineWidth=1.25)
box on
grid, grid minor
fontname(gcf,"Times")
xlabel('$\alpha$','Interpreter','latex')
ylabel('$f(\alpha)$','Interpreter','latex')
xlim([0.85 1.15])
ylim([0.8 1.05])
hold off

% figure(4)
% hold on
% plot(errorrad,'--b',LineWidth=1.25)
% box on
% grid, grid minor
% hold off

errorrec = rms(ftheory - falpharec);
errorpol = rms(ftheory - falphapol);
errortheta = rms(ftheory - falphatheta);
            
        









