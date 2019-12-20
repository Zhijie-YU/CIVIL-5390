clear,clc
E = 2e11; H = 2.5e10; sigY0 = 2.5e8;
mp = [E, H, sigY0];
epN = 0;
stressN = 0;
strainN = 0;
epNL = epN;
stressNL = stressN;
dT = 0.00025;
num = 8:8:32;
flag = -1;
deltaS = [];
for i = 1:size(num, 2)
    flag = flag * (-1);
    for j = 1:num(i)
        deltaS = [deltaS, flag*dT];
    end
end
strainNL = strainN;
for i = 1:size(deltaS,2)
    strainN = strainN + deltaS(i);
    strainNL = [strainNL, strainN];
    [stressN, epN] = IsoHard1D(mp, deltaS(i), stressN, epN);
    stressNL = [stressNL, stressN];
    epNL = [epNL, epN];
end
plot(strainNL*100, stressNL/1e6);
grid on;
xlabel('strain (%)');
ylabel('stress (MPa)');
axis([-0.5,0.5,-500,500]);