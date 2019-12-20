% 1 load step computation for 1D isotropic harding problem
function [stress, ep]=IsoHard1D(mp, deltaS, stressN, epN)
% mp(material parameters): [E,H,sigY0]
% deltaS(delta strain): strain increment
% stressN: current stress 
% epN: current plastic strain
% stress: updated stress after this load step
% ep: updated plastic strain

E = mp(1); H = mp(2); sigY0 = mp(3);
% trial solution
sigTr = stressN + deltaS * E;
sigYn = sigY0 + H * epN;
fTr = abs(sigTr) - sigYn;
% check state and update
if fTr <= 0
    stress = sigTr;
    ep = epN;
else
    deltaEp = fTr/(E + H);
    stress = sigTr - sign(sigTr)*E*deltaEp;
    ep = epN + deltaEp;
end
return


    

