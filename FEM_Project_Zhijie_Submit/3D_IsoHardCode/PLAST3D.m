function PLAST3D(PROP, ETAN, UPDATE, LTAN, NE, NDOF, XYZ, LE)
%% Compute global tangent stiffness matrix and residual force 
% ETAN: elastic stiffness matrix
% UPDATA(Boolean): whether to update SIGMA and XQ
% LTAN(Boolean): whether to compute and update global tangent stiffness matrix GKF
% Based on nodal delta displacement(DISPDD), compute and update SIGMA(UPDATE), XQ(UPDATE), GKF(LTAN), Force 
  global DISPDD DISPTD FORCE GKF SIGMA XQ
  %
  % Integration points and weights (2-point integration)
  XG=[-0.57735026918963D0, 0.57735026918963D0];
  WGT=[1.00000000000000D0, 1.00000000000000D0];
  %
  % Index for history variables (each integration pt)
  INTN=0;
  %
  %LOOP OVER ELEMENTS, THIS IS MAIN LOOP TO COMPUTE K AND F
  for IE=1:NE
    % Nodal coordinates and incremental displacements
    % ELXY: (num of nodes in 1 element * 3)
    ELXY=XYZ(LE(IE,:),:);
    % Local to global mapping
    % eg. local node 3 to global node 9 to global (25,26,27) saved in IDOF
    IDOF=zeros(1,24);
    for I=1:8
      II=(I-1)*NDOF+1;
      IDOF(II:II+2)=(LE(IE,I)-1)*NDOF+1:(LE(IE,I)-1)*NDOF+3;
    end
    DSP=DISPTD(IDOF); % (1,24)
    DSPD=DISPDD(IDOF);
    DSP=reshape(DSP,NDOF,8); % (3,8)
    DSPD=reshape(DSPD,NDOF,8);
    %
    %LOOP OVER INTEGRATION POINTS
    for LX=1:2, for LY=1:2, for LZ=1:2
      % followings are all for a specific Gauss point
      % [E1, E2, E3] coordinates of a point in parent domain
      E1=XG(LX); E2=XG(LY); E3=XG(LZ);
      INTN = INTN + 1;
      % Determinant and shape function derivatives
      % SHPD: gradient of shape functions in physical domain (B) (3*8)
      % det of Jacobian matrix at the given point
      [~, SHPD, DET] = SHAPEL([E1 E2 E3], ELXY);
      FAC=WGT(LX)*WGT(LY)*WGT(LZ)*DET;
      %
      % Previous converged history variables
      % SIGMA: stress at each integration point
      NALPHA=6; 
      STRESSN=SIGMA(1:6,INTN);
      
      ALPHAN=XQ(1:NALPHA,INTN);
      % effective plastic strain
      EPN=XQ(NALPHA+1,INTN);
      % DSPD: delta displacement at nodes of this element
      % DEPS: delta strain at GP (3*3)
      DEPS=DSPD*SHPD'; 
      % DDEPS: delta strain at GP (6*1)
      DDEPS=[DEPS(1,1) DEPS(2,2) DEPS(3,3) ...
             DEPS(1,2)+DEPS(2,1) DEPS(2,3)+DEPS(3,2) DEPS(1,3)+DEPS(3,1)]';

      % Computer stress, back stress & effective plastic strain

      % Infinitesimal plasticity
      [STRESS, ALPHA, EP]=combHard(PROP,ETAN,DDEPS,STRESSN,ALPHAN,EPN);
      % Update plastic variables
      if UPDATE
        SIGMA(1:6,INTN)=STRESS;
        XQ(:,INTN)= [ALPHA; EP];
        continue;
      end
      %
      % Add residual force and tangent stiffness matrix
      % BM:gradient of shape functions 
      BM=zeros(6,24);
      for I=1:8
        COL=(I-1)*3+1:(I-1)*3+3;
        BM(:,COL)=[SHPD(1,I) 0         0;
                   0         SHPD(2,I) 0;
                   0         0         SHPD(3,I);
                   SHPD(2,I) SHPD(1,I) 0;
                   0         SHPD(3,I) SHPD(2,I);
                   SHPD(3,I) 0         SHPD(1,I)];
      end
      % Residual forces
      FORCE(IDOF) = FORCE(IDOF) - FAC*BM'*STRESS;
      % Tangent stiffness
      if LTAN
        DTAN=combHardTan(PROP,ETAN,DDEPS,STRESSN,ALPHAN,EPN);
        EKF = BM'*DTAN*BM;
        GKF(IDOF,IDOF)=GKF(IDOF,IDOF)+FAC*EKF;
      end
    end,        end,        end
  end
end