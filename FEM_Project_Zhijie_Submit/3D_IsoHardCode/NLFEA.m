function NLFEA(ITRA,TOL,ATOL,NTOL,TIMS,NOUT,PROP,EXTFORCE,SDISPT,XYZ,LE)
%% This is the main program for a linear combined isotropic/kinematic hardening model
% beta = 0: purely isotropic; beta = 1: purely kinematic
%% Input variables:
% ITRA: maximum number of convergence iterations
% TOL: residual convergence limit
% ATOL: residual divergence limit
% NTOL: maximum number of bisections
% PROP: material
% EXTFORCE: natural boundary
% SDISPT: essential boundary (displacement)
% XYZ: coordinates
%% global variables:
% DISPDD: nodal displacement increment(num of nodes * 3, 1)
% DISPTD: nodal displacement (num of nodes * 3, 1)
% FORCE: nodal residual force (num of nodes * 3, 1)
% GKF: global tangent stiffness matrix (num of nodes * 3, num of nodes * 3)
% SIGMA: stress at each Gauss integration point (6, num of GP)=(6, num of GP for 1 element * num of elements)
% XQ: history variables(back stress alpha(6)+effective plastic strain(1)) at each GP (7, num of GP)

global DISPDD DISPTD FORCE GKF			%Global variables
%
[NUMNP, NDOF] = size(XYZ);				% Analysis parameters
NE = size(LE,1);
NEQ = NDOF*NUMNP;
%
DISPTD=zeros(NEQ,1); DISPDD=zeros(NEQ,1);	% Nodal displacement & increment
ETAN=PLSET(PROP, NE);	% Initialize material properties
%
ITGZONE(XYZ, LE, NOUT);					% Check element connectivity
%
% Load increments [Start End Increment InitialLoad FinalLoad]
NLOAD=size(TIMS,2);
ILOAD=1;						% First load increment
TIMEF=TIMS(1,ILOAD);				% Starting time
TIMEI=TIMS(2,ILOAD);				% Ending time
DELTA=TIMS(3,ILOAD);				% Time increment
CUR1=TIMS(4,ILOAD);					% Starting load factor
CUR2=TIMS(5,ILOAD);					% Ending load factor
DELTA0 = DELTA;					% Saved time increment
TIME = TIMEF;						% Starting time
TDELTA = TIMEI - TIMEF;				% Time interval for load step
ITOL = 1;						% Bisection level
TARY=zeros(NTOL,1);					% Time stamps for bisections
%
% Load increment loop
%----------------------------------------------------------------------
ISTEP = -1; % control whether to output the results of this substep
FLAG10 = 1;
while(FLAG10 == 1)	
  % control whether to end the whole program
  FLAG10 = 0; 
  FLAG11 = 1; 
  FLAG20 = 1;
  %
  CDISP = DISPTD; 					% Store converged displacement
  %
  if(ITOL==1) 					% No bisection
    DELTA = DELTA0;
    TARY(ITOL) = TIME + DELTA;
  else							% Recover previous bisection
    ITOL = ITOL-1;					% Reduce the bisection level
    DELTA = TARY(ITOL)-TARY(ITOL+1);		% New time increment
    TARY(ITOL+1) = 0;				% Empty converged bisection level
    ISTEP = ISTEP - 1;				% Decrease load increment
  end
  TIME0 = TIME;					% Save the current time
  %
  % Update stresses and history variables from the previous step
  % stress and history variables are not updated during NR iteration
  UPDATE=true; LTAN=false;
  PLAST3D(PROP, ETAN, UPDATE, LTAN, NE, NDOF, XYZ, LE); 
  % Print results
  if(ISTEP>=0)
      PROUT(NOUT, TIME, NUMNP, NE, NDOF)
  end
  %
  TIME = TIME + DELTA;				% Increase time
  ISTEP = ISTEP + 1;
  %
  % Check time and control bisection
  while(FLAG11 == 1)				% Bisection loop start
    FLAG11 = 0;
    if ((TIME-TIMEI)>1E-10)			% Time passed the end time
      if ((TIMEI+DELTA-TIME)>1E-10)		% One more at the end time
        DELTA=TIMEI+DELTA-TIME;			% Time increment to the end
        DELTA0=DELTA;				% Saved time increment
        TIME=TIMEI;					% Current time is the end
      else
        ILOAD=ILOAD+1;				% Progress to next load step
        if(ILOAD>NLOAD)				% Finished final load step
          FLAG10 = 0;				% Stop the program
          break;
        else						% Next load step
          TIME=TIME-DELTA;
          DELTA=TIMS(3,ILOAD);
          DELTA0=DELTA;
          TIME = TIME + DELTA;
          TIMEF = TIMS(1,ILOAD);
          TIMEI = TIMS(2,ILOAD);
          TDELTA = TIMEI - TIMEF;
          CUR1 = TIMS(4,ILOAD);
          CUR2 = TIMS(5,ILOAD);
        end
      end
    end
    %
    % Load factor and prescribed displacements
    FACTOR = CUR1 + (TIME-TIMEF)/TDELTA*(CUR2-CUR1);
    SDISP = DELTA*SDISPT(:,3)/TDELTA*(CUR2-CUR1);
    %
    % Start convergence iteration
    %------------------------------------------------------------------
    ITER = 0;
    DISPDD = zeros(NEQ,1);
    while(FLAG20 == 1)
      % Newton-Raphson iteration begins
      FLAG20 = 0;
      ITER = ITER + 1;
      %
      % Initialize global stiffness K and residual vector F
      GKF = sparse(NEQ,NEQ);
      FORCE = sparse(NEQ,1);
      %
      % Assemble K and F
      UPDATE=false; LTAN=true;
      PLAST3D(PROP, ETAN, UPDATE, LTAN, NE, NDOF, XYZ, LE); 

      % Increase external force
      if size(EXTFORCE,1)>0
        LOC = NDOF*(EXTFORCE(:,1)-1)+EXTFORCE(:,2);
        FORCE(LOC) = FORCE(LOC) + FACTOR*EXTFORCE(:,3);
      end
      %
      % Prescribed displacement BC
      NDISP=size(SDISPT,1);
      if NDISP~=0
       FIXEDDOF=NDOF*(SDISPT(:,1)-1)+SDISPT(:,2);
       GKF(FIXEDDOF,:)=zeros(NDISP,NEQ);
       % PROP(1) used as a very large value, we use penalty method to
       % solve the system
       GKF(FIXEDDOF,FIXEDDOF)=PROP(1)*eye(NDISP);
     
       FORCE(FIXEDDOF)=0;
       % initial disp
       if ITER==1, FORCE(FIXEDDOF) = PROP(1)*SDISP(:); end
      end
      % Check convergence
      if(ITER>1)
        FIXEDDOF=NDOF*(SDISPT(:,1)-1)+SDISPT(:,2);
        ALLDOF=1:NEQ;
        FREEDOF=setdiff(ALLDOF,FIXEDDOF);
        RESN=max(abs(FORCE(FREEDOF)));
        OUTPUT(1, ITER, RESN, TIME, DELTA)
        %
        if(RESN<TOL)
          FLAG10 = 1;
          break;
        end
        %
        if ((RESN>ATOL)||(ITER>=ITRA))			% Start bisection
          ITOL = ITOL + 1;
          if(ITOL<=NTOL)
            DELTA = 0.5*DELTA;
            TIME = TIME0 + DELTA;
            TARY(ITOL) = TIME;
            DISPTD=CDISP;
            fprintf(1,'Not converged. Bisecting load increment %3d\n',ITOL);
          else
            fprintf(1,'\t\t *** Max No. of bisection ***\n'); return;
          end
          FLAG11 = 1;
          FLAG20 = 1;
          break;
        end
      end
      %
      % Solve the system equation!!
      if(FLAG11 == 0)
        SOLN = GKF\FORCE;
        DISPDD = DISPDD + SOLN;
        DISPTD = DISPTD + SOLN;
        FLAG20 = 1;
      else
        FLAG20 = 0;
      end
      if(FLAG10 == 1), break; end
    end 							%20 Convergence iteration
  end 								%11 Bisection
end 								%10 Load increment
%
% Successful end of program
fprintf(1,'\t\t *** Successful end of program ***\n');
fprintf(NOUT,'\t\t *** Successful end of program ***\n');
end

function OUTPUT(FLG, ITER, RESN, TIME, DELTA)
%% Print convergence iteration history
  if FLG == 1
    if ITER>2
      fprintf(1,'%27d %14.5e \n',ITER,full(RESN));
    else
      fprintf(1,'\n \t Time  Time step   Iter \t  Residual \n');
      fprintf(1,'%10.5f %10.3e %5d %14.5e \n',TIME,DELTA,ITER,full(RESN));
    end
  end
end

function PROUT(NOUT, TIME, NUMNP, NE, NDOF)
%% Print converged displacements and stresses

  global SIGMA DISPTD
  %
  fprintf(NOUT,'\r\n\r\nTIME = %11.3e\r\n\r\nNodal Displacements\r\n',TIME);
  fprintf(NOUT,'\r\n Node          U1          U2          U3');
  for I=1:NUMNP
    II=NDOF*(I-1);
    fprintf(NOUT,'\r\n%5d %11.3e %11.3e %11.3e',I,DISPTD(II+1:II+3));
  end
  fprintf(NOUT,'\r\n\r\nElement Stress\r\n');
  fprintf(NOUT,'\r\n        S11         S22         S33         S12         S23         S13');
  for I=1:NE
    fprintf(NOUT,'\r\nElement %5d',I);
    II=(I-1)*8;
    fprintf(NOUT,'\r\n%11.3e %11.3e %11.3e %11.3e %11.3e %11.3e',SIGMA(1:6,II+1:II+8));
  end
  fprintf(NOUT,'\r\n\r\n');
  disp('Node 12 stress & disp')
  disp(SIGMA(3,1))
  disp(DISPTD(36))
end

function ETAN=PLSET(PROP, NE)
%% calculate elastic stiffness matrix and initialize stress and history
%% variables at GP
% Initialize history variables and elastic stiffness matrix
% XQ    : 1-6 = Back stress alpha, 7 = Effective plastic strain
% SIGMA : Stress at GP
% ETAN  : Elastic stiffness matrix

  global SIGMA XQ
  %
  LAM=PROP(1);
  MU=PROP(2);
  %
  N = 8*NE;

  SIGMA=zeros(6,N);
  XQ=zeros(7,N);
  ETAN=[LAM+2*MU LAM      LAM      0  0  0;
        LAM      LAM+2*MU LAM      0  0  0;
        LAM      LAM      LAM+2*MU 0  0  0;
        0        0        0        MU 0  0;
        0        0        0        0  MU 0;
        0        0        0        0  0  MU];

end

function ITGZONE(XYZ, LE, NOUT)
%% Check element connectivity to avoid negative Jacobian
%% the order of nodes in ELXY should be consistent with that in XNODE(in function SHAPEL)!!!
  EPS=1E-7;
  NE = size(LE,1);
  for I=1:NE
    ELXY=XYZ(LE(I,:),:);
    [~, ~, DET] = SHAPEL([0 0 0], ELXY);
    DVOL = 8*DET;
    if DVOL < EPS
      fprintf(NOUT,'\n??? Negative Jacobian ???\nElement connectivity\n');
      fprintf(NOUT,'%5d',LE(I,:));
      fprintf(NOUT,'\nNodal Coordinates\n');
      fprintf(NOUT,'%10.3e %10.3e %10.3e\n',ELXY');
      error('Negative Jacobian');
    end
  end
end