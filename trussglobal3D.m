function [qun,Qun] = trussglobal3D(nel,nSDOF,neDOF)
% nel = number of elements;
% nSDOF = number of DOF's (total);
% neDOF = number of DOF's per element;
load coords.txt; load displacements.txt; load DOFcon.txt; load EA.txt;
load forces.txt; load ID.txt; load ncon.txt;

% Search through ID array and order/renumber the unconstrained DOF
n = 1;
for i = 1:nSDOF
    if ID(i) == 0;
        ID(i) = n;
        n = n+1;
    end
end;

sizek11=n-1; % determine the size of q_qun, Q_k

% Search through ID array again and this time order/renumber the constrained DOFs
for i = 1:nSDOF
    if ID(i)==0.5;
        ID(i)=n;
        n=n+1;
    end
end
ID

% Update the DOFcon matrix based on the new DOF numbering
for i = 1:neDOF
    j = 1:nel;
    p=DOFcon(i,j);
    DOFcon(i,j)=ID(p);
end
DOFcon;

% Reorder the forces and displacment vectors consistent with the newDOF numbering 
% such that the constrained and unconstrained DOF's are grouped together, respectively

% Store the original forces and displacements vectors for reordering
forcesold = forces;
for i = 1:nSDOF
    for j=1:nSDOF
        if ID(j)==i
            forces(i)=forcesold(j);
        end
    end
end
forces;

displacementsold = displacements;
for i = 1:nSDOF
    for j=1:nSDOF
        if ID(j)==i
            displacements(i)=displacementsold(j);
        end
    end
end
displacements;

KG = zeros(nSDOF,nSDOF);
for i = 1:nel
    E=EA(i,1);
    A=EA(i,2);
    nN=ncon(1,i);
    nF=ncon(2,i);
    x1=coords(nN,1);
    y1=coords(nN,2);
    z1=coords(nN,3);
    x2=coords(nF,1);
    y2=coords(nF,2);
    z2=coords(nF,3);
    kq=getkq(x1,y1,z1,x2,y2,z2,E,A);
    for j=1:neDOF;
        map(j)=DOFcon(j,i);
    end

    for i = 1:neDOF
        ii=map(i);
        for j=1:neDOF
            jj=map(j);
            KG(ii,jj)=KG(ii,jj)+kq(i,j);
        end
    end
end

%partition the displacement and load vector as well as the stiffness matrix
K11=KG(1:sizek11,1:sizek11);
K12=KG(1:sizek11,sizek11+1:nSDOF);
K21=KG(sizek11+1:nSDOF,1:sizek11);
K22=KG(sizek11+1:nSDOF,sizek11+1:nSDOF);

Qk=forces(1:sizek11); %known forces
qk=displacements(sizek11+1:nSDOF); %known displacements
qun = inv(K11) * (Qk - K12 * qk); %unknown displacements
Qun=K21*qun+K22*qk; %unknown forces

% Function that assembles kq
function [kq] = getkq(x1,y1,z1,x2,y2,z2,E,A)
L = ((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^0.5; %txt files have been configured to have consistent units of kN,m,m^2
stiff = E*A/L;

%% stiffness matrix for 3D. Set some rows and columns to 0 to emphasize axial forces over shear forces.
kp = [1 0 0 -1 0 0; 
      0 0 0 0 0 0; 
      0 0 0 0 0 0; 
      -1 0 0 1 0 0; 
      0 0 0 0 0 0; 
      0 0 0 0 0 0] * stiff;

% Set certain direction cosines to 1 to reduce complexity in calculating rotations
lx = (x2-x1)/L; 
ly = 1;
lz = 1;
mx = (y2-y1)/L; 
my = 1;
mz = 1;
nx = (z2-z1)/L; 
ny = 1;
nz = 1;

% Transformation matrix for 3D
bQP = [lx, mx, lz, 0, 0, 0; 
       mx, my, mz, 0, 0, 0; 
       nx, ny, nz, 0, 0, 0; 
       0, 0, 0, lx, ly, lz; 
       0, 0, 0, mx, my, mz; 
       0, 0, 0, nx, ny, nz];

kq=bQP*kp*bQP';