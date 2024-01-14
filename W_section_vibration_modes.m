clear all 
close all
clc

location=['C:\Users\cmkra\Dropbox (CEDJHU)\Thinwalled Group\Rajshri\mfiles\fsm_vibration_git']; %change this line to your own path
addpath([location,'\abaqusmaker']);
addpath([location,'\analysis']);
addpath([location,'\analysis\cFSM']);
addpath([location,'\cutwp']);
addpath([location,'\helpers']);
addpath([location,'\holehelper']);
addpath([location,'\interface']);
addpath([location,'\plotters']);

inp_file = 'cross_secn_details_Ibeam';
[span, d, bf, tf, tw] = feval(str2func(inp_file));

l=0;
lengths_reg=logspace(0,3,100);
lengths_lcr=[span];
lengths=[lengths_reg, lengths_lcr];
lengths=(sort(lengths))';
nlengths = length(lengths);
BC = 'S-S';

for i=1:length(lengths)
    m_all{i}=[1];
end

while l<nlengths
l = l+1;
a = lengths(l);
m_a = m_all{l};
totalm = length(m_a);
end

[m_all]=msort(m_all);

prop=[100 29000.00 29000.00 0.0 0.0 11346.15 7.5000e-07];
%


node=[1 0.00 0.00 1 1 1 1 0
    2 bf/4 0 1 1 1 1 0
    3 bf/2 0 1 1 1 1 0
    4 (3*bf/4) 0 1 1 1 1 0
    5 bf 0 1 1 1 1 0
    6 bf/2 ((d-tf)/4) 1 1 1 1 0
    7 bf/2 ((d-tf)/2) 1 1 1 1 0
    8 bf/2 (3*(d-tf)/4) 1 1 1 1 0
    9 0 (d-tf) 1 1 1 1 0
    10 bf/4 (d-tf) 1 1 1 1 0
    11 bf/2 (d-tf) 1 1 1 1 0
    12 (3*bf/4) (d-tf) 1 1 1 1 0
    13 bf (d-tf) 1 1 1 1 0];
%
%Elements
elem=[1 1 2 tf 100
    2 2 3 tf 100
    3 3 4 tf 100
    4 4 5 tf 100
    5 3 6 tw 100
    6 6 7 tw 100
    7 7 8 tw 100
    8 8 11 tw 100
    9 9 10 tf 100
    10 10 11 tf 100
    11 11 12 tf 100
    12 12 13 tf 100];

%-----------------------------------------------------------------
%------------TWEAKING MODEL USING OTHER CUFSM FEATURES------------
%-----------------------------------------------------------------
%Features available in the GUI may also be used in these batch programs,
%for instance the mesh in this default file is rather course, let's double
%the number of elements
%[node,elem]=doubler(node,elem);
%
%In this example the stresses (last column of nodes) are already defined,
%but it is common to use the properties page to define these values instead
%of entering in the nodal stresses. Right now this problem applies a
%reference bending moment, let's apply a reference compressive load using
%the subroutines normally used on the properties page of CUFSM
%
%first calculate the global properties
[A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,J,Xs,Ys,Cw,B1,B2,w] = cutwp_prop2(node(:,2:3),elem(:,2:4));
thetap=thetap*180/pi; %degrees...
Bx=NaN; By=NaN;
%
%second set the refernce stress
fy=50;
%
%third calculate the P and M associated with the reference stress
unsymmetric=0; %i.e. do a restrained bending calculation
%HISTORY
%June 2010, complete update to new Boundary conditions, Z. Li, B. Schafer
%
%INPUTS
%prop: [matnum Ex Ey vx vy G] 6 x nmats
%node: [node# x z dofx dofz dofy dofrot stress] nnodes x 8;
%elem: [elem# nodei nodej t matnum] nelems x 5;
%lengths: [L1 L2 L3...] 1 x nlengths; lengths to be analyzed
%could be half-wavelengths for signiture curve
%or physical lengths for general b.c.
%springs: [node# d.o.f. kspring kflag] where 1=x dir 2= z dir 3 = y dir 4 = q dir (twist) flag says if k is a foundation stiffness or a total stiffness
%constraints:: [node#e dofe coeff node#k dofk] e=dof to be eliminated k=kept dof dofe_node = coeff*dofk_nodek
%GBTcon: GBTcon.glob,GBTcon.dist, GBTcon.local, GBTcon.other vectors of 1's
%  and 0's referring to the inclusion (1) or exclusion of a given mode from the analysis,
%  GBTcon.ospace - choices of ST/O mode
%         1: ST basis
%         2: O space (null space of GDL) with respect to K
%         3: O space (null space of GDL) with respect to Kg
%         4: O space (null space of GDL) in vector sense
%   GBTcon.norm - code for normalization (if normalization is done at all)
%         0: no normalization, 
%         1: vector norm
%         2: strain energy norm
%         3: work norm
%  GBTcon.couple - coupled basis vs uncoupled basis for general B.C. especially for non-simply supported B.C.
%         1: uncoupled basis, the basis will be block diagonal
%         2: coupled basis, the basis is fully spanned
%  GBTcon.orth - natural basis vs modal basis
%         1: natural basis
%         2: modal basis, axial orthogonality
%         3: modal basis, load dependent orthogonality
%BC: ['S-S'] a string specifying boundary conditions to be analyzed:
%'S-S' simply-pimply supported boundary condition at loaded edges
%'C-C' clamped-clamped boundary condition at loaded edges
%'S-C' simply-clamped supported boundary condition at loaded edges
%'C-F' clamped-free supported boundary condition at loaded edges
%'C-G' clamped-guided supported boundary condition at loaded edges
%m_all: m_all{length#}=[longitudinal_num# ... longitudinal_num#],longitudinal terms m for all the lengths in cell notation
% each cell has a vector including the longitudinal terms for this length
%neigs - the number of eigenvalues to be determined at length (default=10)

%OUTPUTS
%curve: buckling curve (load factor) for each length
%curve{l} = [ length mode#1
%             length mode#2
%             ...    ...
%             length mode#]
%shapes = mode shapes for each length
%shapes{l} = mode, mode is a matrix, each column corresponds to a mode.

%FUNCTIONS CALLED IN THIS ROUTINE
% \analysis\addspring.m : add springs to K
% \analysis\assemble.m : asseemble global K, Kg
% \analysis\elemprop.m : element properties
% \analysis\kglocal.m : element kg matrix
% \analysis\klocal.m : element k matrix
% \analysis\trans.m : trasnform k, kg matrix
% \analysis\msort.m : clean up 0's, multiple longitudinal terms. or out-of-order terms
% \analysis\constr_BCFlag.m : determine flags for user constraints and internal (at node) B.C.'s
% \analysis\cFSM\base_column : cFSM base vectors (natural basis, ST)
% \analysis\cFSM\base_update.m' : cFSM base vectors with selected basis,
%                                 orthogonalization, and normalization
% \analysis\cFSM\constr_user.m : user defined contraints in cFSM style
% \analysis\cFSM\mode_select.m : selection of modes for constraint matrix R
springs=0;

%In this case we have no constraint
% if any, constraint should be defined in this format:
% constraints=[node#e DOFe coeff. noder#k DOFk
%               ...   ...   ...    ...    ... ];
constraints=0;

%set the eigenmode you want to output, optional
neigs=10; %GUI default is 20


%cFSM features to utilize constrained finite strip method. We turned off here.
%Use of the cFSM variables are outside the scope of this batchcode example (at least for now)
%however, the following template for cFSM is ready to use with modification.
nnodes = length(node(:,1));
ndof_m= 4*nnodes;
GBTcon.ospace=1;GBTcon.couple=1;GBTcon.orth=1;GBTcon.norm=1; %see strip for possible solutions
[elprop,m_node,m_elem,node_prop,nmno,ncno,nsno,ndm,nlm,DOFperm]=base_properties(node,elem);
ngm=4;nom=2*(length(node(:,1))-1);
GBTcon.local=zeros(1,nlm);
GBTcon.dist=zeros(1,ndm);
GBTcon.glob=zeros(1,ngm); %can be less than 4 for special sections

% GBTcon.glob=ones(1,ngm); %can be less than 4 for special sections
% GBTcon.glob(1,1)=0;
% GBTcon.glob(1,3)=0;
% GBTcon.glob(1,4)=0;
GBTcon.other=zeros(1,nom);
% GBTcon.glob(3)=1;
%to turn on any cFSM classification use 1's on the GBTcon vectors defined
%above. Leave off unless you want it.
%for example say you wanted distortional only->GBTcon.dist=ones(1,ndm);


%-----------------------------------------------------------------
%---------------RUN THE ANALYSIS----------------------------------
%-----------------------------------------------------------------
%
[curve,shapes]=strip_main(prop,node,elem,lengths,springs,constraints,GBTcon,BC,m_all,neigs);
%

clas = 0;
%clas=classify(prop,node,elem,lengths,shapes,GBTcon,BC,m_all);

curvecell{1}=curve; %GUI expects to get a cell...,
filenamecell{1}=['Batch CUFSM5'];
clascell{1}=clas;
filedisplay=1;
minopt=1; %show min
logopt=1; %semilogx
clasopt=0; %classification stuff off
axescurve=figure(1);
xmin=min(lengths)*10/11;
xmax=max(lengths)*11/10;
%%
modeindex=1;

ymin=0;
    for j=1:max(size(curve));
        curve_sign(j,1)=curve{j}(modeindex,1);
        curve_sign(j,2)=curve{j}(modeindex,2);
    end

    ymax=min([max(curve_sign(:,2)),3*median(curve_sign(:,2))]);
modedisplay=1;
fileindex=1;
%lengthindex=find(lengths(:,1)==span);
lengthindex=90;
%lengthindex=find(diff(curve_sign(:,2))>0,1,'first');
picpoint=[curve{lengthindex}(modeindex,1) curve{lengthindex}(modeindex,2)];
%call the plotter
thecurve3(curvecell,filenamecell,clascell,filedisplay,minopt,logopt,clasopt,axescurve,xmin,xmax,ymin,ymax,modedisplay,fileindex,modeindex,picpoint)
grid on
pbaspect([4 3 1])

%2D buckled shape at a selected index into lengths
%inputs required for buckled shape in 2D
undef=1; %plot undeformed
%node
%elem
mode=shapes{lengthindex}(:,modeindex);
axes2dshapelarge=figure(2);
scale=1;
%springs
m_a=1;
%BC
SurfPos=1/2;
dispshap(undef,node,elem,mode,axes2dshapelarge,scale,springs,m_a,BC,SurfPos);
title(['Length (in.)=',num2str(lengths(lengthindex)),' f_n (Hz)=',num2str(curve{lengthindex}(modeindex,2))]);


for i = 1:nlengths
    vibration_length_freq(i,1) = curve{1,i}(modeindex,1);
    vibration_length_freq(i,2) = (curve{1,i}(modeindex,2));
    vibration_length_freq(i,3) = 1/vibration_length_freq(i,2);
end
%%
undef=1;
ifpatch=1;
L=100;
dispshp2_separated(undef,L,node,elem,mode,axes2dshapelarge,scale,m_a,BC,ifpatch);
% filename = 'vib_freq_lts_W_secn_beam.xlsx';
% sheetName = sprintf('%s %s %s %d', BC, '_m=', '_vibration mode', modeindex);
% col_header = {'Length of member (in.)', 'frequency (rad/sec)', 'Time period (secs)'}; 
% xlswrite(filename, vibration_length_freq, sheetName, 'A2');
% xlswrite(filename, col_header, sheetName, 'A1');

% pathname=[location,'\inputs_cufsm_gui\']
% filename = 'vibration_modes_W_secn_S_S_1'
% saver_fastfloor(pathname,filename,prop,node,elem,lengths,curve,shapes,springs,constraints,GBTcon,clas,BC,m_all)