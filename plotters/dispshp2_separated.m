
function []=dispshp2_separated(undef,L,node,elem,mode,axesnum,scalem,m_a,BC,ifpatch)
%function []=dispshp2(L,node,elem,mode,axesnum,scalem,m_a,BC,ifpatch,Item3D,color4what_3D,Style3D,ifDataTip)

Item3D=1;    % can be 1,2,3 
color4what_3D=1; % can be 1:7
Style3D=1;  % can be 1 2
ifDataTip=1;  %can be 1 0

%BWS
%1998 originated
%2010 BWS and Z.Li improvements for speed and generality
%Dec 2022. S. Jin. 3D surface graphics objects are used for speed
%L=length
%node: [node# x z dofx dofz dofy dofrot stress] nnodes x 8;
%elem: [elem# nodei nodej t] nelems x 4;
%mode: vector of displacements global dof [u1 v1...un vn w1 01...wn 0n]'
%Item3D: the objects to be plotted: 1 - Deformed shape only, 2 - Undeformed shape only, 3 - Deformed shape & undeformed mesh
%color4what_3D: the value used to render the color of the surface.
	%1: Vector sum of Displacement, 2: X-Component of displacement, 3: Y-Component of displacement, 4: Z-Component of displacement
	%5: Y-Component of Normal Strain, 6: In-strip-plane of Shear Strain
	%7: No Color
%Style3D: surface or mesh

%确定横截面上节点的遍历路线 To determine the ergodic path(s) of the cross-section, "Euler path"
listNd=node(:,1);%节点编号列表 list of the nodes
listEdge=elem(:,1:3);%边的列表 list of the edges
nNd=size(listNd,1);%节点数量 number of nodes
nEdge=size(listEdge,1);%边的数量 number of edges
PathNodes={};%遍历所有边而不重复的“欧拉路径”问题 Node sequence on the path(s)
PathElems={};%路径遍历单元编号 element sequence on the path(s)
PathDirects={};%路径单元方向 is the element direction consistent to the path direction? 
while size(listEdge,1)>=1
	curPathNodes=[];
	curPathElems=[];
	curPathDirects=[];
	
	curPathNodes=listEdge(1,2:3);%任选一条边 just get one
	curPathElems=listEdge(1,1);
	curPathDirects=true;
	listEdge(1,:)=[];

	while 1
		%计算与路线终点相邻的剩余边的数量 number of the remainder edge(s) connect to the current end point of the path
		[rowB,colB]=find(listEdge(:,2:3)==curPathNodes(end),2);
		if length(rowB)==1%仅剩一条连接边 there is only one
			curPathNodes(end+1)=listEdge(rowB,4-colB);
			curPathElems(end+1)=listEdge(rowB,1);
			curPathDirects(end+1)=(colB==1);
			listEdge(rowB,:)=[];
			continue
		end
		%计算与路线起点相邻的剩余路径数量 number of the remainder edge(s) connect to the current start point of the path
		[rowA,colA]=find(listEdge(:,2:3)==curPathNodes(1),2);
		if length(rowA)==1%仅剩一条连接边 no other choice
			curPathNodes=horzcat(listEdge(rowA,4-colA),curPathNodes);
			curPathElems=[listEdge(rowA,1),curPathElems];
			curPathDirects=[(colA==2),curPathDirects];
			listEdge(rowA,:)=[];
			continue
		end
		if length(rowB)>=2%终端剩余两条及以上连接边 at least two remainder edges connect to the current end point of the path 
			curPathNodes(end+1)=listEdge(rowB(1),4-colB(1));
			curPathElems(end+1)=listEdge(rowB(1),1);
			curPathDirects(end+1)=(colB(1)==1);
			listEdge(rowB(1),:)=[];
			continue
		end
		if length(rowA)>=2%始端剩余两条及以上连接边 at least 2 edges to the start point
			curPathNodes=horzcat(listEdge(rowA(1),4-colA(1)),curPathNodes);
			curPathElems=[listEdge(rowA(1),1),curPathElems];
			curPathDirects=[(colA(1)==2),curPathDirects];
			listEdge(rowA(1),:)=[];
			continue
		end
		%so far both the end nodes of the current path has no edge connected to them
		PathNodes{end+1}=curPathNodes;
		PathElems{end+1}=curPathElems;
		PathDirects{end+1}=curPathDirects;
		break
	end
end

transSeg=4;%单元横向分段数量 trans segments of an element
longiSeg=4;%纵向每四分之一波分段数量
modeDisp=reshape(mode,4*nNd,[]);
%确定沿长度方向分段数量
%计算对应于各longitudinal terms 的变形向量的模
dispNorm_lgFuncs=vecnorm(modeDisp);
%不小于最大分量的1/50被认为足够显著
km=2*longiSeg*max(m_a(dispNorm_lgFuncs>0.02*max(dispNorm_lgFuncs)))+1;%不小于最大分量的1/50的最大半波数

nM_a=length(m_a);%纵向函数的数量 number of the longitudinal functions
nStrip=nEdge;%条单元数量同图问题中的Edge数量 strips are the edges

nd1=elem(:,2);
nd2=elem(:,3);
X1=node(nd1,2);X2=node(nd2,2);
Z1=node(nd1,3);Z2=node(nd2,3);
bStrip=sqrt((X2-X1).^2+(Z2-Z1).^2);
sinStrip=(Z2-Z1)./bStrip;
cosStrip=(X2-X1)./bStrip;

ptVect=(1:(transSeg-1))'/transSeg;
SubNodeX=[1-ptVect,ptVect]*[X1';X2'];%i-row j-column：X-coordinate of the i-th subnode on the j-th strip
SubNodeZ=[1-ptVect,ptVect]*[Z1';Z2'];%Z-coordinate


%计算变形量 displacement of the nodes
modeU=modeDisp(1:2:2*nNd,:);%the i-row j-column element: the U of the i-th node j-th function
modeV=modeDisp(2:2:2*nNd,:);
modeW=modeDisp(2*nNd+1:2:end,:);
modeR=modeDisp(2*nNd+2:2:end,:);

%计算单元横向等分点变形量 displacement of the subnodes of the strips
SubNodeUW=zeros((transSeg-1)*2,nM_a,nStrip); %the (2i-1)-th row,j-th column,k-th page element: the coefficient of the j-th function of the i-th equidistant subnode's U on the k-th strip

%the i-th row,j-th column,k-th page element: V, Equidistant Subnode i, Function j, Strip k
SubNodeV=zeros(transSeg-1,nM_a,nStrip); %i行j列k页元素：第k个单元第i个等分点的第j函数的系数

%逐单元计算等分点的变形量 calculating the displacement of the subnodes, loop on strips
for i=1:nStrip
	UWR2uwr=eye(6);
	UWR2uwr(1:2,1:2)=[cosStrip(i),sinStrip(i);-sinStrip(i),cosStrip(i)];
	UWR2uwr(4:5,4:5)=[cosStrip(i),sinStrip(i);-sinStrip(i),cosStrip(i)];
	
	uwrNd2uwSd=zeros(2*(transSeg-1),6);
	uwrNd2uwSd(1:2:2*(transSeg-1),1)=1-ptVect;
	uwrNd2uwSd(1:2:2*(transSeg-1),4)=ptVect;
	uwrNd2uwSd(2:2:2*(transSeg-1),2)=  ptVect.^3.*2-ptVect.^2.*3.+1;
	uwrNd2uwSd(2:2:2*(transSeg-1),3)=( ptVect.^3   -ptVect.^2.*2+ptVect).*bStrip(i);
	uwrNd2uwSd(2:2:2*(transSeg-1),5)= -ptVect.^3.*2+ptVect.^2.*3;
	uwrNd2uwSd(2:2:2*(transSeg-1),6)=( ptVect.^3   -ptVect.^2).*bStrip(i);
	
	uwSd2UWsd=zeros(2*(transSeg-1));
	for iSd=1:transSeg-1
		uwSd2UWsd(2*iSd-1:2*iSd,2*iSd-1:2*iSd)=[cosStrip(i),-sinStrip(i);sinStrip(i),cosStrip(i)];
	end
	SubNodeUW(:,:,i)=uwSd2UWsd*uwrNd2uwSd*UWR2uwr*[modeU(elem(i,2),:);modeW(elem(i,2),:);modeR(elem(i,2),:);modeU(elem(i,3),:);modeW(elem(i,3),:);modeR(elem(i,3),:)];
	
	SubNodeV(:,:,i)=[1-ptVect,ptVect]*[modeV(elem(i,2),:);modeV(elem(i,3),:)];
end

%确定函数在各分段点的取值
%number of cross-section：km
%total segments along the length：(km-1)
y2L=linspace(0,1,km);%y/L

if strcmp(BC,'S-S')
	funcUWR=sin(pi.*m_a'*y2L);
	funcV=cos(pi.*m_a'*y2L);
elseif strcmp(BC,'C-C')
	funcUWR=sin(pi.*m_a'*y2L).*sin(pi.*y2L);
	funcV=cos(pi.*m_a'*y2L).*sin(pi.*y2L)+sin(pi.*m_a'*y2L).*cos(pi.*y2L)./m_a';
elseif strcmp(BC,'S-C')||strcmp(BC,'C-S')
	funcUWR=sin(pi.*(m_a'+1)*y2L)+(m_a'+1)./m_a'.*sin(pi.*m_a'*y2L);
	funcV=(m_a'+1)./m_a'.*(cos(pi.*(m_a'+1)*y2L)+cos(pi.*m_a'*y2L));
elseif strcmp(BC,'F-C')||strcmp(BC,'C-F')
	funcUWR=1-cos(pi.*(m_a'-0.5)*y2L);
	funcV=(m_a'-0.5)./m_a'.*sin(pi.*(m_a'-0.5)*y2L);
elseif strcmp(BC,'G-C')||strcmp(BC,'C-G')
	funcUWR=sin(pi.*(m_a'-0.5)*y2L).*sin(pi/2.*y2L);
	funcV=(m_a'-.5)./m_a'.*cos(pi.*(m_a'-.5)*y2L).*sin(pi/2.*y2L)+sin(pi.*(m_a'-.5)*y2L).*cos(pi/2.*y2L)./m_a'./2;
else
	fprintf('\nError: Unrecognized boundary conditions.');
end

%Determine a scaling factor for the displaced shape
scale=scalem*1/3*max(max(node(:,2:3)));

watchon;
%
%
%
%axes(axesnum);
figure
hold on
cla 
axis off
hold on
colormap jet
nPath=length(PathNodes);


%绘制edge
if Item3D==3
	plot3([node(:,2)';node(:,2)'],[zeros(1,nNd);ones(1,nNd)*L],[node(:,3)';node(:,3)'],'Color',[.7 .7 .7],'LineStyle',':');
	for iPath=1:nPath
		plot3(node(PathNodes{iPath},2),zeros(length(PathNodes{iPath}),1) ,node(PathNodes{iPath},3),'Color',[.7 .7 .7],'LineStyle',':');
		plot3(node(PathNodes{iPath},2),ones(length(PathNodes{iPath}),1)*L,node(PathNodes{iPath},3),'Color',[.7 .7 .7],'LineStyle',':');
	end
end

%Generally, the color4data might be continuously distributed in the cross-section, or might be not
%When it is continuously distributed, the surface can be constructed according to the "path", which include several strips.
%When it is discontinuously across nodes, each surface should contain only one strip.
if color4what_3D==6 %绘制剪应变数据 for the in-strip-plan shear strain
	if strcmp(BC,'S-S')
		d_funcUWR=pi/L.*m_a'.*cos(pi.*m_a'*y2L);
	elseif strcmp(BC,'C-C')
		d_funcUWR=pi/L.*m_a'.*cos(pi.*m_a'*y2L).*sin(pi.*y2L)+pi/L.*sin(pi.*m_a'*y2L).*cos(pi.*y2L);
	elseif strcmp(BC,'S-C')||strcmp(BC,'C-S')
		d_funcUWR=pi/L.*(m_a'+1).*(cos(pi.*(m_a'+1)*y2L)+cos(pi.*m_a'*y2L));
	elseif strcmp(BC,'F-C')||strcmp(BC,'C-F')
		d_funcUWR=pi/L.*(m_a'-0.5).*sin(pi.*(m_a'-0.5)*y2L);
	elseif strcmp(BC,'G-C')||strcmp(BC,'C-G')
		d_funcUWR=pi/L.*(m_a'-.5).*cos(pi.*(m_a'-.5)*y2L).*sin(pi/2.*y2L)+pi*.5/L.*sin(pi.*(m_a'-.5)*y2L).*cos(pi/2.*y2L);
	else
		fprintf('\nError: Unrecognized boundary conditions.');
	end
	%逐 strip 计算剪应变
	shearStrainStripMesh=zeros(km,transSeg+1,nStrip);
	sVec=linspace(0,1,transSeg+1)';
	V1=modeV(nd1,:);%各板条节点1的V值列表
	V2=modeV(nd2,:);
	U1=modeU(nd1,:);
	U2=modeU(nd2,:);
	W1=modeW(nd1,:);
	W2=modeW(nd2,:);
	u1=[diag(cosStrip),diag(sinStrip)]*[U1;W1];
	u2=[diag(cosStrip),diag(sinStrip)]*[U2;W2];
	dV2ds_All=(V2-V1)./bStrip;%各板条横截面倾斜角
	xStripMesh=zeros(km,transSeg+1,nStrip);
	zStripMesh=zeros(km,transSeg+1,nStrip);
	yStripMesh=zeros(km,transSeg+1,nStrip);
	VStripMesh=zeros(km,transSeg+1,nStrip);
	UStripMesh=zeros(km,transSeg+1,nStrip);
	WStripMesh=zeros(km,transSeg+1,nStrip);
	for iStrip=1:nStrip
		dV2ds_curStrip=repmat(dV2ds_All(iStrip,:),(transSeg+1),1);
		uMesh_curStrip=[1-sVec,sVec]*[u1(iStrip,:);u2(iStrip,:)];
		shearStrainStripMesh(:,:,iStrip)=funcV'*dV2ds_curStrip'+d_funcUWR'*uMesh_curStrip';
		xStripMesh(:,:,iStrip)=repmat([X1(iStrip),X2(iStrip)]*[1-sVec';sVec'],km,1);
		zStripMesh(:,:,iStrip)=repmat([Z1(iStrip),Z2(iStrip)]*[1-sVec';sVec'],km,1);
		yStripMesh(:,:,iStrip)=repmat(linspace(0,L,km)',1,transSeg+1);
		VStripMesh(:,:,iStrip)=funcV'*[V1(iStrip,:);V2(iStrip,:)]'*[1-sVec,sVec]';
		UStripMesh(:,:,iStrip)=funcUWR'*[modeU(nd1(iStrip),:)',SubNodeUW(1:2:end,:,iStrip)',modeU(nd2(iStrip),:)'];
		WStripMesh(:,:,iStrip)=funcUWR'*[modeW(nd1(iStrip),:)',SubNodeUW(2:2:end,:,iStrip)',modeW(nd2(iStrip),:)'];
	end
	%确定极值
	[maxItem_value,maxItem_idx]=max(shearStrainStripMesh,[],'all','linear');
	[minItem_value,minItem_idx]=min(shearStrainStripMesh,[],'all','linear');
	[maxItem_sec,maxItem_nd,maxItem_strip]=ind2sub(size(shearStrainStripMesh),maxItem_idx);
	[minItem_sec,minItem_nd,minItem_strip]=ind2sub(size(shearStrainStripMesh),minItem_idx);
	if Style3D==1
		if Item3D==1 || Item3D==3
			%生成骨干mesh
			isMajorSec=false(km,1);
			isMajorSec(1:longiSeg:end)=true;%每四分之一波长为一个骨干截面 one cross-section each 1/2 half-wave
			xStrMesh_Skeleton=xStripMesh;
			xStrMesh_Skeleton(~isMajorSec,2:transSeg,:)=nan;
			for iStrip=1:nStrip
				mesh(UStripMesh(:,:,iStrip).*scale+xStrMesh_Skeleton(:,:,iStrip),VStripMesh(:,:,iStrip).*scale+yStripMesh(:,:,iStrip),WStripMesh(:,:,iStrip).*scale+zStripMesh(:,:,iStrip),zeros(km,transSeg+1),'EdgeColor',[.6 .6 .6],'LineStyle',':','FaceColor','none');
				surf(UStripMesh(:,:,iStrip).*scale+xStripMesh(:,:,iStrip),VStripMesh(:,:,iStrip).*scale+yStripMesh(:,:,iStrip),WStripMesh(:,:,iStrip).*scale+zStripMesh(:,:,iStrip),shearStrainStripMesh(:,:,iStrip),'EdgeColor','none','FaceColor','interp','FaceAlpha',1.0,'FaceLighting','flat');
			end
		else %Undeformed shape
			for iStrip=1:nStrip
				surf(xStripMesh(:,:,iStrip),yStripMesh(:,:,iStrip),zStripMesh(:,:,iStrip),shearStrainStripMesh(:,:,iStrip),'EdgeColor','none','FaceColor','interp','FaceAlpha',1.0,'FaceLighting','flat');
			end
		end
	else
		%生成 Mesh 的 grid
		%isMajorSec=false(km,1);
		%isMajorSec(1:longiSeg/2:end)=true;%每八分之一波长为一个骨干截面
		xStrMesh_Skeleton=xStripMesh;
		%xStrMesh_Skeleton(~isMajorSec,2:2:end,:)=nan;
		for iStrip=1:nStrip
			if Item3D==1 || Item3D==3
				mesh(UStripMesh(:,:,iStrip).*scale+xStrMesh_Skeleton(:,:,iStrip),VStripMesh(:,:,iStrip).*scale+yStripMesh(:,:,iStrip),WStripMesh(:,:,iStrip).*scale+zStripMesh(:,:,iStrip),shearStrainStripMesh(:,:,iStrip),'EdgeColor','interp','FaceColor','none','FaceAlpha',1.0,'FaceLighting','flat');
			else
				mesh(xStrMesh_Skeleton(:,:,iStrip),yStripMesh(:,:,iStrip),zStripMesh(:,:,iStrip),shearStrainStripMesh(:,:,iStrip),'EdgeColor','interp','FaceColor','none','FaceAlpha',1.0,'FaceLighting','flat');
			end
		end
	end
	
	
	if ifDataTip
		if Item3D==1 || Item3D==3
			maxItem_handle=plot3(UStripMesh(maxItem_sec,maxItem_nd,maxItem_strip).*scale+xStripMesh(maxItem_sec,maxItem_nd,maxItem_strip),VStripMesh(maxItem_sec,maxItem_nd,maxItem_strip).*scale+yStripMesh(maxItem_sec,maxItem_nd,maxItem_strip),...
				WStripMesh(maxItem_sec,maxItem_nd,maxItem_strip).*scale+zStripMesh(maxItem_sec,maxItem_nd,maxItem_strip),'k.','markersize',1);
			minItem_handle=plot3(UStripMesh(minItem_sec,minItem_nd,minItem_strip).*scale+xStripMesh(minItem_sec,minItem_nd,minItem_strip),VStripMesh(minItem_sec,minItem_nd,minItem_strip).*scale+yStripMesh(minItem_sec,minItem_nd,minItem_strip),...
				WStripMesh(minItem_sec,minItem_nd,minItem_strip).*scale+zStripMesh(minItem_sec,minItem_nd,minItem_strip),'k.','markersize',1);
		else
			maxItem_handle=plot3(xStripMesh(maxItem_sec,maxItem_nd,maxItem_strip),yStripMesh(maxItem_sec,maxItem_nd,maxItem_strip),...
				zStripMesh(maxItem_sec,maxItem_nd,maxItem_strip),'k.','markersize',1);
			minItem_handle=plot3(xStripMesh(minItem_sec,minItem_nd,minItem_strip),yStripMesh(minItem_sec,minItem_nd,minItem_strip),...
				zStripMesh(minItem_sec,minItem_nd,minItem_strip),'k.','markersize',1);
		end
		maxItem_handle.DataTipTemplate.DataTipRows(1) = dataTipTextRow('MX',maxItem_value);
		maxItem_handle.DataTipTemplate.DataTipRows(2:end) =[];
		datatip(maxItem_handle);
		
		minItem_handle.DataTipTemplate.DataTipRows(1) = dataTipTextRow('MN',minItem_value);
		minItem_handle.DataTipTemplate.DataTipRows(2:end) =[];
		datatip(minItem_handle);
	end
	settings3dPloting;
	return;
end

%Data数组在Cross-section 上不连续的情况前已处理完毕。现在，所需绘制的Data都是连续的了。绘图工作将基于 path 进行
%Now, all the Data are continuously distributed across the cross-section

%未变形的mesh
pathNdX=cell(nPath,1);
pathNdZ=cell(nPath,1);
pathVexClas=cell(nPath,1);
Xmesh=cell(nPath,1);
Ymesh=cell(nPath,1);
Zmesh=cell(nPath,1);
for iPath=1:nPath
	pathNdX{iPath}=node(PathNodes{iPath},2);
	pathNdZ{iPath}=node(PathNodes{iPath},3);
	%Path各顶点属性  true-节点  false-再分点
	pathVexClas{iPath}=false(1,transSeg*length(PathElems{iPath})+1);
	pathVexClas{iPath}(1:transSeg:end)=true;
	
	pathSdX=SubNodeX(:,PathElems{iPath});
	pathSdX(:,~PathDirects{iPath})=flipud(pathSdX(:,~PathDirects{iPath}));
	pathSdX=reshape(pathSdX,[],1);
	
	pathSdZ=SubNodeZ(:,PathElems{iPath});
	pathSdZ(:,~PathDirects{iPath})=flipud(pathSdZ(:,~PathDirects{iPath}));
	pathSdZ=reshape(pathSdZ,[],1);
	
	Xmesh{iPath}=nan(km,transSeg*length(PathElems{iPath})+1);
	Ymesh{iPath}=nan(km,transSeg*length(PathElems{iPath})+1);
	[Xmesh{iPath}(:,pathVexClas{iPath}),Ymesh{iPath}(:,pathVexClas{iPath})]=meshgrid(pathNdX{iPath},linspace(0,L,km));
	[Xmesh{iPath}(:,~pathVexClas{iPath}),Ymesh{iPath}(:,~pathVexClas{iPath})]=meshgrid(pathSdX,linspace(0,L,km));
	
	Zmesh{iPath}=nan(km,transSeg*length(PathElems{iPath})+1);
	Zmesh{iPath}(:,pathVexClas{iPath})=repmat(pathNdZ{iPath}',km,1);
	Zmesh{iPath}(:,~pathVexClas{iPath})=repmat(pathSdZ',km,1);
end


%计算 mesh 网格上的变形量
dispUmesh=cell(nPath,1);
dispWmesh=cell(nPath,1);
dispVmesh=cell(nPath,1);
for iPath=1:nPath
	pathNdU=modeU(PathNodes{iPath},:);
	pathNdW=modeW(PathNodes{iPath},:);
	pathNdV=modeV(PathNodes{iPath},:);
	
	pathSdU=SubNodeUW(1:2:end,:,PathElems{iPath});
	pathSdU(:,:,~PathDirects{iPath})=flipud(pathSdU(:,:,~PathDirects{iPath}));
	pathSdU=reshape(permute(pathSdU,[1,3,2]),[],nM_a);
	
	pathSdW=SubNodeUW(2:2:end,:,PathElems{iPath});
	pathSdW(:,:,~PathDirects{iPath})=flipud(pathSdW(:,:,~PathDirects{iPath}));
	pathSdW=reshape(permute(pathSdW,[1,3,2]),[],nM_a);
	
	pathSdV=SubNodeV(:,:,PathElems{iPath});
	pathSdV(:,:,~PathDirects{iPath})=flipud(pathSdV(:,:,~PathDirects{iPath}));
	pathSdV=reshape(permute(pathSdV,[1,3,2]),[],nM_a);
	
	%变形mesh汇总
	dispUmesh{iPath}=nan(km,transSeg*length(PathElems{iPath})+1);
	dispUmesh{iPath}(:,pathVexClas{iPath})=funcUWR'*pathNdU';
	dispUmesh{iPath}(:,~pathVexClas{iPath})=funcUWR(:,:)'*pathSdU'; % for the equidistant nodes on the real cross-sections
	
	dispWmesh{iPath}=nan(km,transSeg*length(PathElems{iPath})+1);
	dispWmesh{iPath}(:,pathVexClas{iPath})=funcUWR'*pathNdW';
	dispWmesh{iPath}(:,~pathVexClas{iPath})=funcUWR(:,:)'*pathSdW';
	
	dispVmesh{iPath}=nan(km,transSeg*length(PathElems{iPath})+1);
	dispVmesh{iPath}(:,pathVexClas{iPath})=funcV'*pathNdV';
	dispVmesh{iPath}(:,~pathVexClas{iPath})=funcV(:,:)'*pathSdV';
end


isMajorSec=false(km,1);
isMajorSec(1:longiSeg:end)=true;%每四分之一波长为一个骨干截面
XMeshSkeleton=Xmesh;
for iPath=1:nPath
	XMeshSkeleton{iPath}(~isMajorSec,~pathVexClas{iPath})=nan;
end
if (Item3D==1 || Item3D==3) && Style3D ==1 % deformed surface 的 标志网格
	for iPath=1:nPath
		mesh(dispUmesh{iPath}.*scale+XMeshSkeleton{iPath},dispVmesh{iPath}.*scale+Ymesh{iPath},dispWmesh{iPath}.*scale+Zmesh{iPath},zeros(size(XMeshSkeleton{iPath})),'EdgeColor',[0.6 0.6 0.6],'LineStyle',':','FaceColor','none');
	end
end


if color4what_3D==7 %no color data
	for iPath=1:nPath
		if Item3D==1 || Item3D==3
			if Style3D==1
				surf(dispUmesh{iPath}.*scale+Xmesh{iPath},dispVmesh{iPath}.*scale+Ymesh{iPath},dispWmesh{iPath}.*scale+Zmesh{iPath},zeros(size(Xmesh{iPath})),'EdgeColor','none','FaceColor','w','FaceAlpha',1.0,'FaceLighting','flat');
			else
				surf(dispUmesh{iPath}.*scale+Xmesh{iPath},dispVmesh{iPath}.*scale+Ymesh{iPath},dispWmesh{iPath}.*scale+Zmesh{iPath},zeros(size(Xmesh{iPath})),'EdgeColor','k','FaceColor','non','FaceAlpha',1.0,'FaceLighting','flat');
			end
		else
			if Style3D==1
				surf(Xmesh{iPath},Ymesh{iPath},Zmesh{iPath},zeros(size(Xmesh{iPath})),'EdgeColor','none','FaceColor','w','FaceAlpha',1.0,'FaceLighting','flat');
				surf(XMeshSkeleton{iPath},Ymesh{iPath},Zmesh{iPath},zeros(size(Xmesh{iPath})),'EdgeColor',[0.6 0.6 0.6],'LineStyle',':','FaceColor','w','FaceAlpha',1.0,'FaceLighting','flat');
			else
				surf(Xmesh{iPath},Ymesh{iPath},Zmesh{iPath},zeros(size(Xmesh{iPath})),'EdgeColor','k','FaceColor','none','FaceAlpha',1.0,'FaceLighting','flat');
			end
		end
	end
	
	
	
% 	
% 	if Style3D==1
% 		if Item3D==1 || Item3D==3
% 			for iPath=1:nPath
% 				surf(dispUmesh{iPath}.*scale+Xmesh{iPath},dispVmesh{iPath}.*scale+Ymesh{iPath},dispWmesh{iPath}.*scale+Zmesh{iPath},'EdgeColor','k','LineStyle','none','FaceColor','white','FaceAlpha',1.0,'FaceLighting','flat');
% 			end
% 		else
% 			for iPath=1:nPath
% 				surf(Xmesh{iPath},Ymesh{iPath},Zmesh{iPath},'EdgeColor','k','LineStyle','none','FaceColor','white','FaceAlpha',1.0,'FaceLighting','flat');
% 			end
% 		end
% 	else
% 		if Item3D==1 || Item3D==3
% 			for iPath=1:nPath
% 				mesh(dispUmesh{iPath}.*scale+Xmesh{iPath},dispVmesh{iPath}.*scale+Ymesh{iPath},dispWmesh{iPath}.*scale+Zmesh{iPath},zeros(size(Xmesh{iPath})),'LineStyle','-','FaceColor','none');
% 			end
% 		else
% 			for iPath=1:nPath
% 				mesh(Xmesh{iPath},Ymesh{iPath},Zmesh{iPath},'EdgeColor','g','LineStyle','-','FaceColor','none');
% 			end
% 		end
% 	end
	settings3dPloting;
	return;
end

%% 至此， color 所需要表示的数据类型是 1 到 5 类
if color4what_3D==1
	colorItem=cell(nPath,1);
	for iPath=1:nPath
		colorItem{iPath}=sqrt(dispUmesh{iPath}.^2+dispVmesh{iPath}.^2+dispWmesh{iPath}.^2);
	end
elseif color4what_3D==2
	colorItem=dispUmesh;
elseif color4what_3D==3
	colorItem=dispVmesh;
elseif color4what_3D==4
	colorItem=dispWmesh;
else %Item3D==5, 横截面正应变
	if strcmp(BC,'S-S')
		d_funcV=-pi/L.*m_a'.*sin(pi.*m_a'*y2L);
	elseif strcmp(BC,'C-C')
		d_funcV=2*pi/L.*cos(pi.*m_a'*y2L).*cos(pi.*y2L)-pi/L.*(m_a'+1./m_a').*sin(pi.*m_a'*y2L).*sin(pi.*y2L);
	elseif strcmp(BC,'S-C')||strcmp(BC,'C-S')
		d_funcV=-pi/L.*(m_a'+1)./m_a'.*(sin(pi.*(m_a'+1)*y2L).*(m_a'+1)+sin(pi.*m_a'*y2L).*m_a');
	elseif strcmp(BC,'F-C')||strcmp(BC,'C-F')
		d_funcV=pi/L.*(m_a'-0.5)^2./m_a'.*cos(pi.*(m_a'-0.5)*y2L);
	elseif strcmp(BC,'G-C')||strcmp(BC,'C-G')
		d_funcV=pi/L.*(m_a'-.5)./m_a'.*cos(pi.*(m_a'-.5)*y2L).*cos(pi/2.*y2L) - pi/L.*((m_a'-.5).^2+.25)./m_a'.*sin(pi.*(m_a'-.5)*y2L).*sin(pi/2.*y2L);
	else
		fprintf('\nError: Unrecognized boundary conditions.');
	end
	normalStrainYMesh=cell(nPath,1);
	for iPath=1:nPath
		normalStrainYMesh{iPath}=nan(km,transSeg*length(PathElems{iPath})+1);
		normalStrainYMesh{iPath}(:,pathVexClas{iPath})=d_funcV'*modeV(PathNodes{iPath},:)';
		
		pathSdV=SubNodeV(:,:,PathElems{iPath});
		pathSdV(:,:,~PathDirects{iPath})=flipud(pathSdV(:,:,~PathDirects{iPath}));
		pathSdV=reshape(permute(pathSdV,[1,3,2]),[],nM_a);
		normalStrainYMesh{iPath}(:,~pathVexClas{iPath})=d_funcV'*pathSdV';
	end
	colorItem=normalStrainYMesh;
end
%确定最大值和最小值的 Path序号、数值、Mesh坐标
maxItem_value=-Inf;
minItem_value=Inf;
for iPath=1:nPath
	[maxItem_val_thisPath,maxItem_idx_thisPath]=max(colorItem{iPath},[],'all','linear');
	if maxItem_value < maxItem_val_thisPath
		maxItem_value=maxItem_val_thisPath;
		maxItem_path=iPath;
		[maxItem_sec,maxItem_nd]=ind2sub(size(colorItem{iPath}),maxItem_idx_thisPath);
	end
	[minItem_val_thisPath,minItem_idx_thisPath]=min(colorItem{iPath},[],'all','linear');
	if minItem_value > minItem_val_thisPath
		minItem_value=minItem_val_thisPath;
		minItem_path=iPath;
		[minItem_sec,minItem_nd]=ind2sub(size(colorItem{iPath}),minItem_idx_thisPath);
	end
end

% 如果是 Deformed only 或者 Deformd shape + undeformed mesh, 则在deformed surf 上按 item 进行 render
if Item3D==1 || Item3D==3
	for iPath=1:nPath
		if Style3D==1
			surf(dispUmesh{iPath}.*scale+Xmesh{iPath},dispVmesh{iPath}.*scale+Ymesh{iPath},dispWmesh{iPath}.*scale+Zmesh{iPath},colorItem{iPath},'EdgeColor','none','FaceColor','interp','FaceAlpha',1.0,'FaceLighting','flat');
		else
			surf(dispUmesh{iPath}.*scale+Xmesh{iPath},dispVmesh{iPath}.*scale+Ymesh{iPath},dispWmesh{iPath}.*scale+Zmesh{iPath},colorItem{iPath},'EdgeColor','interp','FaceColor','none','FaceAlpha',1.0,'FaceLighting','flat');
		end
	end
	if ifDataTip
		maxItem_handle=plot3(dispUmesh{maxItem_path}(maxItem_sec,maxItem_nd).*scale+Xmesh{maxItem_path}(maxItem_sec,maxItem_nd),dispVmesh{maxItem_path}(maxItem_sec,maxItem_nd).*scale+Ymesh{maxItem_path}(maxItem_sec,maxItem_nd),...
			dispWmesh{maxItem_path}(maxItem_sec,maxItem_nd).*scale+Zmesh{maxItem_path}(maxItem_sec,maxItem_nd),'k.','markersize',1);
		minItem_handle=plot3(dispUmesh{minItem_path}(minItem_sec,minItem_nd).*scale+Xmesh{minItem_path}(minItem_sec,minItem_nd),dispVmesh{minItem_path}(minItem_sec,minItem_nd).*scale+Ymesh{minItem_path}(minItem_sec,minItem_nd),...
			dispWmesh{minItem_path}(minItem_sec,minItem_nd).*scale+Zmesh{minItem_path}(minItem_sec,minItem_nd),'k.','markersize',1);
	end
else %是 Undeformed only, 在 Undeformed surf 上按 item 进行 render
	for iPath=1:nPath
		if Style3D==1
			surf(Xmesh{iPath},Ymesh{iPath},Zmesh{iPath},colorItem{iPath},'EdgeColor','none','FaceColor','interp','FaceAlpha',1.0,'FaceLighting','flat');
		else
			mesh(Xmesh{iPath},Ymesh{iPath},Zmesh{iPath},colorItem{iPath},'EdgeColor','interp','FaceColor','none','FaceAlpha',1.0,'FaceLighting','flat');
		end
	end
	if ifDataTip
		maxItem_handle=plot3(Xmesh{maxItem_path}(maxItem_sec,maxItem_nd),Ymesh{maxItem_path}(maxItem_sec,maxItem_nd),...
			Zmesh{maxItem_path}(maxItem_sec,maxItem_nd),'k.','markersize',1);
		minItem_handle=plot3(Xmesh{minItem_path}(minItem_sec,minItem_nd),Ymesh{minItem_path}(minItem_sec,minItem_nd),...
			Zmesh{minItem_path}(minItem_sec,minItem_nd),'k.','markersize',1);
	end
	
end
if ifDataTip
	maxItem_handle.DataTipTemplate.DataTipRows(1) = dataTipTextRow('MX',maxItem_value);
	maxItem_handle.DataTipTemplate.DataTipRows(2:end) =[];
	datatip(maxItem_handle);
	
	minItem_handle.DataTipTemplate.DataTipRows(1) = dataTipTextRow('MN',minItem_value);
	minItem_handle.DataTipTemplate.DataTipRows(2:end) =[];
	datatip(minItem_handle);
end
settings3dPloting;

end

function settings3dPloting()
light('position',[1,1,2],'style','infinite');
light('position',[-2,-1,0],'style','infinite');
lighting phong
material metal

view(37.5,30)
hZoom=zoom;
setAxes3DPanAndZoomStyle(hZoom,gca,'camera');
hold off
axis equal
watchoff
end