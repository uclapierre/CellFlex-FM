%% Section 1 : Initialization
clear all

Z = 21;                 % Number of stacks
xy_pix_size = 0.3e-6;   % xy pixel size (m)
z_pix_size = 1e-6;      % z pixel size (m)
L = 80e-6;              % Fiber length (m)
nb_nodes = 41;          % Number of nodes for the FE model
n_elem = nb_nodes-1;    % Number of elements for the FE model
param_fiber = 600;      % Exposure time (us) of the fibers (600,650,700,750,800)

[vq,nb_stacks,path] = Read4Dimage(Z,xy_pix_size,z_pix_size); % Loading of the 3D+T image stack with resampling for isometric final voxel size

% Creation of the subfolders
cd(path)
mkdir Seg
mkdir Track
mkdir Def3D
mkdir Cell_Mask
mkdir Force3D
mkdir Lsmooth
pathSeg = [path,'\Seg'];
pathTrack = [path,'\Track'];
pathDef3D = [path,'\Def3D'];
pathCellMask = [path,'\Cell_Mask'];
pathForce3D = [path,'\Force3D'];
pathLsmooth = [path,'\Lsmooth'];
pathInputCell = [path,'\cell'];
pathInputMask = [path,'\cell mask'];

%% Section 2 : Segmentation and Tracking
[L3,MASK] = FunSegFibers3D(vq,nb_stacks,pathSeg);
IDX = FunTrackingFibers3D(L3,MASK,nb_stacks,pathTrack);
save("L3 MASK IDX","L3","MASK","IDX", '-v7.3')
Lsmooth = FunComputeLsmooth(IDX,nb_stacks,pathLsmooth);

%% Section 3 : Intersection and Deflection
vqCell = Read4Dimage(Z,xy_pix_size,z_pix_size);
Mask_mode = 0;
[Inter,Mask_cell] = Calculate_Intersection(Mask_mode,pathCellMask,pathInputMask,vqCell,nb_stacks,Lsmooth,Z,xy_pix_size,z_pix_size); % Calculate intersection between cell mask and fibers
save("Inter MaskCell","Inter","Mask_cell", '-v7.3')

[Deflection,Intersection] = FunCalculateDeflection(Lsmooth,Inter,nb_stacks,Z,xy_pix_size,nb_nodes);
save("Deflection Lsmooth Intersection","Deflection","Lsmooth","Intersection", '-v7.3')

[DEF3D] = FunDeflectionAmplitude3D(pathDef3D,Deflection,Lsmooth,n_elem,nb_nodes);
save("DEF3D","DEF3D", '-v7.3')

%% Section 4 : Traction force recovery
%Inverse x and y column when needed and subsample to 41 nodes
Deflection2 = {};
Intersection2 = Intersection;
for t = 1:nb_stacks
    for j = 1:size(Deflection,1)
        if size(Deflection{j,t})>0
            j
            t
            deflection = Deflection{j,t};
            SUM = sum(deflection,1);
            if abs(SUM(1))>abs(SUM(2))
                deflection(:,2) = deflection(:,1);
                deflection(:,1) = 0;
            else
            end
            Deflection2{j,t} = deflection;
        end
    end
end

%Regularization Parameters
lambda1 = 100;
lambda2 = 1e9;

%Inverse traction recovery
pres = FunConstrainedOptimization(L,param_fiber,Deflection2,Intersection2,nb_nodes,nb_stacks,lambda1,lambda2,1e-13);

[FORCE3D] = FunForceAmplitude3D(pathForce3D,nb_stacks,pres,Lsmooth,n_elem,nb_nodes);
save("pres FORCE3D","FORCE3D","pres", '-v7.3')

