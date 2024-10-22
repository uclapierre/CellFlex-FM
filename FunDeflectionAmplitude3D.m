function [DEF3D] = FunDeflectionAmplitude3D(path,Deflection,Lsmooth,n_elem,n_nodes)

nb_stacks = size(Lsmooth,2);
DEF3D = cell(size(Lsmooth));
for k = 1:nb_stacks
    k
    DEF3D{1,k} = zeros(size(Lsmooth{1}));
        for j = nonzeros(unique(Lsmooth{k}(:)))' 
            if size(Deflection{j,k},1)>0
            mask = zeros(size(Lsmooth{k}));
            M = Lsmooth{k}==j;
            R = regionprops3(M,'Centroid','VoxelList');
            Centr = [R.Centroid];
            [Ix] = R.VoxelList{1}(:,2);
            Ix = sort(unique(Ix));
            [Iy] = R.VoxelList{1}(:,1);
            Iy = sort(unique(Iy));
           % if Centr(:,3)>floor(size(M,3)/2)
           if size(Ix,1)>size(Iy,1)
               j
               "top"
               ratio = n_elem/size(Ix,1);
               xq = 1:ratio:n_nodes;

               if size(Deflection{j,k},1)>0
                   DEF = Deflection{j,k};
                   DEF2 = DEF(:,1:3);
                   SumDEF2 = smooth(sum(abs(DEF2),2));
                   INTERSumDEF2 = griddedInterpolant(SumDEF2,'linear');
                   SUMDEF2Q = INTERSumDEF2({xq});
                   Deflection2{j,k} = SUMDEF2Q;
               end

               for N = 1:size(Ix,1)                   
                   mask(Ix(N),:,:) = Deflection2{j,k}(N).*(M(Ix(N),:,:)>0);
                end
            else
                j
                "bot"
                ratio = n_elem/size(Iy,1);
                yq = 1:ratio:n_nodes;

                        if size(Deflection{j,k},1)>0
                            DEF = Deflection{j,k};
                            DEF2 = DEF(:,1:3);
                            SumDEF2 = smooth(sum(abs(DEF2),2));
                            INTERSumDEF2 = griddedInterpolant(SumDEF2,'linear');
                            SUMDEF2Q = INTERSumDEF2({yq});
                            Deflection2{j,k} = SUMDEF2Q;
                        end
                    
                for N = 1:size(Iy,1)
                   mask(:,Iy(N),:) = Deflection2{j,k}(N).*(M(:,Iy(N),:)>0);
                end
            end
            DEF3D{1,k} = DEF3D{1,k} + mask;
            end
        end
end

for t = 1:nb_stacks
    t
    mini = min(DEF3D{t}(:));
    maxi = 50e-6;
    IM = uint16(2^16*(DEF3D{t}-mini)/(maxi-mini));
    for z = 1:size(DEF3D{1},3)
          imwrite(uint16(IM(:,:,z)),[path '\Def' num2str(t,'%04.f') '_' num2str(z,'%04.f') '.tif']);
    end
end
end


