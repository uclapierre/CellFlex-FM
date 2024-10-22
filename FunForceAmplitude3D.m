function [FORCE3D] = FunForceAmplitude3D(path,nb_stacks,pres,Lsmooth,n_elem,n_nodes)

FORCE3D = cell(size(Lsmooth));
for k = 1:nb_stacks
    k
    FORCE3D{1,k} = zeros(size(Lsmooth{1}));
        for j = 1:max(Lsmooth{k}(:))
            if size(pres{j,k},1)>0
            mask = zeros(size(Lsmooth{k}));
            M = Lsmooth{k}==j;
            R = regionprops3(M,'Centroid','VoxelList');
            Centr = [R.Centroid];
            [Ix] = R.VoxelList{1}(:,2);
            Ix = sort(unique(Ix));
            [Iy] = R.VoxelList{1}(:,1);
            Iy = sort(unique(Iy));
           if size(Ix,1)>size(Iy,1)
               j
               "top"
               ratio = n_elem/size(Ix,1);
               xq = 1:ratio:n_nodes;
               if size(pres{j,k},1)>0
                   FORCE = pres{j,k};
                   FORCE2 = FORCE(:,1:3);
                   SumFORCE2 = smooth(sum(abs(FORCE2),2));
                   INTERSumFORCE2 = griddedInterpolant(SumFORCE2,'linear');

                   SUMFORCE2Q = INTERSumFORCE2({xq});
                   pres2{j,k} = SUMFORCE2Q;
               end

               for N = 1:size(Ix,1)
                   mask(Ix(N),:,:) = pres2{j,k}(N).*(M(Ix(N),:,:)>0);
                end
            else
                j
                "bot"
                
                ratio = n_elem/size(Iy,1);
                yq = 1:ratio:n_nodes;

                        if size(pres{j,k},1)>0
                            FORCE = pres{j,k};
                            FORCE2 = FORCE(:,1:3);
                            SumFORCE2 = smooth(sum(abs(FORCE2),2));
                            INTERSumFORCE2 = griddedInterpolant(SumFORCE2,'linear');

                            SUMFORCE2Q = INTERSumFORCE2({yq});
                            pres2{j,k} = SUMFORCE2Q;
                        end
                                   
                for N = 1:size(Iy,1)
                   mask(:,Iy(N),:) = pres2{j,k}(N).*(M(:,Iy(N),:)>0);
                end
            end
            FORCE3D{1,k} = FORCE3D{1,k} + mask;
            end
        end

end

for t = 1:nb_stacks
    t
    mini = 0;
    maxi = 100e-9;
    IM = uint16(2^16*(FORCE3D{t}-mini)/(maxi-mini));
    for z = 1:size(FORCE3D{1},3)
          imwrite(uint16(IM(:,:,z)),[path '\FORCE' num2str(t,'%04.f') '_' num2str(z,'%04.f') '.tif']);
    end
end
end


