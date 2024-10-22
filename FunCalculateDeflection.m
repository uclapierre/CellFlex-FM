function [Deflection,Intersection] = FunCalculateDeflection(Lsmooth,Inter,nb_stacks,Z,xy_pix_size,nb_nodes)

maxi = 0;
for t = 1:nb_stacks
    if max(Lsmooth{t},[],'all')>maxi
        t
    end
    maxi = max(maxi,max(Lsmooth{t},[],'all'));
end

Deflection = cell(maxi,nb_stacks);

for t = 1:nb_stacks
    L = Lsmooth{t};
    for j = 1:maxi
    t
    j
    tic 
    mask = (L== j);
    if sum(sum(sum(mask>0))) > 0
        %smooth mask fiber
        for ii = 1:size(L,3)
            maskop(:,:,ii) = imopen(mask(:,:,ii),strel('line',10,70)) + imopen(mask(:,:,ii),strel('line',10,90)) + imopen(mask(:,:,ii),strel('line',10,110)) + imopen(mask(:,:,ii),strel('line',10,-20)) + imopen(mask(:,:,ii),strel('line',10,0)) + imopen(mask(:,:,ii),strel('line',20,10));
        end
        smask = bwskel(logical(maskop),'MinBranchLength',20); 

        Intdil = imdilate(Inter{t},strel('sphere',3)).*imdilate(smask,strel('sphere',2));
        smask = double(smask);
        smask(smask.*Intdil==1)=2;

        [Ix,Iy,Iz] = ind2sub(size(L),find(smask));
        if size(Ix,1)>0
            if abs(Ix(end)-Ix(1)) > abs(Iy(end)-Iy(1)) %X oriented Fiber
                [Ix,I] = sort(Ix);
                Iy = Iy(I);
                Iz = Iz(I);
                mask3 = zeros(size(mask));
                mask4 = zeros(size(mask));
                mask3(Ix(1),Iy(1),Iz(1)) = 1;
                mask4(Ix(end),Iy(end),Iz(end)) = 1;
                d = round(bwdist(mask3) + bwdist(mask4));
                m = imregionalmin(d);
                ms = bwskel(m);
        
                [Jx,Jy,Jz] = ind2sub(size(L),find(ms));
                [Jx,J] = sort(Jx);
                Jy = Jy(J);
                Jz = Jz(J);
                [JJx,JJy,JJz] = ind2sub(size(L),find(smask));
                [JJx,JJ] = sort(JJx);
                JJy = JJy(JJ);
                JJz = JJz(JJ);
                step = size(Jx,1)/nb_nodes;
                k=1;
                deflection = zeros(floor(size(Jx,1)/step),3);
                
                parfor k = 1:nb_nodes
                    k
                    i = floor(1+(k-1)*step);
                    if i<=size(Jy,1) && i<=size(JJy,1) && i<=size(Jz,1) && i<=size(JJz,1)
                    deflection(k,:) = [0,JJy(i)-Jy(i),JJz(i)-Jz(i)];
                    intersection(k) = smask(JJx(i),JJy(i),JJz(i)) == 2;
                    else
                    deflection(k,:) = 0;
                    intersection(k) = 0;
                    end
                end
                               
            else %Y oriented Fiber
                [Iy,I] = sort(Iy); 
                Ix = Ix(I);
                Iz = Iz(I);
                mask3 = zeros(size(mask));
                mask4 = zeros(size(mask));
                mask3(Ix(1),Iy(1),Iz(1)) = 1;
                mask4(Ix(end),Iy(end),Iz(end)) = 1;
                d = round(bwdist(mask3) + bwdist(mask4));
                m = imregionalmin(d);
                ms = bwskel(m);
        
                [Jx,Jy,Jz] = ind2sub(size(L),find(ms));
                [Jy,J] = sort(Jy);
                Jx = Jx(J);
                Jz = Jz(J);
                [JJx,JJy,JJz] = ind2sub(size(L),find(smask));
                [JJy,JJ] = sort(JJy);
                JJx = JJx(JJ);
                JJz = JJz(JJ);
                step = size(Jy,1)/nb_nodes;
                k=1;
                deflection = zeros(floor(size(Jy,1)/step),3);
                               
                parfor k = 1:nb_nodes
                    k
                    i = floor(1+(k-1)*step);
                    if i<=size(Jx,1) && i<=size(JJx,1) && i<=size(Jz,1) && i<=size(JJz,1)
                    deflection(k,:) = [JJx(i)-Jx(i),0,JJz(i)-Jz(i)];
                    intersection(k) = smask(JJx(i),JJy(i),JJz(i)) == 2;
                    else
                    deflection(k,:) = 0;
                    intersection(k) = 0;
                    end
                end
                 
            end
            
        end
        Deflection{j,t} = deflection*xy_pix_size;
        Intersection{j,t} = intersection;
        
    end
    toc
    end
end

end


