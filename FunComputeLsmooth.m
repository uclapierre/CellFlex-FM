function [Lsmooth] = FunComputeLsmooth(IDX,nb_stacks,path)

IMA = IDX;
Lsmooth = cell(size(IDX));
maxi = max(IDX{1}(:));

for t = 1:nb_stacks
    L = IMA{t};
    Lsmooth{t} = zeros(size(L));
    tic
    for j = 1:maxi
    t
    j
    mask = (L == j);
    if sum(sum(sum(mask>0))) > 0
        %smooth mask fiber
        parfor ii = 1:size(L,3)
            maskop(:,:,ii) = imopen(mask(:,:,ii),strel('line',10,70)) + imopen(mask(:,:,ii),strel('line',10,90)) + imopen(mask(:,:,ii),strel('line',10,110)) + imopen(mask(:,:,ii),strel('line',10,-20)) + imopen(mask(:,:,ii),strel('line',10,0)) + imopen(mask(:,:,ii),strel('line',20,10));
        end
        smask = bwskel(logical(maskop),'MinBranchLength',20); 
        Lsmooth{t} = Lsmooth{t} + (Lsmooth{t}==0).*(imdilate(smask,strel('sphere',2))*j);
    end
    end
    Lsmooth{t}(~ismember(Lsmooth{t},1:maxi)) = 0; 
    toc
    for z = 1:size(Lsmooth{1},3)
        t
        z
        imwrite(uint16(Lsmooth{t}(:,:,z)),[path '\Lsmooth' num2str(t,'%04.f') '_' num2str(z,'%04.f') '.tif']);
    end
end

end