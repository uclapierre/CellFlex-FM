function [IDX] = FunTrackingFibers3D(L3,MASK,nb_stacks,path)
IMA = L3;
tic
maxi = 0;
for t = 1:nb_stacks
    if max(IMA{t},[],'all')>maxi
        t
    end
    maxi = max(maxi,max(IMA{t},[],'all'));
end
a = maxi;
b = size(IMA,2)-1;
Cost = cell(1,b);
IDX = cell(1,size(IMA,2));
IDX{1} = IMA{1};
for t = 2:size(IMA,2)
    IDX{t} = zeros(size(IDX{1}));
end

for t = 1:size(IMA,2)-1
    t
    for i = 1:a
        i
        mask1 = find(imdilate((IDX{1} == i).*(MASK>0),strel('sphere',3)));
        for j = 1:max(IMA{t+1},[],'all')
            mask2 = find(IMA{t+1}==j);
            Cost{t}(i,j) = size(intersect(mask1,mask2),1);
        end
    end
    for i = 1:max(IMA{t+1},[],'all')
        [l,idx] = max(Cost{t}(:,i));
        IDX{t+1} = IDX{t+1} + double(IMA{t+1}==i) * idx;
    end
end

for t = 1:size(IMA,2)-1
    t
    for z = 1:size(IDX{1},3)
          imwrite(uint16(IDX{t}(:,:,z)),[path,'\Track' num2str(t,'%04.f') '_' num2str(z,'%04.f') '.tif']);
    end
end
toc

end


