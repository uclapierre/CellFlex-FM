function [Inter,Mask_cell] = Calculate_Intersection(bool,path,pathInputMask,vq,nb_stacks,Lsmooth,Z,xy_pix_size,z_pix_size)

if bool == 1
    for t = 1:nb_stacks
        vqeq{t} = imadjustn(vq{t});
        vqmed{t} = medfilt3(vqeq{t});
    end
    
    figure
    MM = max(vqmed{40},[],3);
    imagesc(MM)
    [xthresh,ythresh] = ginput(1);
    thresh = MM(round(ythresh),round(xthresh));
    for t = 1:nb_stacks
        t
        vqthresh = vqmed{t} > thresh;
        vqmed2 = medfilt3(vqthresh,[3,3,3]);
        Mask_cell{t} = vqmed2;
        Ltmp = Lsmooth{t};
        Maskopen = imopen(Mask_cell{t},strel("sphere",1))>1;
        Inter{t} = (Ltmp>0).*Maskopen;
    
        M = Maskopen;
        for z = 1:size(M,3)
              imwrite(uint16(M(:,:,z)),[path '\CellMask' num2str(t,'%04.f') '_' num2str(z,'%04.f') '.tif']);
        end
    end

elseif bool == 0
    Mask_cell = Read4Dimage(Z,xy_pix_size,z_pix_size,pathInputMask);
    for t = 1:nb_stacks
        Ltmp = Lsmooth{t};
        Maskopen = imopen(Mask_cell{t},strel("sphere",2))>1;
        Inter{t} = (Ltmp>0).*(Maskopen>0);
        M = Maskopen;
        for z = 1:size(M,3)
            imwrite(uint16(M(:,:,z)),[path '\CellMask' num2str(t,'%04.f') '_' num2str(z,'%04.f') '.tif']);
        end
    end

end
