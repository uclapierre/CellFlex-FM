function [L3,MASK] = FunSegFibers3D(vq,nb_stacks,path)

for t = 1:nb_stacks
     vqeq{t} = imadjustn(vq{t});
end

%% Segmentation of T0 (ref)
for t = 1
    t
    figure
    MM = max(vqeq{1},[],3);
    imagesc(MM)
    [xthresh,ythresh] = ginput(1);
    thresh = MM(round(ythresh),round(xthresh));
    vqmed = medfilt3(vqeq{t},[3 3 3]);
    vqthresh = vqmed > thresh;
    L2 = bwlabeln(bwareaopen(vqthresh,1000));

%Erase false cell_body remaining
    figure
    imagesc(max(L2,[],3))
    hold on
    w = waitforbuttonpress;
    while w == 0
        r = drawfreehand;
        area = createMask(r);
        for z = 1:size(L2,3)
            Area(:,:,z) = area;
        end
        L2 = bwlabeln(bwareaopen(L2-L2.*Area,1000));
        figure
        imagesc(max(L2,[],3))
        hold on
        w = waitforbuttonpress;
    end
end

z = size(L2,3);
R = regionprops3(L2);
Centr = [R.Centroid];
idx = find(Centr(:,3)>floor(z/2));
idx2 = find(Centr(:,3)<floor(z/2));
L2top0 = bwlabeln(ismember(L2,idx));
L2bot0 = bwlabeln(ismember(L2,idx2));

for i = 1:max(L2bot0,[],'all')
    s0bot(i) = sum(sum(sum(L2bot0==i))); %volume of each labelled object
end

for i = 1:max(L2top0,[],'all')
    s0top(i) = sum(sum(sum(L2top0==i))); %volume of each labelled object
end

%% Masks of fibers endpoints
MASK1 = zeros(size(L2bot0));
MASK2 = zeros(size(L2top0));
Mask = {};
Mask2 = {};
dilMask = {};
dilMask2 = {};

figure
imagesc(max(L2bot0,[],3))
hold on

r1 = drawfreehand;
area1 = createMask(r1);
for z = 1:size(L2,3)
    Area1(:,:,z) = area1;
end

r2 = drawfreehand;
area2 = createMask(r2);
for z = 1:size(L2,3)
    Area2(:,:,z) = area2;
end

for i = 1:max(L2bot0,[],'all')
    i
    mask = (L2bot0==i).*Area1;
    mask2 = (L2bot0==i).*Area2;
    MASK1 = MASK1 + imdilate(mask,strel('sphere',8))*i;
    MASK2 = MASK2 + imdilate(mask2,strel('sphere',8))*i;
    Mask{i} = mask;
    Mask2{i} = mask2;
    dilMask{i} = imdilate(mask,strel('sphere',8));
    dilMask2{i} = imdilate(mask2,strel('sphere',8));

end
dilMask = dilMask(~cellfun(@isempty, Mask));
dilMask2 = dilMask2(~cellfun(@isempty, Mask2));
Mask = Mask(~cellfun(@isempty, Mask));
Mask2 = Mask2(~cellfun(@isempty, Mask2));


MASK3 = zeros(size(L2bot0));
MASK4 = zeros(size(L2top0));
Mask3 = {};
Mask4 = {};
dilMask3 = {};
dilMask4 = {};

figure
imagesc(max(L2top0,[],3))
hold on

r1 = drawfreehand;
area1 = createMask(r1);
for z = 1:size(L2,3)
    Area1(:,:,z) = area1;
end

r2 = drawfreehand;
area2 = createMask(r2);
for z = 1:size(L2,3)
    Area2(:,:,z) = area2;
end

for i = 1:max(L2top0,[],'all')
    i        
    mask3 = (L2top0==i).*Area1;
    mask4 = (L2top0==i).*Area2;
    MASK3 = MASK3 + imdilate(mask3,strel('sphere',8))*i;
    MASK4 = MASK4 + imdilate(mask4,strel('sphere',8))*i;
    Mask3{i} = mask3;
    Mask4{i} = mask4;
    dilMask3{i} = imdilate(mask3,strel('sphere',8));
    dilMask4{i} = imdilate(mask4,strel('sphere',8));

end
dilMask3 = dilMask3(~cellfun(@isempty, Mask3));
dilMask4 = dilMask4(~cellfun(@isempty, Mask4));
Mask3 = Mask3(~cellfun(@isempty, Mask3));
Mask4 = Mask4(~cellfun(@isempty, Mask4));
MASK = MASK1 + MASK2 + MASK3 + MASK4;

%% Segmentation of other timepoints
tic
for t = 1:nb_stacks
    t
    vqmed = medfilt3(vqeq{t},[3 3 3]);
    vqthresh = vqmed > thresh;
    L2 = bwlabeln(bwareaopen(vqthresh,1000));
    R = regionprops3(L2);
    Centr = [R.Centroid];
    idx = find(Centr(:,3)>floor(z/2));
    idx2 = find(Centr(:,3)<floor(z/2));
    L2top = bwlabeln(ismember(L2,idx));
    L2bot = bwlabeln(ismember(L2,idx2));

    BWbot = L2bot >0;
    BWbot2 = BWbot;
    
    for i = 1:max(L2bot,[],'all')
        s(i) = sum(sum(sum(L2bot==i))); %volume of each labelled object
    end
    for i = 1:max(L2bot,[],'all')
        if s(i) < 0.5*median(s0bot)
            BWbot(L2bot == i) = 0;
            BWbot2(L2bot == i) = 0;
        end
        if s(i) >1.3*median(s0bot) %objects that are too big (fused objects) are treated hereafter
            BWbot(L2bot == i) = 0;
            BWbot2(L2bot == i) = 0;
            for j = 1:size(Mask,2)
                if (sum(sum(sum(dilMask{j}.*(L2bot==i)))) > 0) && (sum(sum(sum(dilMask2{j}.*(L2bot==i)))) > 0)
                    j
                    end1 = (L2bot == i) .* dilMask{j};
                    end2 = (L2bot == i) .* dilMask2{j};
                    G = bwdistgeodesic(L2bot>0,logical(end1),'quasi-euclidean');
                    G2 = bwdistgeodesic(L2bot>0,logical(end2),'quasi-euclidean');
                    Gsum = G + G2;
                    Gsum = round(Gsum,2);
                    Gsum(isnan(Gsum))= 10000;
                    P = imregionalmin(Gsum);
                    P = double(P);
                    BWbot2 = BWbot2 + P + end1 + end2;
                end
            end
            for j = 1:size(Mask3,2)
                if (sum(sum(sum(dilMask3{j}.*(L2bot==i)))) > 0) && (sum(sum(sum(dilMask4{j}.*(L2bot==i)))) > 0)
                    j
                    end1 = (L2bot == i) .* dilMask3{j};
                    end2 = (L2bot == i) .* dilMask4{j};
                    G = bwdistgeodesic(L2bot>0,logical(end1),'quasi-euclidean');
                    G2 = bwdistgeodesic(L2bot>0,logical(end2),'quasi-euclidean');
                    Gsum = G + G2;
                    Gsum = round(Gsum,2);
                    Gsum(isnan(Gsum))= 10000;
                    P = imregionalmin(Gsum);
                    P = double(P);
                    BWbot2 = BWbot2 + P + end1 + end2;
                end
            end
        end
    end
    
    L3bot1 = imdilate(bwlabeln(BWbot2 - BWbot).*(L2bot>0),strel('sphere',2)).*(L2bot>0);
    L3bot2 = bwlabeln(BWbot);
    L3bot2(L3bot2>0) = L3bot2(L3bot2>0) + max(L3bot1(:));
    L3bot = L3bot1 + L3bot2;

    BWtop = L2top >0;
    BWtop2 = BWtop;

    if exist("s0top")
        for i = 1:max(L2top,[],'all')
            s(i) = sum(sum(sum(L2top==i))); %volume of each labelled object
        end
        for i = 1:max(L2top,[],'all')
            if s(i) <0.5*median(s0top)

                BWtop(L2top == i) = 0;
                BWtop2(L2top == i) = 0;
            end
            if s(i) >1.3*median(s0top)%objects that are too big (fused objects) are treated hereafter

                BWtop(L2top == i) = 0;
                BWtop2(L2top == i) = 0;
                for j = 1:size(Mask,2)
                    if (sum(sum(sum(dilMask{j}.*(L2top==i)))) > 0) && (sum(sum(sum(dilMask2{j}.*(L2top==i)))) > 0)
                        j
                        end1 = (L2top == i) .* dilMask{j};
                        end2 = (L2top == i) .* dilMask2{j};
                        G = bwdistgeodesic(L2top>0,logical(end1),'quasi-euclidean');
                        G2 = bwdistgeodesic(L2top>0,logical(end2),'quasi-euclidean');
                        Gsum = G + G2;
                        Gsum = round(Gsum,2);
                        Gsum(isnan(Gsum))= 10000;
                        P = imregionalmin(Gsum);
                        P = double(P);
                        BWtop2 = BWbot2 + P + end1 + end2;
                    end
                end
                for j = 1:size(Mask3,2)
                    if (sum(sum(sum(dilMask3{j}.*(L2top==i)))) > 0) && (sum(sum(sum(dilMask4{j}.*(L2top==i)))) > 0)
                        j
                        end1 = (L2top == i) .* dilMask3{j};
                        end2 = (L2top == i) .* dilMask4{j};
                        G = bwdistgeodesic(L2top>0,logical(end1),'quasi-euclidean');
                        G2 = bwdistgeodesic(L2top>0,logical(end2),'quasi-euclidean');
                        Gsum = G + G2;
                        Gsum = round(Gsum,2);
                        Gsum(isnan(Gsum))= 10000;
                        P = imregionalmin(Gsum);
                        P = double(P);
                        BWtop2 = double(logical(BWtop2 + P + end1 + end2));
                    end
                end
            end
        end
    end
    
    L3tot = bwlabeln(BWtop2 + BWbot2);
    L3{t} = L3tot;

    for z = 1:size(vqmed,3)
        imwrite(uint16(L3tot(:,:,z)),[path,'\Seg_' num2str(t,'%04.f') '_' num2str(z,'%04.f') '.tif']);

    end
       
end
toc
end

