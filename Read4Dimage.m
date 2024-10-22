function [vq,nb_stacks,path] = Read4Dimage(Z,pixel_xy,pixel_z,path)

if nargin<4
path = uigetdir;
end
directory = dir([path,'\*.tif']);
nb_stacks=length(directory)/Z;
tmpIM = imread([directory(1).folder '\' directory(1).name]);

[l m] = size(tmpIM);

IM = zeros(l,m,Z);

for t = 1:nb_stacks
    t
    for z = 1:Z
        IM(:,:,z) = imread([directory((t-1)*(Z)+z).folder '\' directory((t-1)*(Z)+z).name]);
    end   
    IMA{t} = IM;
end

mini = min(IM(:));
maxi = max(IM(:));
if isa(tmpIM,'single')
IM = uint16(2^16*(IM-mini/(maxi-mini)));
elseif isa(tmpIM,'uint8')
IM = uint8(IM);
elseif isa(tmpIM,'uint16')
IM = uint16(IM);
end

ratio = pixel_xy / pixel_z;
for t = 1:nb_stacks
    t
    F = griddedInterpolant(double(IMA{t}),'linear');

    [sx,sy,sz] = size(IM);
    xq = (1:sx)';
    yq = (1:sy)';
    zq = (1:ratio:sz)';
    vq{t} = uint16(F({xq,yq,zq}));
end
end

