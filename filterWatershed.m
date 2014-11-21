%Watershed analysis of fibers

%import image
% path = 'M:\Prosjekter\2014 Fiber Image Analysis\data\';
% img = zeros(1024,1024,39);
% for i = 0:38
%     istr = num2str(i);
%     if i < 10
%         istr = ['0' istr];
%     end
%     filename = ['rype.lif_Series058_z' istr '.tif'];
%     tempim = imread([path filename]);
%     img(:,:,i+1) = tempim(:,:,2);
% end

%
%imagesc(img(:,:,20));

% m3DIN = img(550:600,550:600,1:38);
% vout = savitzkyGolay3D_rle_coupling(51,51,38,m3DIN,7,7,7,5);
% figure, imagesc(m3DIN(:,:,10));
% figure, imagesc(vout(:,:,10));

gx = 1:size(m3DIN,1);
gy = 1:size(m3DIN,2);
gz = 1:size(m3DIN,3);
[X,Y,Z] = meshgrid(gx,gy,gz);
%m3DINs = smooth3(m3DIN);
fv = isosurface(X,Y,Z,vout,150);
p = patch(fv, 'FaceColor', 'red', 'EdgeColor', 'none');
view(3), axis tight
camlight 
lighting gouraud
colormap copper

%3D gradient operatro
gx = [-1 -2 -1;0 0 0;1 2 1;];
gy = [-1 0 1;-2 0 2;-1 0 1;];
gy = [-1 0 1;-2 0 2;-1 0 1;];