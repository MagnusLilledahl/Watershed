function vOUT = savitzkyGolay3D_rle_coupling(Nx,Ny,Nz,m3DIN,lengthX,lengthY,lengthZ,order)
%
% input:
% Nx = int, (fastest)number of pixels in the 1st direction
% Ny = int, number of pixels in the 2nd direction
% Nz = int, (slowest)number of pixels in the 3rd direction
% m3DIN(x,y,z) = *IMPORTANT* Nx(row)*Ny(column)*Nz(layer) matrix, 3D input matrix for smoothing
%         real or complex numbers
%                    e.g. experimental data with noise
% lengthX = int, must be *odd*, the size of the window in the 1st direction
% lengthY = int, must be *odd*, the size of the window in the 2nd direction
% lengthZ = int, must be *odd*, the size of the window in the 3rd direction
% order    = int, order of polynomial
% 
% output:
% vOUT     = Nx(row)*Ny(column)*Nz(layer) matrix, 3D smoothed data of m3DIN
% 
% Reference: 
% [1] Abraham Savitzky and Marcel J. E. Golay, 'Smoothing and
% differentiation of data by simplified least sqaure procedures',
% Analytical Chemistry, Vol. 36, No. 8, page 1627-1639, July 1964
%
% Author: Shao Ying HUANG (shaoying.h@gmail.com)
% Date: 9 April 2012
%% 
% window size
Nloc = lengthX*lengthY*lengthZ;
% order
ord = order;
% moving window
x = -(lengthX-1)/2:1:(lengthX-1)/2;
y = -(lengthY-1)/2:1:(lengthY-1)/2;
z = -(lengthZ-1)/2:1:(lengthZ-1)/2;

coor=zeros(Nloc,3);
index3D=zeros(lengthX,lengthY,lengthZ);
count = 1;
for zz = 1:lengthZ
    for yy = 1:lengthY
        for xx = 1:lengthX
            coor(count,1) = x(xx);
            coor(count,2) = y(yy);
            coor(count,3) = z(zz);
            index3D(xx,yy,zz)=count;
            count = count +1;
        end
    end
end

for nn = 1:Nloc
    count = 1;
    for nx = 0:ord
        for ny = 0:ord
            for nz = 0:ord
                if nx+ny+nz<=ord
                    A(nn,count) = coor(nn,1)^nx*coor(nn,2)^ny*coor(nn,3)^nz;
                    count = count+1;
                end
            end
        end
    end
end

AT = A';
AT_A = AT*A;

F = eye(Nloc);    

g = zeros(Nloc,Nloc);
for nn = 1:Nloc %excitation
    CC = AT_A\(AT*F(:,nn));
    %CCM(nn,:) = CC;
    for ii = 1:Nloc % location
        g(ii,nn) = 0;
        count = 1;
        for nx = 0:ord
            for ny = 0:ord
                for nz = 0:ord
                    if nx+ny+nz<=ord
                        g(ii,nn) = g(ii,nn) + CC(count)*coor(ii,1)^nx*coor(ii,2)^ny*coor(ii,3)^nz;
                        count = count+1;
                    end
                end
            end
        end
    end
end

hlengthX = (lengthX-1)/2;
hlengthY = (lengthY-1)/2;
hlengthZ = (lengthZ-1)/2;
%% start filling in vOUT
vOUT= zeros(Nx,Ny,Nz);
%% top-bottom planes
xNBlk = idivide(int32(Nx),int32(lengthX));
xRemain = rem(Nx,lengthX);

zz = 1;
%initial y
yy = 1;
% integer part
for Xcount = 1:xNBlk
    xx = (Xcount-1)*lengthX+1;
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:hlengthZ %get the top
            for yyy = 1:lengthY
                for xxx = 1:lengthX
                    map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                    weightedM= map.*mIN_blk;
                    vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
                end
            end
        end
end
% remaider part
xx = Nx-lengthX+1;
mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));

for zzz = 1:hlengthZ %get the top
    for yyy = 1:lengthY
        for xxx = lengthX-xRemain+1:lengthX
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end

% sweeping y
% integer part
for yy = 2:(Ny-lengthY+1)
    for Xcount = 1:xNBlk
        xx = (Xcount-1)*lengthX+1;
        mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:hlengthZ %get the top
            
            for xxx = 1:lengthX
                yyy =lengthY;
                map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                weightedM= map.*mIN_blk;
                vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
            end
        end
    end
end
% remainder part
xx = Nx-lengthX+1;
for yy = 2:(Ny-lengthY+1)
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
    
    for zzz = 1:hlengthZ %get the top
        for xxx = lengthX-xRemain+1:lengthX
            yyy = lengthY;
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end

zz = Nz-lengthZ+1;
%initial y
yy = 1;
% integer part
for Xcount = 1:xNBlk
    xx = (Xcount-1)*lengthX+1;
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = hlengthZ+2:lengthZ %get the bottom
            for yyy = 1:lengthY
                for xxx = 1:lengthX
                    map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                    weightedM= map.*mIN_blk;
                    vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
                end
            end
        end
end
% remaider part
xx = Nx-lengthX+1;
mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));

for zzz = hlengthZ+2:lengthZ %get the bottom
    for yyy = 1:lengthY
        for xxx = lengthX-xRemain+1:lengthX
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end

% sweeping y
% integer part
for yy = 2:(Ny-lengthY+1)
    for Xcount = 1:xNBlk
        xx = (Xcount-1)*lengthX+1;
        mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = hlengthZ+2:lengthZ %get the bottom
            
            for xxx = 1:lengthX
                yyy =lengthY;
                map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                weightedM= map.*mIN_blk;
                vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
            end
        end
    end
end
% remainder part
xx = Nx-lengthX+1;
for yy = 2:(Ny-lengthY+1)
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
    
    for zzz = hlengthZ+2:lengthZ %get the bottom
        for xxx = lengthX-xRemain+1:lengthX
            yyy = lengthY;
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end
%% front-back planes
zNBlk = idivide(int32(Nz),int32(lengthZ));
zRemain = rem(Nz,lengthZ);
%front plane
yy = 1;
%initial y
xx = 1;
% integer part
for Zcount = 1:zNBlk
    zz = (Zcount-1)*lengthZ+1;
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:lengthZ 
            for yyy = 1:hlengthY %get the front
                for xxx = 1:lengthX
                    map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                    weightedM= map.*mIN_blk;
                    vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
                end
            end
        end
end
% remaider part
zz = Nz-lengthZ+1;
mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));

for zzz = lengthZ-zRemain+1:lengthZ 
    for yyy = 1:hlengthY %get the front
        for xxx = 1:lengthX
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end

% sweeping y
% integer part
for xx = 2:(Nx-lengthX+1)
    for Zcount = 1:zNBlk
        zz = (Zcount-1)*lengthZ+1;
        mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:lengthZ 
            for yyy = 1:hlengthY %get the front
                xxx =lengthX;
                map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                weightedM= map.*mIN_blk;
                vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
            end
        end
    end
end
% remainder part
zz = Nz-lengthZ+1;
for xx = 2:(Nx-lengthX+1)
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
    
    for zzz = lengthZ-zRemain+1:lengthZ
        for yyy = 1:hlengthY %get the front
            xxx = lengthX;
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end
% back-plane
yy = Ny-lengthY+1;
%initial y
xx = 1;
% integer part
for Zcount = 1:zNBlk
    zz = (Zcount-1)*lengthZ+1;
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:lengthZ 
            for yyy = hlengthY+2:lengthY %get the back
                for xxx = 1:lengthX
                    map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                    weightedM= map.*mIN_blk;
                    vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
                end
            end
        end
end
% remaider part
zz = Nz-lengthZ+1;
mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));

for zzz = lengthZ-zRemain+1:lengthZ 
    for yyy = hlengthY+2:lengthY %get the back
        for xxx = 1:lengthX
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end

% sweeping y
% integer part
for xx = 2:(Nx-lengthX+1)
    for Zcount = 1:zNBlk
        zz = (Zcount-1)*lengthZ+1;
        mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:lengthZ 
            for yyy = hlengthY+2:lengthY %get the back
                xxx =lengthX;
                map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                weightedM= map.*mIN_blk;
                vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
            end
        end
    end
end
% remainder part
zz = Nz-lengthZ+1;
for xx = 2:(Nx-lengthX+1)
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
    
    for zzz = lengthZ-zRemain+1:lengthZ
        for yyy = hlengthY+2:lengthY %get the back
            xxx = lengthX;
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end
%% left-right plane
% left plane
xx = 1;
%initial y
yy = 1;
% integer part
for Zcount = 1:zNBlk
    zz = (Zcount-1)*lengthZ+1;
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:lengthZ 
            for yyy = 1:lengthY 
                for xxx = 1:hlengthX %get the left
                    map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                    weightedM= map.*mIN_blk;
                    vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
                end
            end
        end
end
% remaider part
zz = Nz-lengthZ+1;
mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));

for zzz = lengthZ-zRemain+1:lengthZ 
    for yyy = 1:lengthY
        for xxx = 1:hlengthX %get the left
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end

% sweeping y
% integer part
for yy = 2:(Ny-hlengthY-lengthY+1)
    for Zcount = 1:zNBlk
        zz = (Zcount-1)*lengthZ+1;
        mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:lengthZ
            yyy = lengthY;
            for xxx =1:hlengthX;%get the left
                map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                weightedM= map.*mIN_blk;
                vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
            end
        end
    end
end
% remainder part
zz = Nz-lengthZ+1;
for yy = 2:(Ny-hlengthY-lengthY+1)
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
    
    for zzz = lengthZ-zRemain+1:lengthZ
        yyy = lengthY;
        for xxx =1:hlengthX;%get the left
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end
% right-plane
xx = Nx-lengthX+1;
%initial y
yy = 1;
% integer part
for Zcount = 1:zNBlk
    zz = (Zcount-1)*lengthZ+1;
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:lengthZ 
            for yyy = 1:lengthY 
                for xxx = hlengthX+2:lengthX %get the right
                    map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                    weightedM= map.*mIN_blk;
                    vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
                end
            end
        end
end
% remaider part
zz = Nz-lengthZ+1;
mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));

for zzz = lengthZ-zRemain+1:lengthZ 
    for yyy = 1:lengthY
        for xxx = hlengthX+2:lengthX %get the right
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end

% sweeping y
% integer part
for yy = 2:(Ny-hlengthY-lengthY+1)
    for Zcount = 1:zNBlk
        zz = (Zcount-1)*lengthZ+1;
        mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
        
        for zzz = 1:lengthZ
            yyy = lengthY;
            for xxx =hlengthX+2:lengthX %get the right
                map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
                weightedM= map.*mIN_blk;
                vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
            end
        end
    end
end
% remainder part
zz = Nz-lengthZ+1;
for yy = 2:(Ny-hlengthY-lengthY+1)
    mIN_blk = m3DIN(xx:(xx+lengthX-1),yy:(yy+lengthY-1),zz:(zz+lengthZ-1));
    
    for zzz = lengthZ-zRemain+1:lengthZ
        yyy = lengthY;
        for xxx =hlengthX+2:lengthX %get the right
            map = reshape(g(index3D(xxx,yyy,zzz),:),lengthX,lengthY,lengthZ);
            weightedM= map.*mIN_blk;
            vOUT(xx+xxx-1,yy+yyy-1,zz+zzz-1) = sum(sum(sum(weightedM)));
        end
    end
end
%% center block
center = index3D(hlengthX+1,hlengthY+1,hlengthZ+1);
mapCenter = reshape(g(index3D(hlengthX+1,hlengthY+1,hlengthZ+1),:),lengthX,lengthY,lengthZ);
for xx = 1:Nx-lengthX+1
    for yy = 1:Ny-lengthY+1
        for zz = 1:Nz-lengthZ+1
            mIN_blk = m3DIN(xx:xx+lengthX-1,yy:yy+lengthY-1,zz:zz+lengthZ-1);
            
            weightedM= mapCenter.*mIN_blk;
            vOUT(xx+hlengthX,yy+hlengthY,zz+hlengthZ) = sum(sum(sum(weightedM)));
        end
    end
end
%%
vDelta2 = abs(m3DIN-vOUT);