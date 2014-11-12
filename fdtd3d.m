  %compute CEx
    for nx = 1 : Nx
        for ny = 1 : Ny-1
            for nz = 1 : Nz-1
                CEx(nx,ny,nz) = (Ez(nx,ny+1,nz) - Ez(nx,ny,nz))/dy ...
                - (Ey(nx,ny,nz+1) - Ey(nx,ny,nz))/dz;
            end
            CEx(nx,ny,Nz) = (Ez(nx,ny+1,Nz) - Ez(nx,ny,Nz))/dy ...
            - (Ey(nx,ny,1) - Ey(nx,ny,Nz))/dz;
        end
        for nz = 1 : Nz-1
            CEx(nx,Ny,nz) = (Ez(nx,1,nz) - Ez(nx,Ny,nz))/dy ...
            - (Ey(nx,Ny,nz+1) - Ey(nx,Ny,nz))/dz;
        end
        CEx(nx,Ny,Nz) = (Ez(nx,1,Nz) - Ez(nx,Ny,Nz))/dy ...
        - (Ey(nx,Ny,1) - Ey(nx,Ny,Nz))/dz;
    end

    %compute CEy
    for nx = 1 : Nx
        for ny = 1 : Ny-1
            for nz = 1 : Nz-1
                CEy(nx,ny,nz) = (Ex(nx,ny+1,nz) - Ex(nx,ny,nz))/dz ...
                - (Ez(nx,ny,nz+1) - Ez(nx,ny,nz))/dx;
            end
            CEy(nx,ny,Nz) = (Ex(nx,ny+1,Nz) - Ex(nx,ny,Nz))/dz ...
            - (Ez(nx,ny,1) - Ez(nx,ny,Nz))/dx;
        end
        for nz = 1 : Nz-1
            CEy(nx,Ny,nz) = (Ex(nx,1,nz) - Ex(nx,Ny,nz))/dz ...
            - (Ez(nx,Ny,nz+1) - Ez(nx,Ny,nz))/dx;
        end
        CEy(nx,Ny,Nz) = (Ex(nx,1,Nz) - Ex(nx,Ny,Nz))/dz ...
        - (Ez(nx,Ny,1) - Ez(nx,Ny,Nz))/dx;
    end
    
    %compute CEz
    for nx = 1 : Nx
        for ny = 1 : Ny-1
            for nz = 1 : Nz-1
                CEz(nx,ny,nz) = (Ey(nx,ny+1,nz) - Ey(nx,ny,nz))/dx ...
                - (Ex(nx,ny,nz+1) - Ex(nx,ny,nz))/dy;
            end
            CEz(nx,ny,Nz) = (Ey(nx,ny+1,Nz) - Ey(nx,ny,Nz))/dx ...
            - (Ex(nx,ny,1) - Ex(nx,ny,Nz))/dy;
        end
        for nz = 1 : Nz-1
            CEz(nx,Ny,nz) = (Ey(nx,1,nz) - Ey(nx,Ny,nz))/dx ...
            - (Ex(nx,Ny,nz+1) - Ex(nx,Ny,nz))/dy;
        end
        CEz(nx,Ny,Nz) = (Ey(nx,1,Nz) - Ey(nx,Ny,Nz))/dx ...
        - (Ex(nx,Ny,1) - Ex(nx,Ny,Nz))/dy;
    end   
    