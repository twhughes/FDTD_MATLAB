anim_fields = true;
min_lambda = 5;
dmin = 3000.48;
dc = 3000.48;
Nres = 1;
Ndres = 1;
%fundamental constants
nmax = 1;
nsrc = 1;
c0 = 2.99179e8;    %speed of light (m/s)
e0 = 8.854e-12;    %epsilon_0
m0 = 1.2566e-6;    %mu_0

disp( 'initializing.....' );


%compute grid resultion
nbc = 1;
dz1 = min_lambda/nmax/Nres ;
dz2 = dmin/Ndres;
dz  = min(dz1, dz2);
%snap to grid critical dimenstions
N  = (1+floor(dc/dz)); %device cells
dz = dc/N/30;
Nz = N;
Nx = Nz;   %hack for now
Ny = Nz;
dx = dz;
dy = dz;

epsx = ones(Nx,Ny);
epsy = ones(Nx,Ny);
mux = ones(Nx,Ny);
muy = ones(Nx,Ny);

%put in slab
thickness = Nz/20;
nz1 = floor(Nz/4);
nz2 = nz1+thickness;
%mux(nz1:nz2) = 1.0*2.0;
%epsy(nz1:nz2) = 1.0*6.0;
%nmax = sqrt(max(mux)*max(epsy));
%nmin = sqrt(min(mux)*min(epsy));
nmax = 1;
nmin = 1;

%time step 
dt = nbc*dz/(2.0*c0);
fc = c0/nmax/min_lambda;          %highest frequency needed

%dsfsdf
Nbounces = 5;
T = dt*5000;                     % 6/fc + Nbounces*nmax*Nz*dz/c0    %time for n bounces
steps = 1 + floor(T/dt);         %approximate number of steps needed
trange = (0:dt:T);

%calculate source
loc_source = [20,20];
deltat = 3*dt/2;

A = -1.0*sqrt(epsy(loc_source)/mux(loc_source));   %incorporate
GEx = pulse(trange, fc, 0);  
GEy = pulse(trange, fc, 0);  
GHx = pulse(trange, fc, deltat);
GHy = pulse(trange, fc, deltat);
for i = (1:length(GHx))
    GHx(i) = A*GHx(i);  %make faster
end

mEx = dt*c0./epsx;
mEy = dt*c0./epsy;
mHx = dt*c0./mux;
mHy = dt*c0./muy;

%TIME LOOP
steps = floor(T/dt);

%initialize curls and fields
CEx = zeros(Nx, Ny);
CEy = zeros(Nx, Ny);
CEz = zeros(Nx, Ny);
CHx = zeros(Nx, Ny);
CHy = zeros(Nx, Ny);
CHz = zeros(Nx, Ny);
Ex = zeros(Nx, Ny);
Ey = zeros(Nx, Ny);
Ez = zeros(Nx, Ny);
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
Hz = zeros(Nx, Ny);

for t = (1:steps)

    %compute CEx
    for nx = 1 : Nx
        for ny = 1 : Ny-1
            for nz = 1 : Nz-1
                CEx(nx,ny) = (Ez(nx,ny+1) - Ez(nx,ny))/dy ...
                - (Ey(nx,ny+1) - Ey(nx,ny))/dz;
            end
            CEx(nx,ny) = (Ez(nx,ny+1) - Ez(nx,ny))/dy ...
            - (Ey(nx,ny,1) - Ey(nx,ny))/dz;
        end
        for nz = 1 : Nz-1
            CEx(nx,Ny) = (Ez(nx,1) - Ez(nx,Ny))/dy ...
            - (Ey(nx,Ny+1) - Ey(nx,Ny))/dz;
        end
        CEx(nx,Ny) = (Ez(nx,1) - Ez(nx,Ny))/dy ...
        - (Ey(nx,Ny,1) - Ey(nx,Ny))/dz;
    end

    %compute CEy
    for nx = 1 : Nx
        for ny = 1 : Ny-1
            CEy(nx,ny) = - (Ez(nx,ny,1) - Ez(nx,ny))/dx;
        end
        CEy(nx,Ny) = - (Ez(nx,Ny,1) - Ez(nx,Ny))/dx;
    end
    
    
    
    Hx(1:Nz-1) = Hx(1:Nz-1) + mHx(1:Nz-1).*(Ey(2:Nz) - Ey(1:Nz-1))/dz;
    Hx(Nz) = Hx(Nz) + mHx(Nz)*(e3 - Ey(Nz))/dz;

    %update H source
    Hx(loc_source) = Hx(loc_source) - mHx(loc_source)/dz*GHx(t);

    %record H at boundary
    h3 = h2;
    h2 = h1;
    h1 = Hx(1);

    %Update E from H
    Ey(1) = Ey(1) + mEy(1)*(Hx(1) - h3)/dz;
    %for i = (2:Nz)
    %    Ey(i) = Ey(i) + mEy(i)*(Hx(i) - Hx(i-1))/dz;
    %end
    Ey(2:Nz) = Ey(2:Nz) + mEy(2:Nz).*(Hx(2:Nz) - Hx(1:Nz-1))/dz;
    
    
    %Update E source
    Ey(loc_source) = Ey(loc_source) - mEy(loc_source)/dz*GEy(t);

    %record E at boundaries
    e3 = e2;
    e2 = e1;
    e1 = Ey(Nz);

    %Update Fourier Tranforms
    FS = FS + GEy(t).*(K.^t)*dt;
    if (t > 0)
        FT = FT + Ey(Nz).*(K.^t)*dt;
        FR = FR + Ey(1).*(K.^t)*dt;
    end
    if anim_fields
        figure(1);
        clf;
        hold on;
        basevalue = -1;
        ha = area([nz12], [1,1]);
        ha2 = area([nz12], [-1,-1]);    
        plot((1:Nz), Ey, 'g');
        plot((1:Nz), Hx, 'b');
        ylim([-1,1]);
   %     annotation('textbox', [.75,0.75,.05,0.05],...
   %        'String', int2str(t));
        hold off;
    else
        figure(2);
        clf;
        hold on;
        plot(frange,abs(FT), 'r');
        plot(frange,abs(FR), 'g');
        plot(frange,abs(FS), 'b');
        %xlim(frange);
        %ylim([-1,1]);
        hold off;
    end
    %{
    clf;
    hold on;
    plot((1:length(FT)), abs(FT), 'r');
    plot((1:length(FT)), abs(FR), 'k');
    xlim([0]);
    ylim([-1,1]);
    pause(0.01);
    hold off;
    %}
    %visualize solution
    Esum = [Esum, (sum(Ey))];


    %Fourier post processing
end
%%
%FT = FT.*dt;
%FR = FR.*dt;
%FS = FS.*dt;
clf;
hold on;
plot(frange, (abs(FT).^2)./(abs(FS).^2), 'g');
plot(frange, (abs(FR).^2)./(abs(FS).^2), 'r');
plot(frange, ((abs(FR).^2)+(abs(FT).^2))./(abs(FS).^2), 'r');
hold off;