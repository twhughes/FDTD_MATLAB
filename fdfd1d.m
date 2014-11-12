anim_fields = true;
%pause(10);
min_lambda = 2;
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
Ey = zeros(1,Nz);
Hx = zeros(1,Nz);
epsy = ones(1,Nz);
mux = ones(1,Nz);

%put in slab
thickness = Nz/50;
nz1 = floor(Nz/4);
nz2 = nz1+thickness;
figure(1);
areas = [];
for n = (floor(Nz/3):floor(thickness):floor(Nz*2/3))
    mux(n:n+floor(thickness/2)) = 1.0*2.0;
    epsy(n:n+floor(thickness/2)) = 1.0*6.0;
    areas = [areas, area([n,n+floor(thickness/2)], [1,1])];
end
nmax = sqrt(max(mux)*max(epsy));
nmin = sqrt(min(mux)*min(epsy));

%time step 
dt = nbc*dz/(2.0*c0);
fc = c0/nmax/min_lambda;          %highest frequency needed

%dsfsdf
Nbounces = 5;
T = dt*500000;                     % 6/fc + Nbounces*nmax*Nz*dz/c0    %time for n bounces
steps = 1 + floor(T/dt);         %approximate number of steps needed
trange = (0:dt:T);

%calculate source
loc_source = 2;  %1 unit cells from boundary |()(s)()()()()|
deltat = 3*dt/2;
A = -1.0*sqrt(epsy(loc_source)/mux(loc_source));   %???? normalize H and E to be same order where do eps and mu come from?

GEy = pulse(trange, fc, 0);  
GHx = pulse(trange, fc, deltat);
for i = (1:length(GHx))
    GHx(i) = A*GHx(i);
end

%Initialize Fourier Tranforms and kernels:

frange = (0:fc/1000.0:fc);

K = zeros(1,length(frange));
for i = (1:length(frange))
    K(i) = exp(-1j*2*pi*frange(i)*dt);
end

FT = zeros(1, length(frange));
FR = zeros(1, length(frange));
FS = zeros(1, length(frange));

mEy = dt*c0./epsy;
mHx = dt*c0./mux;

%Initialize Fields
%    Ey = zeros((Nz), float)    
%    Hx = zeros((Nz), float)    
Ey = zeros(1, Nz);
Hx = zeros(1, Nz);
h3 = 0;
h2 = 0;
h1 = 0;
e3 = 0;
e2 = 0;
e1 = 0;

%TIME LOOP
steps = floor(T/dt);
Esum = [];
Hsum = [];
EyX = [];
HxX = [];
FTX = [];
FRX = [];
FSX = [];
figure(1);
M = struct([getframe(gcf)]);
steps_ani = 0;

fprintf( 'Steps:     ', steps);
fprintf( 'Tmax:      ', T);
fprintf( 'max(mEy):  ', max(mEy));
fprintf( 'max(mHx):  ', max(mHx));
fprintf( 'dz:        ', dz);
fprintf( 'dt:        ', dt);
fprintf( 'Nz:        ', Nz);
fprintf( 'domain:    ', Nz*dz, ' meters');
fprintf( 'thickness: ', thickness*dz, ' meters');
fprintf( 'starting main time loop....');
fprintf( 'update F:  ', max(mEy)/dz);
fprintf( 'update H:  ', max(mHx)/dz);

disp( nbc*dz/2.0/c0 );
percentage = -1;
w = .04;
for t = (1:steps)
    %GEy(t) = sin(w*t);
    %GHx(t) = -sin(w*t);
    %disp(gEY)
%perfect boundary conditions                                                              

    %update H from E
                Hx(1:Nz-1) = Hx(1:Nz-1) + mHx(1:Nz-1).*(Ey(2:Nz) - Ey(1:Nz-1))/dz;
   % for i = (1:Nz-1)
   %      Hx(i) = Hx(i) + mHx(i)*(Ey(i+1) - Ey(i))/dz;
   % end
    Hx(Nz) = Hx(Nz) + mHx(Nz)*(e3 - Ey(Nz))/dz;

    %update H source
    Hx(loc_source) = Hx(loc_source) - mHx(loc_source)/dz*GHx(t);

    %record H at boundary
    h3 = h2;
    h2 = h1;
    h1 = Hx(1);
    %h1 = 0;%

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
    %e1 = 0;

    %Update Fourier Tranforms
    FS = FS + GEy(t).*(K.^t)*dt;
    if (t > 1000)
        FT = FT + Ey(Nz).*(K.^t)*dt;
        FR = FR + Ey(1).*(K.^t)*dt;
    end
    if anim_fields
        figure(1);
        for i = (1:length(areas))
            a = areas(i);
        end
        clf;
        hold on;
        basevalue = -1;   
        plot((1:Nz), Ey, 'g');
        plot((1:Nz), Hx, 'b');
        ylim([-1,1]);
   %     annotation('textbox', [.75,0.75,.05,0.05],...
   %        'String', int2str(t));
        hold off;
    else
        figure(2);
      %  M(i)=getframe(gcf);
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
    xlim([0,Nz]);
    ylim([-1,1]);
    pause(0.01);
    hold off;
    %}
    %visualize solution
    Esum = [Esum, (sum(Ey))];
    Hsum = [Hsum, (sum(Hx))];

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
plot(frange, ((abs(FR).^2)+(abs(FT).^2))./(abs(FS).^2), 'b');
hold off;