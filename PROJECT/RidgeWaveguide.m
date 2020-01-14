function [el2no,no2xy,noExt] = RidgeWaveguide(Hmax,Lx,Ly,LxRidge,LyRidge)

Lx = 2e-2;  % WG width
Ly = 1e-2;  % WG height

 LxRidge = 1e-2; % Ridge width
 LyRH = 0.2e-2;  % Ridge heigth
LxWGf = 0.5*(Lx-LxRidge); % width of WG next to ridge


% Draw Polyshape Geometry
shape1 = polyshape([0 0 LxWGf LxWGf],[Ly 0 0 Ly]);
shape2 = polyshape([LxWGf LxWGf LxWGf+LxRidge LxWGf+LxRidge],[Ly-LyRidge LyRidge LyRidge Ly-LyRidge]);
shape3 = polyshape([Lx-LxWGf Lx-LxWGf Lx Lx],[Ly 0 0 Ly]);

shape12 = union(shape1,shape2);
RWG = union(shape12,shape3);

% Triangulation of Polyshape
tr = triangulation(RWG);
% triplot(tr)

% Convert Polyshape to pde object
model = createpde;

tnodes = tr.Points';
telements = tr.ConnectivityList';

geometryFromMesh(model,tnodes,telements);
% pdegplot(model);

% figure
% pdemesh(model);

% Generate mesh from pde object
MESH1 = generateMesh(model,'GeometricOrder','linear','Hmax',Hmax);
no2xy = MESH1.Nodes;
el2no = MESH1.Elements;
noExt = boundary(no2xy',1);

%TEST:
xx = no2xy(1,:);
yy = no2xy(2,:);
figure;plot(xx(noExt),yy(noExt));

% MESH2 = generateMesh(model,'GeometricOrder','linear','Hmax',0.5e-3)

figure
pdeplot(MESH1); title('MESH1');
xlabel('x [m]'); ylabel('y [m]')
% figure
% pdeplot(MESH2); title('MESH2');

end