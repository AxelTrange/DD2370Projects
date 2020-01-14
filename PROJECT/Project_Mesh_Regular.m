%Fixing meshing of regular waveguide cross section:

function [el2no,no2xy,noInt,noExt] = Project_Mesh_Regular(NrPtsX,NrPtsY)

%Dimensions:
a = 0.02; %Width 2 cm
b = 0.01; %Height 1 cm

RegX = linspace(0,a,NrPtsX);
RegY = linspace(0,b,NrPtsY);

XList_Reg = [];
YList_Reg = [];

idx = 0;
for x = RegX
    for y = RegY
        idx = idx+1;
        XList_Reg(idx) = x;
        YList_Reg(idx) = y;
    end
end

%TRI = delaunay(XList_Reg,YList_Reg);

%object with TR.Points and TR.ConnectivityList
TR = delaunayTriangulation(XList_Reg', YList_Reg');
Edges = freeBoundary(TR);
Edge_1 = Edges(:,1);

%figure; triplot(TRI,XList_Reg,YList_Reg); hold on;
figure; triplot(TR.ConnectivityList,XList_Reg,YList_Reg); hold on;
xlabel('x [m]'); ylabel('y [m]');
axis equal;

%el2no = TRI';
el2no = TR.ConnectivityList';
no2xy = [XList_Reg; YList_Reg];
noInt = []; %No internal points or edges
noExt = Edge_1; %... Need to find these.

end