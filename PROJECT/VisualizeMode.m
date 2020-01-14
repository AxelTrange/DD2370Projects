function VisualizeMode(no2xy, el2no, vtr_Fz,label)

% Function:
%   VisualizeMode(no2xy, el2no, vtr_Fz)
% Arguments:
%   no2xy = the coordinates of the nodes
%   el2no = the nodes of the triangles
%   vtr_Fz = z-component of the eigenmode
% Returns:
%   -
% Comments:
%   Visualizes the z-component of the eigenmode together with its
%   curl. 

noNum = size(no2xy,2);  % size of column
elNum = size(el2no,2);

for elIdx = 1:elNum
    no = el2no(:,elIdx);
    xy = no2xy(:,no);
end

%plot grid
xy1 = no2xy(:,el2no(1,:));
xy2 = no2xy(:,el2no(2,:));
xy3 = no2xy(:,el2no(3,:));
xy = 100*[xy1; xy2; xy3; xy1; NaN*xy1];

x = 100*no2xy(1,:);
y = 100*no2xy(2,:);

% Find curl of eigenmode
vtr_Fx = zeros(length(vtr_Fz),1);
vtr_Fy = zeros(length(vtr_Fz),1);
% [curlx,curly,curlz,cav] = curl(vtr_Fx,vtr_Fy,vtr_Fz)
tri = el2no';
figure;trisurf(tri,x,y,abs(vtr_Fz)) % Error is from indexing!!!!
axis equal;
xlabel('x [m]'); ylabel('y [m]');title(label)

end
