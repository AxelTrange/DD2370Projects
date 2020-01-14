function Visualize_Curl(no2xy,el2no,Mode,Fz,Field)


DerivX = zeros(1,size(no2xy(1,:),2));
DerivY = zeros(1,size(no2xy(1,:),2));
for no = 1:size(no2xy,2)
   
    %Here no is the element number!
    no_xy = no2xy(:,no);
    
    Common_no = ismember(el2no,no); %Check if nodes are shared
    NrCommon = sum(Common_no); %Nr of common nodes
    Ind_Common = find(NrCommon == 1); %indices
    
    Adj_tri_no = el2no(:,Ind_Common);
    Adj_tri_xy = no2xy(:,Adj_tri_no);

    Uni_el = unique(Adj_tri_no);
    Uni_xy = unique(Adj_tri_xy','rows');
    
    Diff = no_xy - Uni_xy';
    
    DeltaX = Diff(1,:);
    DeltaY = Diff(2,:);
    d = sqrt(sum(Diff.^2));
    %f1_H = Fz(no,Mode); %for first mode
    f1_H = Fz(no); %for first mode
    %f2_H = Fz(Uni_el,Mode); %for first mode
    f2_H = Fz(Uni_el); %for first mode
    
    DerivAllX = (f1_H -  f2_H')./d .* 1./(abs(DeltaX)+abs(DeltaY)).*DeltaX;
    DerivX(no) = abs(sum(DerivAllX(~isnan(DerivAllX))));
    DerivAllY = (f1_H -  f2_H')./d .* 1./(abs(DeltaX)+abs(DeltaY)).*DeltaY;
    DerivY(no) = abs(sum(DerivAllY(~isnan(DerivAllY))));
    
end

Curl_mag_H = sqrt(DerivX.^2 + DerivY.^2);
%Et_field = Curl_mag_H/(2*pi*f*eps0); We just normalize it instead.
figure;trisurf(el2no',no2xy(1,:),no2xy(2,:),Curl_mag_H/max(Curl_mag_H));
xlabel('x [m]');ylabel('y [m]');title(Field); axis equal;

end