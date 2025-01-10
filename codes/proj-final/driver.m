clc;
clear;

run('quarter_plate.m') %导入文件获得节点坐标

E = 1E9;                %Young's modular
nu = 0.3;               %Poison's ratio

%exact stress in polar coor
sigma_rr = @(Tx,R,r,theta) Tx./2.*(1 - R.^2./r.^2) + Tx./2.*(1 - 4.*R.^2./r.^2 + 3.*R.^4./r.^4).*cos(2.*theta);
sigma_tt = @(Tx,R,r,theta) Tx./2.*(1 + R.^2./r.^2) - Tx./2.*(1 + 3.*R.^4./r.^4).*cos(2.*theta);
sigma_rt = @(Tx,R,r,theta) -Tx./2.*(1 + 2.*R.^2./r.^2 - 3.*R.^4./r.^4).*sin(2.*theta);

nid = msh.POS(:,1) == 1; % find the boundary element
x_coor = zeros(length(nid),1);
y_coor = x_coor;
for i =1 : length(nid)
    if nid(i) == 1
        x_coor(i) = msh.POS(i,1);
        y_coor(i) = msh.POS(i,2);    
    end
end


Tx = 1e4; % Boundary condition 
R = 0.3;
r = sqrt((x_coor-(-1)).^2 + (y_coor-(-1)).^2);
theta = atan((y_coor-(-1))./((x_coor-(-1))));
%get sigma_rr/sigma_tt/sigma_rt
sigma_r = sigma_rr(Tx,R,r,theta);
sigma_t = sigma_tt(Tx,R,r,theta);
sigma_r_t = sigma_rt(Tx,R,r,theta);
%convert into Cartesian coor,which is f
sigma_xx = sigma_r .* cos(theta).^2 + sigma_t .* sin(theta).^2 + 2.*sigma_r_t.*sin(theta).*cos(theta);
sigma_yy = sigma_r.*sin(theta).^2 + sigma_t.*cos(theta).^2 - 2.*sigma_r_t.*sin(theta).*cos(theta);
sigma_xy = (sigma_t-sigma_r).*sin(theta).*cos(theta) + sigma_r_t.*(cos(theta).^2-sin(theta).^2);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

n_np = msh.nbNod;               %number of nodal point
n_el = length(msh.QUADS(:,1));  %number of element
n_en = 4;                       %number of point in an element

x_coor = zeros(n_np,1);
y_coor = x_coor;

for i = 1 : n_np                %组件x，y坐标
    x_coor(i,1) = msh.POS(i,1);
    y_coor(i,1) = msh.POS(i,2);
end
%IEN array
IEN = zeros(n_el,n_en);
for i = 1 : n_el
    for ii = 1 : n_en
        IEN(i,ii) = msh.QUADS(i,ii);
    end
end
%ID array
ID = zeros(2,n_np);
counter = 0;
for i = 1 : 2
    for ii = 1 : n_np
       if msh.POS(ii,i) == -1
           ID(i,ii) = 0;
       else
           counter = counter + 1;
           ID(i,ii) = counter;
       end
    end
end
n_eq = counter;

LM = zeros(2,n_el,n_en);
for i = 1:2
    for ii = 1:n_el
        for iii = 1:n_en
            LM(i,ii,iii) = ID(i,IEN(ii,iii));
        end
    end
end

%loop over element to assembly the matrix and vector
D = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];

F = zeros(n_eq, 1);
K = zeros(n_eq, n_eq);

for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en*2, n_en*2); % element stiffness matrix
  f_ele = zeros(n_en*2, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    

    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      
      %f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
      f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * sigma_xx(IEN(ee,aa)) * Na;
      f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * sigma_yy(IEN(ee,aa)) * Na;

      B = zeros(3,2);
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
        B(1,2) = Nb_x;
        B(2,2) = Nb_y;
        B(3,1) = Nb_y;
        B(3,2) = Nb_x;
        
        k_ele(2*aa-1:2*aa,2*bb-1:2*bb) = B'*D*B;
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
 
  for aa = 1 : n_en
    PPx = LM(1,ee, aa);
    PPy = LM(2,ee, aa);
    if PPx > 0
      F(PPx) = F(PPx) + f_ele(2*aa-1);
      
      for bb = 1 : n_en
        QQx = LM(1,ee, bb);
        QQy = LM(2,ee, bb);

        if QQx > 0
          K(PPx, QQx) = K(PPx, QQx) + k_ele(2*aa-1, 2*bb-1);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
        if QQy > 0
           K(PPx, QQy) = K(PPx, QQy) + k_ele(2*aa-1, 2*bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end

    if PPy > 0
       F(PPy) = F(PPy) + f_ele(2*aa);

       for bb = 1 : n_en
           QQx = LM(1,ee, bb);
           QQy = LM(2,ee, bb);
           if QQx > 0
              K(PPy, QQx) = K(PPy, QQx) + k_ele(2*aa, 2*bb-1);
           else
                    % modify F with the boundary data 
                    % here we do nothing because the boundary data g is zero or
                    % homogeneous
           end
           if QQy > 0
             K(PPy, QQy) = K(PPy, QQy) + k_ele(2*aa, 2*bb);
           else
                    % modify F with the boundary data
                    % here we do nothing because the boundary data g is zero or
                    % homogeneous
           end
        end
     end
  end
end


% solve the stiffness matrix
dn = K \ F;

% % insert dn back into the vector for all nodes
% disp = zeros(n_np, 1);
% 
% for ii = 1 : n_np
%   index = ID(ii);
%   if index > 0
%     disp(ii) = dn(index);
%   else
%     % modify disp with the g data. Here it does nothing because g is zero
%   end
% end