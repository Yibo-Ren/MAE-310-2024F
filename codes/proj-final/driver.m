clc;
clear;

run('quarter_plate.m') %导入文件获得节点坐标

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
ID = zeros(2,n_np,1);
counter = 1;
for i = 1 : 2
    for ii = 1 : n_np
       if msh.POS(ii,i) == -1
           ID(i,ii) = 0;
       else
           ID(i,ii) = counter;
           counter = counter + 1;
       end
    end
end
