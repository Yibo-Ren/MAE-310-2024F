clc;
clear;

%导入mesh文件生成网格

filename = 'quarter-plate-with-hole-quad-12112216任奕博.msh';
fid = fopen(filename, 'r');
if fid == -1
    error('无法打开文件');
end

% 读取文件内容
line = fgetl(fid);
while ischar(line)
    if contains(line, '$Nodes')             % 读取节点数据
        n_np = str2double(fgetl(fid));
        nodes = zeros(n_np, 4);             % nodal index, x, y, z
        for i = 1:n_np
            nodes(i, :) = sscanf(fgetl(fid), '%d %f %f %f')';
        end
    elseif contains(line, '$Elements')      % 读取单元数据
        n_el = str2double(fgetl(fid));
        elements = cell(n_el, 1);
        for i = 1:n_el
            elements{i} = sscanf(fgetl(fid), '%d')';
        end
        disp('单元数据已读取');
    end
    line = fgetl(fid);
end
fclose(fid);

% 显示结果
disp('节点:');
disp(nodes);

disp('单元:');
disp(elements);
