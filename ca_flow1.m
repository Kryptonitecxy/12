clc;
clear;
%% 加载case33数据
mpc = loadcase('case33bw');

% 运行潮流计算
results = runpf(mpc);
 
% 显示结果
printpf(results);

% 提取节点数据
bus = results.bus;
branch = results.branch;
 
% 提取节点位置（这里假设有预定义的位置数据）
% 可以根据实际情况设定节点的 (x, y) 坐标
% 例如：bus_positions = [x1, y1; x2, y2; ...];
bus_positions = [0,0;1,0;2,0;3,0;4,0;5,0;6,0;7,0;8,0;9,0;10,0;11,0;12,0;13,0;14,0;15,0;16,0;17,0;
                 2,-1;3,-1;4,-1;5,-1;
                 3,2;4,2;5,2;
                 6,1;7,1;8,1;9,1;10,1;11,1;12,1;13,1];
 
% 绘制节点
figure;
hold on;
scatter(bus_positions(:, 1), bus_positions(:, 2), 'filled');
text(bus_positions(:, 1), bus_positions(:, 2), num2str(bus(:, 1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
 
% 绘制线路并显示潮流分布
for k = 1:size(branch, 1)
    fbus = branch(k, 1);  % 起始节点
    tbus = branch(k, 2);  % 终止节点
    x = [bus_positions(fbus, 1), bus_positions(tbus, 1)];
    y = [bus_positions(fbus, 2), bus_positions(tbus, 2)];
    
    % 获取潮流数据
    Pf = branch(k, 14);  % 起始节点的有功功率
    Qf = branch(k, 15);  % 起始节点的无功功率
    Pt = branch(k, 16);  % 终止节点的有功功率
    Qt = branch(k, 17);  % 终止节点的无功功率
    
    % 绘制线路
    plot(x, y, 'k');
    
    % 绘制有功潮流箭头
    quiver(x(1), y(1), x(2)-x(1), y(2)-y(1), 'MaxHeadSize', 0.1, 'Color', 'r');
    
    % 在箭头上标注潮流数据
%     mid_x = (x(1) + x(2)) / 2;
%     mid_y = (y(1) + y(2)) / 2;
%     text(mid_x, mid_y, sprintf('%.2f MW', Pf), 'Color', 'r', 'FontSize', 8);
end
 
% 设置图形属性
xlabel('X 坐标');
ylabel('Y 坐标');
title('电力系统潮流分布图');
grid on;
hold off;