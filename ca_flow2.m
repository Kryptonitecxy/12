clc;
clear;
%% 加载case33数据
mpc0 = loadcase('case33bw');

% 运行潮流计算
results1 = runpf(mpc0);
 
% 显示结果
% printpf(results);
% plot(results1.bus(:,8));

mpc = loadcase('case33bw');
A=0.7293.*sum(mpc.bus(:,3))
mpc.bus(3,3) = mpc.bus(3,3)-3.715-0.1529;

% mpc.bus(17,3) = mpc.bus(17,3)-4.8;
% mpc.bus(28,3) = mpc.bus(28,3)-2.2;
% for j = 1:33
%     if mpc.bus(j,3) ~= 0
%         tan(j)= mpc0.bus(j,4)./mpc0.bus(j,3);
%         mpc.bus(j,3) = normrnd(1.05*mpc0.bus(j,3),1.05*mpc0.bus(j,3)*0.05);
%         mpc.bus(j,4) = mpc.bus(j,3).*tan(j);
%     end
% end
results2 = runpf(mpc);
B=sum(results2.branch(:,14)+results2.branch(:,16))
C=sum(results1.branch(:,14)+results1.branch(:,16))
plot(results1.bus(:,8));
hold on;
plot(results2.bus(:,8));
legend('原始','风电并网');

