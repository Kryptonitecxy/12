%% 配电网风机和储能优化配置（决策变量：风机、储能位置和容量(台数)，并网逆变器容量）
clc;
clear;
close all;

powerWind;
clearvars -except P_wind_pre pi %仅保留n个场景的风机出力数据及对应概率
global P_wind pi gxbest1 Power_st all_Cost;
parameter; %基本参数定义

%变量初始化x=(adwind1,adwind2,...,adst1,adst2,...,numwd1,numwd2,...,numst1,numst2,...,Inverters)
%对应为(风机1选址、风机2选址...储能1选址、储能2选址...风机1台数、风机2台数...储能1台数、储能2台数...并网逆变器容量)
for i=1:N    
    %位置初始化(不允许重复)
    A=ad_wd(randperm(numel(ad_wd),Wd));    %numel(ad_wd)风机待选节点总数
    B=ad_st(randperm(numel(ad_st),St));
    %台数初始化(允许重复)
    C=num_wd(randi(numel(num_wd),1,Wd));   %randi可重复，randperm不可重复
    d=num_st(randi(numel(num_st),1,St));   %范围在 1 到 numel(num_wd) 之间。
    %并网逆变器容量初始化(kW)
    E=1200.*rand();
    x(i,:)=horzcat(A,B,C,d,E);
end

%开始迭代,上层pso返回各类最优数据
[uu,pg,c,g,po,ac]=up_pso(x);

%数据处理
Wind_power=[]; %运行阶段风机出力
SOC=[]; %运行阶段储能SOC变化
P_SOC=[]; %运行阶段储能出力变化
Inv_power=[]; %运行阶段并网逆变器出力
for i=1:K*Wd
    Wind_power(i,:)=g(1+T*(i-1):T*i);
end
for i=1:K*St
    SOC(i,:)=po(1+(T+1)*(i-1):(T+1)*i);
    P_SOC(i,:)=-po((T+1)*K*St+1+T*(i-1):(T+1)*K*St+T*i);
end
for i=1:K
    Inv_power(i,:)=po((T+1)*K*St+T*K*St+1+T*(i-1):(T+1)*K*St+T*K*St+T*i);
end
disp('总运维成本:(元)');
Cost_yw=ac(1)
disp('总购售电成本:(元)');
Cost_gou=ac(2)
disp('总弃风成本:(元)');
Cost_qi=ac(3)
disp('总失电成本:(元)');
Cost_loss=ac(4)
disp('总网损:');
Ploss=ac(5)
disp('电压偏移量:');
Vp=ac(6)
disp('总配置成本(万元):');
Cost_pei=c
disp('电网可靠性');
V_re=Ploss+Vp
disp('总运行成本(万元)');
Cost_yunxing=(Cost_yw+Cost_gou+Cost_qi+Cost_loss)*20/10000
mpc=loadcase('case33bw');
for i=1:T
    PL(i)=pl(i).*sum(mpc.bus(:,3)); %负荷  
end

%绘图
figure(4);
plot(uu(:,1),'Color',[0.54118 0.16863 0.88627],'MarkerFaceColor','k','LineWidth',2);
title('总目标函数迭代收敛图');
xlabel('迭代次数');
ylabel('目标函数大小');

figure(5);
plot3(uu(:,2),uu(:,3),uu(:,4));
xlabel('总运行成本')
ylabel('电网可靠性')
zlabel('总配置成本')


    
    





