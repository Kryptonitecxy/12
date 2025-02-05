clear,clc
parameter;
times=200;%场景数
global P_wind pi
shuju=xlsread('data.xlsx'); %把一天划分为24小时
PV_a1=shuju(4,:);%光伏功率24个不同时刻下的参数1
PV_b1=shuju(5,:);%光伏功率24个不同时刻下的参数2
PV_a=PV_a1';%光伏功率24个不同时刻下的参数1
PV_b=PV_b1';
WT_c1=shuju(6,:);%风电参数 24个不同时刻下的参数1
WT_k1=shuju(7,:);%风电参数 24个不同时刻下的参数2
WT_c=WT_c1';
WT_k=WT_k1';

for t=1:24
Ppv_samp=zeros(1,times);  %time 是蒙特卡洛仿真的次数，也是场景数
a_pv=PV_a(t); 
b_pv=PV_b(t);  
 if a_pv>0
% 光伏发电相关参数：光伏额定容量P_CAP、光强基准值Bg、额定光强Ge(kW/m2)
P_CAP=275;%随意给出，自己调节
Bg=2;%随意给出，自己调节
Ge=0.7;%随意给出，自己调节
% 光伏有功出力样本
pv_samp(1,:)=betarnd(a_pv,b_pv,1,times);%生成形状为(1,times,)，参数为a_pv,b_pv的Beta分布
 Ppv_samp(1,:)=pv_samp(1,:)*P_CAP*Bg/Ge;% 光伏功率=Beta分布*光伏额定容量*光强基准值/额定光强
 else
    Ppv_samp=zeros(1,times) ;
 end %对应的是if a_pv>0  
 
 
 %%==风电%%%%%
k_wt=WT_k(t);%风电参数 24个不同时刻下的参数1
c_wt=WT_c(t);%风电参数 24个不同时刻下的参数2
wt_samp =wblrnd(c_wt,k_wt,1,times);  % 风速 产生服从weibull分布的样本 ，形状为（(1,times,)
vci=3 ;%随意给出，自己调节; 切入风速
vN=13;%随意给出，自己调节  ： 额定风速
vco=25 ;%随意给出，自己调节 切出风速
    for i=1:times  %得到风电出力样本
        if wt_samp(i)<vci %如果风速小于切入风速
            Pwt_samp(i)=0; %风机功率为0
        end %对应if wt_samp(i)<vci
         if wt_samp(i)>vci&&wt_samp(i)<vN %如果风速大于切入风速，同时小于额定风速
            Pwt_samp(i)=(wt_samp(i)-vci)/(vN-vci)*s_wd;%式子2-4中的2   ，这是当成了线性的关系？ 
            if   Pwt_samp(i)>s_wd %如果风电功率大于额定功率
                 Pwt_samp(i)=s_wd; %则风电功率等于额定功率
            end %对应 if   Pwt_samp(i)>PN_wt
         end
         if wt_samp(i)>vN&&wt_samp(i)< vco %如果风速大于额定风速 同时小于切出风速
            Pwt_samp(i)=s_wd;              %风电功率等于额定功率
         end                               %对应 if wt_samp(i)>vN&&wt_samp(i)< vco
         if wt_samp(i)> vco                %如果风速大于切出风速
            Pwt_samp(i)=0;                 %风电功率等于0
         end                               %对应 if wt_samp(i)> vco
    end%对应for i=1:times%得到风电出力样本
    for i=1:times
     Pwt_samp(i)= Pwt_samp(i) ; %存储每个场景下 ，t时刻 风电发功率
    end
    
 P_pv1(t,:)=  Ppv_samp;%存储t时刻 所有场景的光伏功率
 P_wd1(t,:)=  Pwt_samp;%存储t时刻 所有场景的风电功率

end%对应的是for t=1:24

   


%电热负荷服从正态分布
D_fuhe1=shuju(2,:);
R_fuhe1=shuju(3,:);
Wz= 2000;   %误差量
D_fangcha=zeros(1,24);
R_fangcha=zeros(1,24);
[m,T]=size(D_fuhe1);
for i=1:T
    D_fangcha(1,i)=D_fuhe1(1,i)/5+Wz/50; %计算方差 
    R_fangcha(1,i)=R_fuhe1(1,i)/5+Wz/50;
end
D=zeros(200,24);
R=zeros(200,24);
for i=1:T
    D (:,i)= normrnd(0,D_fangcha(1,i),1,200);%生成200组服从均值为0，方差为  的正态分布随机数
    R (:,i)= normrnd(0,R_fangcha(1,i),1,200);
end
for i=1:200
    D_fuhe(i,:)=D_fuhe1(1,:)+D(i,:);%预测+误差
    R_fuhe(i,:)=R_fuhe1(1,:)+R(i,:);
end

P_wind_pre=[P_pv1'  P_wd1' D_fuhe R_fuhe]; %存储光伏 负荷数据


%场景削减
pi=1/200*ones(1,200);%每个场景的概率1/200
number=4;%缩减至场景数
N=times;%蒙特卡洛数场景数
while (N~=number)% 如果数不等于场景数
    %计算每个样本的欧氏距离
    for i=1:N %遍历每一个场景
        l(i,i)=inf; %第i个场景 和第i个场景的距离为无穷大
        for j=1:N %遍历每一个场景
            if(i~=j)%如果两个场景编号不一样
                l(i,j)=0;%初始化两个场景距离为0
                for k=1:24
                    l(i,j)=l(i,j)+(P_wind_pre(i,k)-P_wind_pre(j,k))^2;%计算场景i与场景j的距离
                end %对应for k=1:24
                l(i,j)=sqrt(l(i,j));%距离开平方
            end%对应if(i~=j)%如果两个场景编号不一样
        end%对应 for j=1:N %遍历每一个场景
    end%对应 for i=1:N %遍历每一个场景
     
    %计算样本的相应值
    for i=1:N%遍历每一个场景
        [m(i),n(i)]=min(l(i,:));%m(i)为最小值的值，n(i)为最小值存在的索引号
        d(i)=pi(i)*m(i);
    end %对应  %计算样本的相应值 for i=1:N
     
    %找出样本最小的
    [s,I]=min(d);
    %更新
    pi(n(I))=pi(n(I))+pi(I);%场景概率
    P_wind_pre(I,:)=[];
    pi(I)=[];
    l(I,:)=[];
    l(:,I)=[];
    m(I)=[];
    n(I)=[];
    d(I)=[];
    N=N-1;
end%对应while (N~=number)% 如果数不等于场景数

global pi; %概率 
P_pv2=zeros(number,24);
P_wd2=zeros(number,24);
D_fuhe2=zeros(number,24);
R_fuhe2=zeros(number,24);
for j=1:number
    for i=1:24
         P_pv2(j,i)=P_pv2(j,i)+pi(j)*P_wind_pre(j,i);
         P_wd2(j,i)=P_wd2(j,i)+pi(j)*P_wind_pre(j,i+24);
         D_fuhe2(j,i)=D_fuhe2(j,i)+pi(j)*P_wind_pre(j,i+48);
         R_fuhe2(j,i)=R_fuhe2(j,i)+pi(j)*P_wind_pre(j,i+72);
    end
end
P_PV=zeros(1,24);
P_WD=zeros(1,24);
OP_load_e=zeros(1,24);
OP_load_h=zeros(1,24);

for i=1:24
    for j=1:number
         P_PV(1,i)=P_PV(1,i)+P_pv2(j,i);
         P_WD(1,i)=P_PV(1,i)+P_wd2(j,i);
         OP_load_e(1,i)=OP_load_e(1,i)+D_fuhe2(j,i);
         OP_load_h(1,i)=OP_load_h(1,i)+R_fuhe2(j,i);
    end
end

%% 作图--概率 
figure(1); 
plot(pi); 
xlabel('场景数'); 
ylabel('概率');
title('各场景下概率');
%%
% %% 作图--光(原始数据)
% figure(2);
% for i=1:200
%     x=1:24;
%     y=i*ones(1,24);
%     z=P_pv1(:,i);
%     plot3(x,y,z);hold on;
%     grid on
%     xlabel('时间/h');
%     xlim([1 24])
%     ylabel('场景');
%     zlabel('光伏出力/kW');
%     title('200个场景下光伏出力');
% end
% 
% %% 作图--光(缩减后)
% figure(3);
% x=1:24;
% y1=ones(1,24);
% y2=2*ones(1,24);
% y3=3*ones(1,24);
% y4=4*ones(1,24);
% 
% p1=P_wind_pre(1,1:24);
% p2=P_wind_pre(2,1:24);
% p3=P_wind_pre(3,1:24);
% p4=P_wind_pre(4,1:24);
% 
% plot3(x,y1,p1);hold on;
% plot3(x,y2,p2);hold on;
% plot3(x,y3,p3);hold on;
% plot3(x,y4,p4);hold on;
% 
% grid on
% xlabel('时间/h');
% xlim([1 24])
% ylabel('场景');
% set(gca,'ytick',[1,2,3,4]);
% zlabel('光伏出力/kW');
% legend('场景一','场景二','场景三','场景四')
% title('缩减后4个场景下光伏出力');
% 
%% 作图--风(原始数据)
figure(2);
for i=1:200
    x=1:24;
    y=i*ones(1,24);
    z=P_wd1(:,i);
    plot3(x,y,z);hold on;
    grid on
    xlabel('时间/h');
    xlim([1 24])
    ylabel('场景');
    zlabel('风力发电/kW');
    title('200个场景下风力发电');
end
%% 作图--风(缩减后)
figure(3);
x=1:24;
y1=ones(1,24);
y2=2*ones(1,24);
y3=3*ones(1,24);
y4=4*ones(1,24);

p1=P_wind_pre(1,25:48);
p2=P_wind_pre(2,25:48);
p3=P_wind_pre(3,25:48);
p4=P_wind_pre(4,25:48);
P_wind=[p1;p2;p3;p4];

plot3(x,y1,p1);hold on;
plot3(x,y2,p2);hold on;
plot3(x,y3,p3);hold on;
plot3(x,y4,p4);hold on;

grid on
xlabel('时间/h');
xlim([1 24])
ylabel('场景');
set(gca,'ytick',[1,2,3,4]);
zlabel('风力发电/kW');
legend('场景一','场景二','场景三','场景四')
title('缩减后4个场景下风力发电');

% %% 作图--电负荷(原始数据)
% figure(6);
% for i=1:200
%     x=1:24;
%     y=i*ones(1,24);
%     z=R_fuhe(i,:);
%     plot3(x,y,z);hold on;
%     grid on
%     xlabel('时间/h');
%     xlim([1 24])
%     ylabel('场景');
%     zlabel('电负荷/kW');
%     title('200个场景下电负荷');
% end
% %% 作图--电负荷(缩减后)
% figure(7);
% x=1:24;
% y1=ones(1,24);
% y2=2*ones(1,24);
% y3=3*ones(1,24);
% y4=4*ones(1,24);
% 
% p1=P_wind_pre(1,49:72);
% p2=P_wind_pre(2,49:72);
% p3=P_wind_pre(3,49:72);
% p4=P_wind_pre(4,49:72);
% 
% plot3(x,y1,p1);hold on;
% plot3(x,y2,p2);hold on;
% plot3(x,y3,p3);hold on;
% plot3(x,y4,p4);hold on;
% 
% grid on
% xlabel('时间/h');
% xlim([1 24])
% ylabel('场景');
% set(gca,'ytick',[1,2,3,4]);
% zlabel('电负荷/kW');
% legend('场景一','场景二','场景三','场景四')
% title('缩减后4个场景下电负荷');

% %% 作图--热负荷(原始数据)
% figure(8);
% for i=1:200
%     x=1:24;
%     y=i*ones(1,24);
%     z=R_fuhe(i,:);
%     plot3(x,y,z);hold on;
%     grid on
%     xlabel('时间/h');
%     xlim([1 24])
%     ylabel('场景');
%     zlabel('热负荷/kW');
%     title('200个场景下热负荷');
% end
% %% 作图--热负荷(缩减后)
% figure(9);
% x=1:24;
% y1=ones(1,24);
% y2=2*ones(1,24);
% y3=3*ones(1,24);
% y4=4*ones(1,24);
% y5=5*ones(1,24);
% p1=P_wind_pre(1,73:96);
% p2=P_wind_pre(2,73:96);
% p3=P_wind_pre(3,73:96);
% p4=P_wind_pre(4,73:96);
% p5=P_wind_pre(5,73:96);
% plot3(x,y1,p1);hold on;
% plot3(x,y2,p2);hold on;
% plot3(x,y3,p3);hold on;
% plot3(x,y4,p4);hold on;
% plot3(x,y5,p5);hold on;
% grid on
% xlabel('时间/h');
% xlim([1 24])
% ylabel('场景');
% set(gca,'ytick',[1,2,3,4,5]);
% zlabel('热负荷/kW');
% legend('场景一','场景二','场景三','场景四','场景五')
% title('缩减后5个场景下热负荷');

% figure(10)
% plot(P_PV); 
% xlabel('时间/t'); 
% ylabel('功率/kW');
% title('光伏不确定性出力');
% figure(11)
% plot(P_WD); 
% xlabel('时间/t'); 
% ylabel('功率/kW');
% title('风电不确定性出力');

% figure(12)
% plot(OP_load_e); 
% xlabel('时间/t'); 
% ylabel('功率/kW');
% title('电负荷不确定性需求功率');
% figure(13)
% plot(OP_load_h); 
% xlabel('时间/t'); 
% ylabel('功率/kW');
% title('热负荷不确定性需求功率');