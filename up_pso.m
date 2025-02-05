function [uu,pg,c,g,po,ac]=up_pso(x)
global gxbest1 Power_st all_Cost
parameter; %基本参数定义
%粒子群参数初始化
Max_Dt=1; %最大迭代次数
w_max=0.9; %最大惯性系数
w_min=0.4; %最小惯性系数
v_max=2; 
p=[]; %保存每个粒子的适应度值
Pbest=inf; %求最小适应度值，初代为无穷大
v=zeros(N,2*D+1);
c=0; %记录最优配置成本
g=[]; %记录运行层最优风机出力
po=[]; %记录运行层储能信息
ac=[]; %记录运行层的各类成本

%计算各个粒子的适应度，并初始化Pi和Pg****************
for i=1:N
    p(i)=inf; %求最小适应度值，初代为无穷大
    y(i,:)=x(i,:); %每个粒子的个体寻优值
end
pg=x(1,:); %Pg为全局最优
uu=[];

%进入主循环*****************************************
for t=1:Max_Dt
    disp(['上层迭代次数：',num2str(t)]);
    for i=1:N
%         disp(['粒子：',num2str(i)]);
        %计算粒子的适应度
        [H(i),C,gxbest1,Power_st,all_Cost]=uplayer(x(i,:));
        if H(i)<p(i)
            p(i)=H(i);
            y(i,:)=x(i,:);
        end
        if p(i)<Pbest %记录全局最优的数据
            Pbest=p(i);
            pg=y(i,:);
            c=C;
            g=gxbest1;
            po=Power_st;
            ac=all_Cost;
        end
        w=w_max-(w_max-w_min)*t/Max_Dt; %惯性权重更新
        c1=(0.5-2.5)*t/Max_Dt+2.5; %自我学习因子
        c2=(2.5-0.5)*t/Max_Dt+0.5; %群体学习因子     
        v(i,:)=w*v(i,:)+c1*rand()*(y(i,:)-x(i,:))+c2*rand()*(pg-x(i,:));
        for m=1:2.*D+1
            if(v(i,m)>v_max)
                v(i,m)=v_max;
            elseif(v(i,m)<-v_max)
                v(i,m)=-v_max;
            end
        end
        x(i,:)=x(i,:)+v(i,:);
        x(i,1:size(x,2)-1)=round(x(1,1:size(x,2)-1));
        %处理粒子边界(即每个变量的范围)
        for n=1:2.*D+1
            if n<(Wd+1)
                if x(i,n)>max(ad_wd)
                     x(i,n)=max(ad_wd);
                     v(i,n)=-v(i,n); 
                elseif x(i,n)<min(ad_wd)
                     x(i,n)=min(ad_wd);
                     v(i,n)=-v(i,n); 
                end
            elseif n>=(Wd+1) && n<(Wd+St+1)
                if x(i,n)>max(ad_st)
                     x(i,n)=max(ad_st);
                     v(i,n)=-v(i,n); 
                elseif x(i,n)<min(ad_st)
                     x(i,n)=min(ad_st);
                     v(i,n)=-v(i,n); 
                end
            elseif n>=(D+1) && n<(D+Wd+1)
                if x(i,n)>max(num_wd)
                     x(i,n)=max(num_wd);
                     v(i,n)=-v(i,n); 
                elseif x(i,n)<min(num_wd)
                     x(i,n)=min(num_wd);
                     v(i,n)=-v(i,n); 
                end
            elseif n>=(D+Wd+1) && n<(2*D+1)
                if x(i,n)>max(num_st)
                     x(i,n)=max(num_st);
                     v(i,n)=-v(i,n); 
                elseif x(i,n)<min(num_st)
                     x(i,n)=min(num_st);
                     v(i,n)=-v(i,n); 
                end
            else
                if x(i,n)>1200
                   x(i,n)=1200;
                     v(i,n)=-v(i,n); 
                elseif x(i,n)<0
                     x(i,n)=0;
                     v(i,n)=-v(i,n); 
                end 
            end
        end
    end
    uu(t,1)=Pbest; %绘制迭代曲线
    uu(t,2)=sum(ac(1:4));
    uu(t,3)=sum(ac(5:6));
    uu(t,4)=c;
end 
end
 
 
 
 
 