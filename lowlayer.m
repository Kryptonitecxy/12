function [gxbest1,Power_st,all_Cost,fymin1]=lowlayer(upx) %下层多目标优化，包括运行成本、电压偏移量

global P_wind;
parameter;
gama=0.1; %风机出力下限，上限为1

maxFun=2; %2个目标函数
fff=[0,800000;0,90000]; %各个目标函数的最小值和最大值，即绝对正理想解和绝对负理想解
n=50; %初始种群个数
d=96*Wd; %决策变量个数(4个场景风机出力和储能出力)
maxIterations=50; %最大迭代次数
wmax=0.9; %maximum of inertia factor，最大惯性系数
wmin=0.4; %minimum of inertia factor，最小惯性系数
Xmax=zeros(1,96*Wd);
Xmin=zeros(1,96*Wd);

%运行阶段4个场景风机和储能出力的决策空间设置
for i=1:Wd
    Xmax(1,96*(i-1)+1:96*(i-1)+96)=[upx(D+i).*P_wind(1,:),upx(D+i).*P_wind(2,:),upx(D+i).*P_wind(3,:),upx(D+i).*P_wind(4,:)];
    Xmin(1,96*(i-1)+1:96*(i-1)+96)=[gama.*upx(D+i).*P_wind(1,:),gama.*upx(D+i).*P_wind(2,:),gama.*upx(D+i).*P_wind(3,:),gama.*upx(D+i).*P_wind(4,:)];
end
%储能的运行策略在储能约束中，此处只需要决策风机的出力

dX=Xmax-Xmin;
Vmax=dX;
Vmin=ones(1,96*Wd);

%种群位置与速度初始化
X=repmat(Xmin,n,1)+repmat((Xmax-Xmin),n,1).*rand(n,d); %初始种群的位置 
V=repmat(Vmin,n,1)+repmat((Vmax-Vmin),n,1).*rand(n,d); % 初始种群的速度

%变量维数(总共96*Wd)
x=X;
v=V;

%个体及种群历史最优初始化
for iterations=1:maxIterations %达到迭代次数退出循环
    for m=1:maxFun
        fxmin(iterations,m)=inf;
    end
end
for b=1:n                 %第b个粒子
    px(b,:)=x(b,:);       %初始化个体的历史最优位置px，n*d矩阵，d---粒子维数
    for m=1:maxFun        %第m个目标函数
        fxm(b,m)=inf;     %初始化个体的历史最佳适应度，因求最小值，初代设为无穷大
    end
end
for m=1:maxFun
    fym(1,m)=inf;        %初始化群体的历史最佳适应度
end

%外部归档集初始化
archive=[]; %pareto临时最优断面在种群更新完成之后，MOPSO进行了三轮筛选。
%首先，根据支配关系进行第一轮筛选，将劣解去除，剩下的加入到存档中。
%其次，在存档中根据支配关系进行第二轮筛选，将劣解去除，并计算存档粒子在网格中的位置。
%最后，若存档数量超过了存档阀值，则根据自适应网格进行筛选，直到阀值限额为止。重新进行网格划分
PP=[]; %运行中的出力等信息
Vall=[]; %运行中的成本等信息

%群体更新
iterations=1;
% disp('下层优化开始：');
while(iterations<=maxIterations)
    %计算种群适应度
    disp(['下层迭代次数：',num2str(iterations)]);
    for b=1:n
        [fx(b,1),fx(b,2),Cy,Cg,Cq,Cdan,Closs,Vp,E_soc,P_st,Inv1,P_lizi]=low_obj(x(b,:),upx);
        P(b,:)=[E_soc P_st Inv1 fx(b,1) fx(b,2) zeros(1,1)]; %保存运行信息
        Val(b,:)=[Cy Cg Cq Cdan Closs Vp fx(b,1) fx(b,2) zeros(1,1)];
        P_Li(b,:)=P_lizi;
    end

    %将粒子的“位置”信息与“适应度f1，f2”信息结合，生成Stem1矩阵
    %前d列存放m个粒子的位置信息，第d+1列、d+2列分别对应存放双目标的适应度信息f1，f2
    stem1=x;                         %存放位置信息
    stem1(:,d+1)=fx(:,1);            %存放适应度信息f1
    stem1(:,d+2)=fx(:,2);            %存放适应度信息f2
    %--update px--%更新个体极值
    for b=1:n
        if((stem1(b,d+1)<=fxm(b,1))&&(stem1(b,d+2)<=fxm(b,2))&&(~((stem1(b,d+1)==fxm(b,1))&&(stem1(b,d+2)==fxm(b,2)))))
            %--new point dominate old pbest--%新粒子支配历史最佳粒子(新粒子较优)
            fxm(b,1)=stem1(b,d+1);
            fxm(b,2)=stem1(b,d+2);
            px(b,:)=stem1(b,1:d);
        elseif((stem1(b,d+1)>=fxm(b,1))&&(stem1(b,d+2)>=fxm(b,2))&&(~((stem1(b,d+1)==fxm(b,1))&&(stem1(b,d+2)==fxm(b,2)))))
            %--old pbest dominate new point--%历史最佳粒子支配新粒子
        else
            %--non-dominate between new point and old pbest--%
            if(rand>0.5)
                fxm(b,1)=stem1(b,d+1);
                fxm(b,2)=stem1(b,d+2);
                px(b,:)=stem1(b,1:d);
            end
            %--注：还可以通过拥挤距离机制选择新的pbest--%
        end
    end
    
%         plot(fx(:,1),fx(:,2),'*');
%  xlabel('目标1');ylabel('目标2');
%  grid on;
%  title('帕累托解集');
% drawnow
    %--update px finished--%
    
     %--update gx--%更新群体极值
%-----外部归档集初始化；外部归档集以Stem1的形式保存粒子信息（位置；适应度；拥挤距离信息）------%
    for b=1:n  %每次只考虑将一个粒子加到外部归档集
        insert_flag=5; %设标记位初始值，目的：每次循环重置标记位的值
        if(isempty(archive))        %第一个粒子直接存入archive
            if stem1(b,d+1)~=inf
                archive(1,1:d+2)=stem1(b,:);
                PP(1,:)=P(b,:);
                Vall(1,:)=Val(b,:);
            end
        else                        %第二个粒子开始，通过和archive中的所有粒子进行比较，进而判断是否将粒子存入archive
        %%--for循环作用：①判断第i个粒子是否可加入档案（“支配”或“非支配”所有archive中粒子） ②标记archive中被粒子i支配的粒子，准备删除
            for k=1:size(archive,1) % 返回archive行数，即粒子数。依次与archive中的所有粒子进行比较，
                if((stem1(b,d+1)<=archive(k,d+1))&&(stem1(b,d+2)<=archive(k,d+2))&&(~((stem1(b,d+1)==archive(k,d+1))&&(stem1(b,d+2)==archive(k,d+2)))))
                    %--i dominate k--%待加入粒子i支配k
                    %archive(k,1)=1000000;%mark delete k，标记待删除粒子
                    archive(k,d+1)=inf;%mark delete k，标记待删除粒子
                    insert_flag=1; 
                elseif((stem1(b,d+1)>=archive(k,d+1))&&(stem1(b,d+2)>=archive(k,d+2))&&(~((stem1(b,d+1)==archive(k,d+1))&&(stem1(b,d+2)==archive(k,d+2)))))
                    %--k dominate i--%待加入粒子i被k支配
                    insert_flag=0;     %粒子i被k支配则不可能支配pareto解集中的其他粒子
                    break     % 跳出for循环
                else
                        %--non-dominate between i,k--%待加入粒子与k无支配关系
                    insert_flag=1;
                end
            end
       %%----for循环结束------------------------------------%%
       %----搜索archive中被标记的粒子，并删除----%
            k1=1;
            while(k1<=size(archive,1))
                if(archive(k1,d+1)==inf)
                    archive(k1,:)=[]; %删除此行
                    PP(k1,:)=[];
                    Vall(k1,:)=[];
                    k1=k1-1;
                end
                k1=k1+1;
            end
        %--删除结束--%
      
        %--添加第i个粒子到archive中--%
           if(insert_flag==1)
               if stem1(b,d+1)~=inf
                   archive((size(archive,1)+1),1:d+2)=stem1(b,:);
                   PP((size(PP,1)+1),:)=P(b,:);
                   Vall((size(Vall,1)+1),:)=Val(b,:);
               end
           end
        %--添加第i个粒子结束--%
        end
    end
    %-----外部归档集初始化结束-------%
    if(isempty(archive))    %空数组返回逻辑1
        display('程序出错');
        pause;
    end
    
    %-----------计算archive中各粒子的拥挤距离-------------%
    %--按照拥挤距离从大到小的顺序排序--%
    %limite archive's size%
    archive(:,d+5)=zeros(size(archive,1),1);
    for m=1:maxFun
        archive=sortrows(archive,d+m); %升序排序，对第两个目标函数排序
        PP=sortrows(PP,size(PP,2)-3+m);
        Vall=sortrows(Vall,size(Vall,2)-3+m);
        fmin=archive(1,d+m); %保留最小的一组
        archive(1,d+5)=inf;
        fmax=archive(size(archive,1),d+m); %保留最大的一组
        archive(size(archive,1),d+5)=inf;
        if size(archive,1)>2
            for i=2:size(archive,1)-1
                archive(i,d+5)=archive(i,d+5)+(archive(i+1,d+m)-archive(i-1,d+m))/(fmax-fmin); %计算拥挤距离
            end
        end
    end
    PP(:,size(PP,2))=archive(:,d+5);
    Vall(:,size(Vall,2))=archive(:,d+5);
    archive=sortrows(archive,-(d+5)); %以拥挤距离降序排序
    PP=sortrows(PP,-size(PP,2));
    Vall=sortrows(Vall,-size(Vall,2));
    if size(archive,1)>100 %控制pareto解集规模
        archive=archive(1:100,:);
        PP=PP(1:100,:);
        Vall=Vall(1:100,:);
    end
    %向上取整
    ddy=ceil(size(archive,1)*0.1*rand); %选取拥挤距离前10%里的任意一组作为最优粒子
    gx=archive(ddy,1:d);
    yy_best(iterations,:)=gx;
    yy_fitness(iterations,:)=archive(ddy,d+1:d+2);
    archive(:,d+3:d+4)=zeros(size(archive,1),2);
    for m=1:maxFun
        archive(:,d+2+m)=(archive(:,d+m)-fff(m,1))/(fff(m,2)-fff(m,1)); %目标函数无量纲化
    end
    %按照拥挤距离从大到小的顺序排序结束
    if iterations==maxIterations
        %基于信息熵确立权重的TOPSIS法开始
        Ad=archive(:,d+1:d+maxFun); %目标函数
        AA=Ad;
        zhi=0;
        for m=1:maxFun
            Ad(:,m)=(Ad(:,m)-fff(m,1))/(fff(m,2)-fff(m,1)); %目标函数无量纲化
            AA(:,m)=Ad(:,m)/sum(Ad(:,m));
        end
        for m=1:maxFun
            for i=1:size(AA,1)
                lAA(i,m)=AA(i,m).*log(AA(i,m));
                if AA(i,m)==0
                   lAA(i,m)=0;
                end
            end
            shang(1,m)=-1/log(maxFun)*sum(lAA(:,m));
            zhi=zhi+(1-shang(1,m));
        end
        for m=1:maxFun
            quan(1,m)=(1-shang(1,m))/zhi;
        end
        B=ones(1,2); %目标函数参考值，最大，即绝对负理想解
        C=zeros(1,2); %目标函数参考值，最大，即绝对正理想解
        SF1=sqrt(sum(bsxfun(@times,bsxfun(@minus,Ad,B),quan).*bsxfun(@times,bsxfun(@minus,Ad,B),quan),2));
        SG1=sqrt(sum(bsxfun(@times,bsxfun(@minus,Ad,C),quan).*bsxfun(@times,bsxfun(@minus,Ad,C),quan),2));
        SC1=SF1./(SF1+SG1);  %贴合度
        [s1,dy1]=max(SC1);
        gxbest1=archive(dy1,1:d); %种群最优
        Power_st=PP(dy1,1:end-1);
        all_Cost=Vall(dy1,1:end-1);
        fymin1=archive(dy1,d+1:d+maxFun);
        %基于信息熵确立权重的TOPSIS法结束
    else
        %1.惯性函数最大值2.惯性函数最小值，3.当前迭代次数4.总迭代次数5.粒子个数6.粒子7.速度.8.pb 9.gb 10
        %粒子的速度更新
        v=updata_v(wmax,wmin,iterations,maxIterations,n,x,v,px,gx,Vmax,d,P_Li);
        %粒子的位置更新
        x=x+v;       
        x=pop_limit(x,n,Xmax,Xmin,d,v);
    end
    iterations=iterations+1;
end
