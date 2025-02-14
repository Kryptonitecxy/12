%运行成本和运行可靠性多目标
function [C1,C2,Cy,Cg,Cq,Cdan,Closs,Vp,E_soc,P_st,Inv1]=low_obj(x,upx)
parameter;
upx(1:end-1)=round(upx(1:end-1));
global P_wind pi;
%  预加载模型与配置***********************************
mpc_base = loadcase('case33bw');% 预加载配电网模型
pf_options = mpoption('OUT_ALL',0,'VERBOSE',0,'PF_ALG',2);
% 预先获取风机和储能的节点位置
wind_nodes = upx(1:Wd);
storage_nodes = upx(Wd+1:Wd+St);
%****************************************************
%运行成本=运维成本+购售电成本

%风机出力矩阵调整
wd1_s = reshape(x(1:K*Wd*T), [K*Wd, T]);%调整为K*Wd, T的矩阵
%****************************************************

%**************
Cy=0; %总运维成本
Cy_1=0; %风机运维成本
Cy_2=0; %储能运维成本
Cy_3=0; %并网逆变器运维成本
%购售电成本
Cgc=zeros(K,T);
%网损
Clossc=zeros(K,T);
%电压偏移量
Vpc=zeros(K,T);
%弃风惩罚
Cq=0;
%失电惩罚
Cdan=0;
%离子更新偏差
%p_lizi=zeros(Wd,K*T);
%P_lizi=zeros(1,Wd*K*T);
E_s=zeros(K, St*(T+1));
E_soc=zeros(1, K*St*(T+1));
P_s=zeros(K,T*St);%
P_st=zeros(1,St*T*K);%
inv2=zeros(K,T);
Inv1=zeros(1, K*T);
%****************************************************
% flag=1;
% for i=1:Wd
%     for j=1:K
%         Cy_1=Cy_1+crun_wd.*(pi(j).*sum(wd1_s(flag,:))); %crun_wd风机运维费用
%         flag=flag+1;
%     end
% end

% 计算每一行（每个风机在 T 个时间段的出力）的和
rowSums = sum(wd1_s, 2);  % 结果为 (K*Wd)×1 向量
% 构造每行对应的概率索引：对于 flag = 1,2,...,K*Wd, 对应的 j = mod(flag-1, K) + 1
flag_idx = mod((1:(K*Wd))-1, K) + 1;  % 1×(K*Wd) 行向量
% 按照权重乘以每行的和，并求总和
weightedSums = rowSums .* (pi(flag_idx))';  % (K*Wd)×1 向量
Cy_1 = crun_wd * sum(weightedSums);

%储能出力约束初始化***********************
%E_en=zeros(St,1); %储能额定容量
E_POW=zeros(St,1); %储能最大充放电出力(为储能额定容量的1/4)
E_SOC=zeros(St,25); %储能状态矩阵
P_St=zeros(St,24); %储能出力矩阵
Inv2=zeros(1,24); %并网逆变器出力矩阵
% for a=1:St
%    E_en(a,1)=upx(D+Wd+a).*s_st;
%    E_SOC(a,1)=0.5.*E_en(a,1); %初始SOC为50%       
%    E_POW(a,1)=0.25.*E_en(a,1);
% end
% 在循环外预计算电价时段掩码*******************
time_vector=1:24;
peak_hours = (time_vector >= 10 & time_vector <= 11) | (time_vector >= 19 & time_vector <= 21);
valley_hours = (time_vector >= 1 & time_vector <= 6) | (time_vector >= 23 & time_vector <= 24);
%********************************************
% ==== 4. 全局变量转局部变量 ====
local_pi = pi; % 解除global依赖
%local_P_wind = P_wind;
%********************************************
parfor i=1:K 
    % ==== 每个迭代的独立变量 ====
    local_Cgc = zeros(1,T);
    local_Clossc = zeros(1,T);
    local_Vpc = zeros(1,T);
 
  
    % ==== 场景独立变量初始化 ====
    local_E_SOC = zeros(St,25); % 储能SOC矩阵
    local_P_St = zeros(St,24);  % 储能出力矩阵
    local_Inv2 = zeros(1,24);   % 逆变器出力矩阵
    % ==== 储能参数初始化 ====
    E_en = upx(D+Wd+1:D+Wd+St).*s_st; % 直接向量化计算
    local_E_SOC(:,1) = 0.5 * E_en;    % 初始SOC
    E_POW = 0.25 * E_en;             % 最大充放电功率

    for t=1:T
       % s=0; 
        cx=mpc_base;
        cx.bus(:,3)=pl(t).*cx.bus(:,3); %cx.bus(:,3)有功负荷kw
        cx1=sum(cx.bus(:,3));%%%
         %计算系统功率盈余还是亏损  
        %s=sum(wd1_s((0: Wd-1)*K+i,t)/1000)-cx1;
        s = sum(x((i-1)*Wd*T + (1:Wd)*T - (T-t))) / 1000 - cx1; % 直接索引风机出力
       
        s=s*1000;
        %% **********储能约束**********%             
        %电池约束
        if s>0 %功率盈余,储能充电，无法充电，则向电网售电
              % 计算所有储能的可用充电容量（向量化）
              available_charge = socmax * E_en(:,1) - local_E_SOC(:,t); % St×1 向量
              % 总可用充电容量
              total_available = sum(available_charge);
              if total_available > 0
              % 计算每个储能的充电比例（归一化）(根据剩余可充电容量大小分配，剩余容量大的多分配功率)
                  charge_ratio = available_charge / total_available; % St×1 向量
                  % 计算理论充电功率（向量化）
                  P_St_theoretical = s * charge_ratio; % St×1 向量
                  % 限制最大充电功率（向量化）
                  local_P_St(:,t) = min(P_St_theoretical, E_POW(:,1)); % St×1 向量
                  local_E_SOC(:,t+1) = local_E_SOC(:,t) + eta * local_P_St(:,t);% 更新SOC（向量化）  
                  exceed_mask = local_E_SOC(:,t+1) > socmax * E_en(:,1);% 处理SOC越界（向量化）
                    if any(exceed_mask) 
                    local_P_St(exceed_mask,t) = (socmax * E_en(exceed_mask,1) - local_E_SOC(exceed_mask,t)) / eta;% 计算修正后的充电功率
                    local_E_SOC(exceed_mask,t+1) = socmax * E_en(exceed_mask,1);% 重新计算SOC
                    end
                  s = s - sum(local_P_St(:,t));% 更新剩余功率
              else % 无可用充电容量
                  local_P_St(:,t) = 0;
                  local_E_SOC(:,t+1) = local_E_SOC(:,t);
              end
        elseif s<0 %功率亏损，储能放电，无法放电，则向电网购电
              available_discharge = local_E_SOC(:,t) - socmin * E_en(:,1); % St×1 向量
              total_available = sum(available_discharge);
              if total_available > 0
                  discharge_ratio = available_discharge / total_available;
                  P_St_theoretical = s * discharge_ratio;
                  local_P_St(:,t) = max(P_St_theoretical, -E_POW(:,1));
                  local_E_SOC(:,t+1) = local_E_SOC(:,t) + local_P_St(:,t) / eta;
                  exceed_mask = local_E_SOC(:,t+1) < socmin * E_en(:,1);
                  if any(exceed_mask)
                      local_P_St(exceed_mask,t) = (local_E_SOC(exceed_mask,t) - socmin * E_en(exceed_mask,1)) * eta;
                      local_E_SOC(exceed_mask,t+1) = socmin * E_en(exceed_mask,1);
                  end
                   s = s - sum(local_P_St(:,t));
              else
                  local_P_St(:,t) = 0;
                  local_E_SOC(:,t+1) = local_E_SOC(:,t);
              end
        end
        Cy_2=Cy_2+crun_st*local_pi(i)*sum(abs(local_P_St(:,t))); %加入储能运维费用统计
        %%
        %更新选址节点处的负荷
        %cx.bus(wind_nodes,3)=cx.bus(wind_nodes,3)-wd1_s((0: Wd-1)*K+i,t)/1000;%upx(j1)是风机j1的选址点
        cx.bus(wind_nodes,3)=cx.bus(wind_nodes,3)- x((i-1)*Wd*T + (1:Wd)*T - (T-t))/1000;%upx(j1)是风机j1的选址点
        
        cx.bus(storage_nodes,3)=cx.bus(storage_nodes,3)+local_P_St(:,t)/1000;
       % disp(['启动第',num2str(t),'次潮流计算']);%调试用
        cc=runpf(cx,pf_options);%潮流计算
        
        %购售电成本(购电为负，售电为正)
        if s<=0 && s>=-upx(2*D+1)*econ %表示从主电网购电，有限制
            if peak_hours(t)
                local_Cgc(i,t)=s*pbuy_top*local_pi(i); %购电成本,不同时间段售价不同
            elseif valley_hours(t)
                local_Cgc(i,t)=s*pbuy_down*local_pi(i);
            else
                local_Cgc(i,t)=s*pbuy_ave*local_pi(i);
            end
            Cy_3=Cy_3+abs(s)*crun_con*local_pi(i); %加入并网逆变器运维费用
            local_Inv2(1,t)=-s;
            %s=0;
        elseif s>0 && s<=upx(2*D+1) %表示向主电网售电，有限制
            if peak_hours(t)
                local_Cgc(i,t)=s*econ*psell_top*local_pi(i); %售电成本
            elseif valley_hours(t)
                local_Cgc(i,t)=s*econ*psell_down*local_pi(i);
            else
                local_Cgc(i,t)=s*econ*psell_ave*local_pi(i);
            end
            Cy_3=Cy_3+abs(s)*crun_con*local_pi(i); %加入并网逆变器运维费用
            local_Inv2(1,t)=-s;
            %s=0;
        elseif s>upx(2*D+1) %超过并网逆变器功率
            if peak_hours(t)
                local_Cgc(i,t)=upx(2*D+1)*econ*pbuy_top*local_pi(i); %售电成本
            elseif valley_hours(t)
                local_Cgc(i,t)=upx(2*D+1)*econ*pbuy_down*local_pi(i);
            else
                local_Cgc(i,t)=upx(2*D+1)*econ*pbuy_ave*local_pi(i);
            end
            Cy_3=Cy_3+upx(2*D+1)*crun_con*local_pi(i); %加入并网逆变器运维费用
            local_Inv2(1,t)=-upx(2*D+1);
            s=(s-upx(2*D+1))/1000;
            Cq=Cq+cq_wd*abs(s)*abs(s)*local_pi(i)/cx1; %加入弃风惩罚成本
        else
            if peak_hours(t)
                local_Cgc(i,t)=-upx(2*D+1)*psell_top*local_pi(i); %购电成本
            elseif valley_hours(t)
                local_Cgc(i,t)=-upx(2*D+1)*psell_down*local_pi(i);
            else
                local_Cgc(i,t)=-upx(2*D+1)*psell_ave*local_pi(i);
            end
            Cy_3=Cy_3+upx(2*D+1)*crun_con*local_pi(i); %加入并网逆变器运维费用
            local_Inv2(1,t)=upx(2*D+1);
            s=(s+upx(2*D+1))/1000;
            Cdan=Cdan+cdan*abs(s)*abs(s)*local_pi(i)/cx1; %加入失电惩罚成本
        end
              
        %运行可靠性=网损＋电压偏移量+弃风惩罚+失电惩罚
        local_Clossc(i,t)=local_pi(i)*sum(cc.branch(:,14)+cc.branch(:,16)); %网损
        voltage_mean=mean(cc.bus(:,8));
        local_Vpc(i,t)=local_pi(i)*sum(abs(cc.bus(:,8)-voltage_mean)); %电压偏移量
        %local_p_lizi(:,t) = s/cx1;
    end
    %disp(['第',num2str(i),'个场景的24小时循环结束']);
    % for j=1:Wd
    %     p_lizi(j,(i-1)*24+1:(i-1)*24+T)=p_lizi(j,(i-1)*24+1:(i-1)*24+T)./max(p_lizi(j,(i-1)*24+1:(i-1)*24+T));
    % end
    Cgc(i,:) = local_Cgc;
    Clossc(i,:) = local_Clossc;
    Vpc(i,:) = local_Vpc;
    E_s(i,:)=reshape(local_E_SOC',[],1)'; %储能状态矩阵
    P_s(i,:)=reshape(local_P_St',[],1)'; %储能出力矩阵
    inv2(i,:)=local_Inv2;%并网逆变器出力矩阵
    
end
E_soc(1,:)=reshape(E_s',[],1)'; %拉直储能状态矩阵
P_st(1,:)=reshape(P_s',[],1)'; %拉直储能出力矩阵（充电为正，放电为负）
Inv1(1,:)=reshape(inv2',[],1)'; %拉直并网逆变器出力矩阵
%P_lizi(1,:)=reshape(p_lizi',[],1)';
Cy=Cy_1+Cy_2+Cy_3;
Cg=-sum(sum(Cgc));
Closs=sum(sum(Clossc));
Vp=sum(sum(Vpc));


C1=Cy+Cg; %总运行成本=总运维成本+购售电成本(购电-售电)
C2=Closs+Vp+Cq+Cdan; %网损+电压偏移量+弃风惩罚+失电惩罚
end





