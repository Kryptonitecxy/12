%运行成本和运行可靠性多目标
function [C1,C2,Cy,Cg,Cq,Cdan,Closs,Vp,E_soc,P_st,Inv1,P_lizi]=low_obj(x,upx)
parameter;
upx(1:end-1)=round(upx(1:end-1));
global P_wind pi;
%运行成本=运维成本+购售电成本
%风机出力矩阵调整
wd1_s=zeros(K*Wd,T);
for i=1:K*Wd
    wd1_s(i,:)=x(1+T*(i-1):T*i);
end

Cy=0; %总运维成本
Cy_1=0; %风机运维成本
Cy_2=0; %储能运维成本
Cy_3=0; %并网逆变器运维成本
E_s=zeros();
E_soc=zeros();
P_s=zeros();
P_st=zeros();
Inv1=zeros();
Inv2=zeros();
flag=1;
for i=1:Wd
    for j=1:K
        Cy_1=Cy_1+crun_wd.*(pi(j).*sum(wd1_s(flag,:)));
        flag=flag+1;
    end
end

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
p_lizi=zeros(Wd,K*T);
P_lizi=zeros(1,Wd*K*T);

for i=1:K
    %储能出力约束初始化
    E_en=zeros(St,1); %储能额定容量
    E_POW=zeros(St,1); %储能最大充放电出力(为储能额定容量的1/4)
    E_SOC=zeros(St,25); %储能状态矩阵
    P_St=zeros(St,24); %储能出力矩阵
    Inv2=zeros(1,24); %并网逆变器出力矩阵
    for a=1:St
        E_en(a,1)=upx(D+Wd+a).*s_st;
        E_SOC(a,1)=0.5.*E_en(a,1); %初始SOC为50%
        E_POW(a,1)=0.25.*E_en(a,1);
    end
    for t=1:T
        s=0; 
        cx=loadcase('case33bw');
        cx.bus(:,3)=pl(t).*cx.bus(:,3); %负荷
        cx1=sum(pl(t).*cx.bus(:,3));
        for j=1:Wd
            s=s+wd1_s(i+4*(j-1),t)/1000; 
        end
        s=s-sum(cx.bus(:,3)); %计算系统功率盈余还是亏损
        for j=1:Wd
            p_lizi(j,(i-1)*24+t)=s/sum(cx.bus(:,3));
        end
        s=s*1000;
        %% **********储能约束**********%             
        %电池约束
        if s>0 %功率盈余,储能充电，无法充电，则向电网售电
            for b=1:St
                if sum(socmax.*E_en(b:end,1)-E_SOC(b:end,t)) > 0 %有可用于充电的容量
                    E_re=(socmax.*E_en(b,1)-E_SOC(b,t))/sum(socmax.*E_en(b:end,1)-E_SOC(b:end,t)); %剩余可充电容量比例
                    P_St(b,t)=s.*E_re; %计算充电的功率(根据剩余可充电容量大小分配，剩余容量大的多分配功率)
                else
                    P_St(b,t)=0;
                    E_SOC(b,t+1)=E_SOC(b,t);
                    continue;
                end
                if P_St(b,t)>E_POW(b,1) %最大出力限制
                    P_St(b,t)=E_POW(b,1);
                end
                E_SOC(b,t+1)=E_SOC(b,t)+eta.*P_St(b,t); %计算SOC
                if E_SOC(b,t+1)>=socmax*E_en(b,1) %SOC限额(<=90%)
                    P_St(b,t)=-(E_SOC(b,t)-socmax*E_en(b,1))/eta; 
                    E_SOC(b,t+1)=E_SOC(b,t)+eta.*P_St(b,t);
                end
                s=s-P_St(b,t);
            end
        elseif s<0 %功率亏损，储能放电，无法放电，则向电网购电
            for b=1:St
                if sum(E_SOC(b:end,t)-socmin.*E_en(b:end,1)) > 0 %有可用于放电的容量
                    E_re=(E_SOC(b,t)-socmin.*E_en(b,1))/sum(E_SOC(b:end,t)-socmin.*E_en(b:end,1)); %剩余可放电容量比例
                    P_St(b,t)=s.*E_re; %计算充电的功率(根据剩余可放电容量大小分配，剩余容量大的多功率)
                else
                    P_St(b,t)=0;
                    E_SOC(b,t+1)=E_SOC(b,t);
                    continue;
                end
                if abs(P_St(b,t))>E_POW(b,1) %最大出力限制
                    P_St(b,t)=-E_POW(b,1);
                end
                E_SOC(b,t+1)=E_SOC(b,t)+P_St(b,t)./eta; %计算SOC
                if E_SOC(b,t+1)<=socmin*E_en(b,1) %SOC限额(<=90%)
                    P_St(b,t)=-(E_SOC(b,t)-socmin*E_en(b,1))*eta; 
                    E_SOC(b,t+1)=E_SOC(b,t)+P_St(b,t)./eta;
                end
                s=s-P_St(b,t);
            end
        end
        Cy_2=Cy_2+crun_st*pi(i)*sum(abs(P_St(:,t))); %加入储能运维费用统计
        %%
        for j1=1:Wd
            cx.bus(upx(j1),3)=cx.bus(upx(j1),3)-wd1_s((j1-1)*4+i,t)/1000;
        end
        for j2=1:St
            cx.bus(upx(Wd+j2),3)=cx.bus(upx(Wd+j2),3)+P_St(j2,t)/1000;
        end
        cc=runpf(cx,mpoption('OUT_ALL',0,'VERBOSE',0,'PF_ALG',3));
        
        %购售电成本(购电为负，售电为正)
        if s<=0 && s>=-upx(2*D+1)*econ %表示从主电网购电，有限制
            if (t>=10 && t<=11) || (t>=19 && t<=21)
                Cgc(i,t)=s*pbuy_top*pi(i); %购电成本
            elseif (t>=1 && t<=6) || (t>=23 && t<=24)
                Cgc(i,t)=s*pbuy_down*pi(i);
            else
                Cgc(i,t)=s*pbuy_ave*pi(i);
            end
            Cy_3=Cy_3+abs(s)*crun_con*pi(i); %加入并网逆变器运维费用
            Inv2(1,t)=-s;
            s=0;
        elseif s>0 && s<=upx(2*D+1) %表示向主电网售电，有限制
            if (t>=10 && t<=11) || (t>=19 && t<=21)
                Cgc(i,t)=s*econ*psell_top*pi(i); %售电成本
            elseif (t>=1 && t<=6) || (t>=23 && t<=24)
                Cgc(i,t)=s*econ*psell_down*pi(i);
            else
                Cgc(i,t)=s*econ*psell_ave*pi(i);
            end
            Cy_3=Cy_3+abs(s)*crun_con*pi(i); %加入并网逆变器运维费用
            Inv2(1,t)=-s;
            s=0;
        elseif s>upx(2*D+1) %超过并网逆变器功率
            if (t>=10 && t<=11) || (t>=19 && t<=21)
                Cgc(i,t)=upx(2*D+1)*econ*pbuy_top*pi(i); %售电成本
            elseif (t>=1 && t<=6) || (t>=23 && t<=24)
                Cgc(i,t)=upx(2*D+1)*econ*pbuy_down*pi(i);
            else
                Cgc(i,t)=upx(2*D+1)*econ*pbuy_ave*pi(i);
            end
            Cy_3=Cy_3+upx(2*D+1)*crun_con*pi(i); %加入并网逆变器运维费用
            Inv2(1,t)=-upx(2*D+1);
            s=(s-upx(2*D+1))/1000;
            Cq=Cq+cq_wd*abs(s)*abs(s)*abs(s)*pi(i)/cx1; %加入弃风惩罚成本
        else
            if (t>=10 && t<=11) || (t>=19 && t<=21)
                Cgc(i,t)=-upx(2*D+1)*psell_top*pi(i); %购电成本
            elseif (t>=1 && t<=6) || (t>=23 && t<=24)
                Cgc(i,t)=-upx(2*D+1)*psell_down*pi(i);
            else
                Cgc(i,t)=-upx(2*D+1)*psell_ave*pi(i);
            end
            Cy_3=Cy_3+upx(2*D+1)*crun_con*pi(i); %加入并网逆变器运维费用
            Inv2(1,t)=upx(2*D+1);
            s=(s+upx(2*D+1))/1000;
            Cdan=Cdan+cdan*abs(s)*abs(s)*abs(s)*pi(i)/cx1; %加入失电惩罚成本
        end
              
        %运行可靠性=网损＋电压偏移量+弃风惩罚+失电惩罚
        Clossc(i,t)=pi(i)*sum(cc.branch(:,14)+cc.branch(:,16)); %网损
        Vpc(i,t)=pi(i)*sum(abs(cc.bus(:,8)-mean(cc.bus(:,8)))); %电压偏移量
    end
    E_s(i,1:numel(E_SOC))=reshape(E_SOC',[],1)'; %储能状态矩阵
    P_s(i,1:numel(P_St))=reshape(P_St',[],1)'; %储能出力矩阵
    inv2(i,1:1:numel(Inv2))=reshape(Inv2',[],1)';
    for j=1:Wd
        p_lizi(j,(i-1)*24+1:(i-1)*24+T)=p_lizi(j,(i-1)*24+1:(i-1)*24+T)./max(p_lizi(j,(i-1)*24+1:(i-1)*24+T));
    end
end
E_soc(1,1:numel(E_s))=reshape(E_s',[],1)'; %拉直储能状态矩阵
P_st(1,1:numel(P_s))=reshape(P_s',[],1)'; %拉直储能出力矩阵（充电为正，放电为负）
Inv1(1,1:numel(inv2))=reshape(inv2',[],1)'; %拉直并网逆变器出力矩阵
P_lizi(1,1:numel(p_lizi))=reshape(p_lizi',[],1)';
Cy=Cy_1+Cy_2+Cy_3;
Cg=-sum(sum(Cgc));
Closs=sum(sum(Clossc));
Vp=sum(sum(Vpc));
% Cq=Cq/100;
% Cdan=Cdan/100;
C1=Cy+Cg; %总运行成本=总运维成本+购售电成本(购电-售电)
C2=Closs+Vp+Cq+Cdan; %网损+电压偏移量+弃风惩罚+失电惩罚
end





