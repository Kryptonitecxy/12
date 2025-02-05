function [result,C,gxbest1,Power_st,all_Cost]=uplayer(x) 

global gxbest1 Power_st all_Cost
parameter;
%上层目标：配置成本=投资成本+安装费用-剩余价值
C1=sum(x(D+1:D+Wd))*(cbuy_wd*(r*(1+r)^year1)/((1+r)^year1-1)+cins_wd-crem_wd); %风机配置成本
C2=sum(x(D+Wd+1:2*D))*(cbuy_st*(r*(1+r)^year2)/((1+r)^year2-1)+cins_st-crem_st); %储能配置成本
C3=x(2*D+1)*cbuy_con*(r*(1+r)^year3)/((1+r)^year3-1); %并网逆变器投资成本
C=C1+C2+C3; %总配置成本

%带入参数调用下层运行层程序
[gxbest1,Power_st,all_Cost,fymin1]=lowlayer(x); 
% disp('下层优化结束，进入上层：');
result=C+sum(fymin1);
% result=C;
end