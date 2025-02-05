function [V]=updata_v(wmax,wmin,index_i,maxIterations,sizepop,pop,v,pbest,gbest,vmax,dimpop,P_Li)
c1=((0.5-2.5)*index_i/maxIterations)+2.5;         %个体学习因子 随迭代次数增加，减小，防止过早期早收敛
c2=((2.5-0.5)*index_i/maxIterations)+0.5;         %群体学习因子 随迭代次数增加，增大，增加后期收敛速度
w=wmax-(wmax-wmin)*(index_i)^2/(maxIterations)^2; % 随迭代次数增加减少，减少自身在迭代次数的影响因素，加快后期收敛速度，与精度
% w=-w.*P_Li;
for index_j = 1:sizepop
    for index_k=1:dimpop
        %             %% 速度更新
        r1=2*rand(1)-1;                                             %-1到1随机值
        r2=2*rand(1)-1;
        %             %为增加计算速度，此处r1r2用之前的
        %             gailv=sign(((r1+r2)/4+0.5-0.2)-dijian*0.8);%最开始有80%的概率大于零，最后大于零的概率为0.
        %             dijian_k=((index_k-1)^2/(dimpop-1)^2);
        %             gailv_k=sign(((r1+r2)/4+0.5-0.2)-dijian_k*0.8);
        v(index_j,index_k) = ((w*v(index_j,index_k) + r1*c1*(pbest(index_j,index_k) - pop(index_j,index_k)) + r2*c2*(gbest(index_k) - pop(index_j,index_k))));
        %             if dis(index_k-1)*pop(index_j,index_k-1)>0||gailv>0||gailv_k>0||dis(index_k-1)*v(index_j,index_k-1)>0
        %                 v(index_j,index_k)=-v(index_j,index_k);
        %             end
        %% 限幅处理
        if(v(index_j,index_k)>vmax(1,index_k))
            v(index_j,index_k)=vmax(1,index_k);       %容量速度超上限
        elseif(v(index_j,index_k)<-vmax(1,index_k))
            v(index_j,index_k)=vmax(1,index_k);        %容量速度超下限
        end
    end
end
V=v;
end