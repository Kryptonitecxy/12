function [pop]=pop_limit(pop,sizepop,xmax,xmin,dimpop,v)
parameter;
for index_j=1:sizepop
    for index_k=1:dimpop
        if pop(index_j,index_k)>xmax(1,index_k)
            pop(index_j,index_k)=xmax(1,index_k) ;%容量位置选择超上限
            v(index_j,index_k)=-v(index_j,index_k);
        elseif pop(index_j,index_k)<xmin(1,index_k)
            pop(index_j,index_k)=xmin(1,index_k) ;%容量位置超下限
            v(index_j,index_k)=-v(index_j,index_k);
        end
    end
end
end

