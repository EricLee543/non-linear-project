% 2N中选出N个个体组成新种群
function new_pop=elitist_strategy(population_inter_sorted,front,pop_size)
index=0;
rank=1;
% 先按rank，再按拥挤度选择
while index < pop_size
    l_f=length(front(rank).fr);
    if index+l_f < pop_size 
        new_pop(index+1:index+l_f,:)= population_inter_sorted(index+1:index+l_f,:);
        index=index+l_f;
    else
        temp=population_inter_sorted(index+1:index+l_f,:);
        temp=sortrows(temp,size(temp,2));
        new_pop(index+1:pop_size,:)= temp(l_f-(pop_size-index)+1:l_f,:);
        index=index+l_f;
    end
    rank=rank+1;
end
end
