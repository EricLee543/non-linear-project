% 锦标赛选择
function [parent_selected] = tour_selection(population)
% 从population中随机选择两个个体进行比较，选择等级高（如果同等级选择聚集距离大的）
% 有放回抽样（可能会选到相同的个体）
[pop_size, distance]=size(population);
rank=distance-1;
candidate=[randperm(pop_size);randperm(pop_size)]';
[~,col] = size(population);
parent_selected = zeros(pop_size, col);
for i = 1:pop_size
    parent=candidate(i,:);
    if population(parent(1),rank)~=population(parent(2),rank)       % 个体等级不同，选择等级较高的
        if population(parent(1),rank)<population(parent(2),rank)
            parent_selected(i,:)=population(parent(1),:);
        else % pool(parent(1),rank)>pool(parent(2),rank)
            parent_selected(i,:)=population(parent(2),:);
        end
    else                                                % 个体等级相同，选择距离大的
        if population(parent(1),distance)>population(parent(2),distance)
            parent_selected(i,:)=population(parent(1),:);
        elseif population(parent(1),distance)< population(parent(2),distance)
            parent_selected(i,:)=population(parent(2),:);
        else                                            % 距离相同，随机选择一个
            temp=randperm(2);
            parent_selected(i,:)=population(parent(temp(1)),:);
        end
    end
end
end