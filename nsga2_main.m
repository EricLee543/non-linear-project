function nsga2_main()
pop_size = 50;  % 种群数量
V = 3;          % (x1,x2...xv)的维度，若翼形为16维描述符V=16
M = 3;          % 目标函数个数
no_runs = 5;    % 过程循环次数
gen_max = 500;    % 最大迭代次数
lb = [0.5 -2 -2]; %初始化的上下界
ub = [1 2 2];
p = parpool('local',8); %cpu核数，主要对 no_run进行并行计算
combined_paretoset=[];% 保存帕累托前沿
parfor run=1:no_runs
    x = lb+((ub-lb).*rand(pop_size, V)); %初始化种群
    [~, col_num] = size(x);
    f_val = zeros(pop_size,col_num);
    for i =1:pop_size
        [f_val(i,:) ,err(i,:)] =f_t(x(i,:));
    end
    error_norm=normalisation(err);
    population_init=[x f_val error_norm];
    [population, ~]=NDS_CD_cons(population_init,V,M); 
    for gen_cnt=1:gen_max
        % selection (Parent Pt of 'N' pop size)
        parent_selected=tour_selection(population);  
        % 生成子代
        child_offspring  = crossover(parent_selected(:,1:V),lb,ub,V);
        % 计算子代value
        [~, col_num] = size(child_offspring);
        child_f_val = zeros(pop_size,col_num);
        err = zeros(pop_size,col_num);
        for j = 1:pop_size
            [child_f_val(j,:) ,err(j,:)]=f_t(child_offspring(j,:));
        end        
        error_norm=normalisation(err);                                  
        child_offspring=[child_offspring child_f_val error_norm];
        
        % 子代与父代合并
        population_inter=[population(:,1:V+M+1) ; child_offspring(:,1:V+M+1)];
        [population_inter_sorted,front]=NDS_CD_cons(population_inter,V,M);
        
        % 精英保留策略选出N个个体，组成新一代种群
        new_pop=elitist_strategy(population_inter_sorted,front,pop_size);
        population=new_pop;
    end
    new_pop=sortrows(new_pop,V+1);
    %paretoset(run).trial=new_pop(:,1:V+M+1);
    combined_paretoset = [combined_paretoset;new_pop(:,1:V+M+1)];
end

if no_runs==1
    plot3(combined_paretoset(:,V+1),combined_paretoset(:,V+2),combined_paretoset(:,V+3),'*')
    grid on
else
    [pareto_filter,~]=NDS_CD_cons(combined_paretoset,V,M);               % Applying non domination sorting on the combined Pareto solution set
    rank1_index=find(pareto_filter(:,V+M+2)==1);        % Filtering the best solutions of rank 1 Pareto
    pareto_rank1=pareto_filter(rank1_index,1:V+M);
    plot3(pareto_rank1(:,V+1),pareto_rank1(:,V+2),pareto_rank1(:,V+3),'*')   % Final Pareto plot
    grid on
end
xlabel('objective function 1')
ylabel('objective function 2')
zlabel('objective function 3')
delete(p)
end


function err_norm  = normalisation(error_pop)
[N,nc]=size(error_pop);
con_max=0.001+max(error_pop);
con_maxx=repmat(con_max,N,1);
cc=error_pop./con_maxx;
err_norm=sum(cc,2);                % finally sum up all violations
end

function [y, err] = f_t(x)
%返回 多目标函数值y
%error的size为 [pop_size, num_of_cons]
%此处num_of_cons 为3，因为x1,x2,x3分别有一个限制条件
y(1) = x(1);
y(2) = 1/x(1)*(1+(x(2).^2+x(3).^2)^(0.25)*(sin(50*(x(2)^2+x(3)^2)^0.1).^2+1));
y(3) = x(2);
if x(1) <=1 && x(1) >= 0.5
    err1 = 0;
else
    err1 = min(abs(x(1)-0.5), abs(x(1)-1));
end

if x(2) <=2 && x(2) >= -2
    err2 = 0;
else
    err2 = min(abs(x(2)+2), abs(x(2)-2));
end

if x(3) <=2 && x(3) >= -2
    err3 = 0;
else
    err3 = min(abs(x(3)+2), abs(x(3)-2));
end
err = [err1 err2 err3];
end
