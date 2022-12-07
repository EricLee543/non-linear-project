function [chromosome_NDS_CD,front] = NDS_CD_cons(population,V,M)
problem_type  = 0;

chromosome_NDS_CD1=[];
infpop=[];
front.fr=[];
struct.sp=[];
rank=1;

if all(population(:,V+M+1)==0)
    problem_type=0;
    chromosome=population(:,1:V+M);                         % All Feasible chromosomes;    
    pop_size1=size(chromosome,1);
elseif all(population(:,V+M+1)~=0)
    problem_type=1;
    pop_size1=0;
    infchromosome=population;                               % All InFeasible chromosomes;       
else
    problem_type=0.5;
    feas_index=find(population(:,V+M+1)==0);
    chromosome=population(feas_index,1:V+M);                % Feasible chromosomes;    
    pop_size1=size(chromosome,1);
    infeas_index=find(population(:,V+M+1)~=0);
    infchromosome=population(infeas_index,1:V+M+1);         % infeasible chromosomes;    
end

rank = 1;
if problem_type==0 || problem_type==0.5
    pop_size1 = size(chromosome,1);
    f = chromosome(:,V+1:V+M);

    n = zeros(1,pop_size1);
    for p=1:pop_size1
        struct(p).sp = [];
        for q=1:pop_size1
            p_ind = f(p,:);
            q_ind = f(q,:);
            domination = judge_dominate(p_ind, q_ind);
            if domination == 1
                struct(p).sp = [struct(p).sp ;q];
            elseif domination == -1
                n(p) = n(p) + 1;
            end
        end
    end

    front(1).fr=find(n==0);

    % 完成quick sort
    while (~isempty(front(rank).fr))
        front_indiv=front(rank).fr;
        n(front_indiv)=inf;
        chromosome(front_indiv,V+M+1)=rank;
        rank=rank+1;
        front(rank).fr=[];
        for i = 1:length(front_indiv)
            temp=struct(front_indiv(i)).sp;
            n(temp)=n(temp)-1;
        end 
        q=find(n==0);
        front(rank).fr=[front(rank).fr q];
    end
    
    chromosome_sorted = sortrows(chromosome,V+M+1);%V+M+1 列存储rank
    rowsindex = 1; %从chromosome的第一行开始计算
    %为每一个层计算拥挤距离    
    for i=1:length(front)-1
        %计算每一层rank的个体数量
        l_f=length(front(i).fr);
        if l_f > 2
            %根据每一个目标进行排序，记录排序结果
            %sortedf的尺寸为[l_f x M]，第j列的值是根据第j个目标排列同rank个体的结果
            %sorted_indf记录排序前的index，可以称之为个体的编号
            for j=1:M
                [val,ind] = sortrows(chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+j));
                fmin = val(1);
                fmax = val(end);
                chromosome_sorted(ind(1)+rowsindex-1,V+M+1+j)=inf;
                chromosome_sorted(ind(end)+rowsindex-1,V+M+1+j)=inf;
                for k = 2:l_f-1
                    if fmax == fmin
                        chromosome_sorted(ind(k)+rowsindex-1,V+M+1+j)=inf;
                    else
                        chromosome_sorted(ind(k)+rowsindex-1,V+M+1+j)=(val(k+1)-val(k-1))/(fmax-fmin);
                    end

                end
            end

        else
            chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+M+2:V+M+1+M)=inf;
        end
        rowsindex = rowsindex + l_f;
    end
    %添加一列，计算每个个体在每个目标下的拥挤距离之和
    chromosome_sorted(:,V+M+1+M+1) = sum(chromosome_sorted(:,V+M+2:V+M+1+M),2);
    chromosome_NDS_CD1 = [chromosome_sorted(:,1:V+M) zeros(pop_size1,1) chromosome_sorted(:,V+M+1) chromosome_sorted(:,V+M+1+M+1)];
end

if problem_type==1 || problem_type==0.5
infpop=sortrows(infchromosome,V+M+1);
infpop=[infpop(:,1:V+M+1) (rank:rank-1+size(infpop,1))' inf*(ones(size(infpop,1),1))];
    for kk = (size(front,2)):(size(front,2))+(length(infchromosome))-1
        front(kk).fr= pop_size1+1;
    end
end
chromosome_NDS_CD = [chromosome_NDS_CD1;infpop]; 
end

function dom = judge_dominate(y1,y2)
tmp_res = y1 - y2;
if tmp_res <= 0 
    if tmp_res~= 0 
        dom = 1;
    else
        dom = 0;
    end
elseif tmp_res >= 0
    dom = -1;
else
    dom = 0;
end

end
