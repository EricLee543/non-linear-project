function child_offspring  = crossover(parent_selected,lb,ub,V)
[N] = size(parent_selected,1);
etac = 20;    % distribution index for crossover
lb1=lb';
ub1=ub';
rc=randperm(N);
for i=1:(N/2)
    parent1=parent_selected((rc(2*i-1)),:);
    parent2=parent_selected((rc(2*i)),:);
    if (isequal(parent1,parent2))==1 && rand(1)>0.5
        child1=parent1;
        child2=parent2;
    else 
        for j = 1: V  
            if parent1(j)<parent2(j)
               beta(j)= 1 + (1/(parent2(j)-parent1(j)))*(min((parent1(j)-lb1(j)),(ub1(j)-parent2(j))));
            else
               beta(j)= 1 + (1/(parent1(j)-parent2(j)))*(min((parent2(j)-lb1(j)),(ub1(j)-parent1(j))));
            end   
        end
        u=rand(1,V);
        alpha=2-beta.^-(etac+1);
        beta_q=(u<=(1./alpha)).*(u.*alpha).^(1/(etac+1))+(u>(1./alpha)).*(1./(2 - u.*alpha)).^(1/(etac+1));
        child1=0.5*(((1 + beta_q).*parent1) + (1 - beta_q).*parent2);
        child2=0.5*(((1 - beta_q).*parent1) + (1 + beta_q).*parent2);
    end
    child_offspring((rc(2*i-1)),:)=child1;
    child_offspring((rc(2*i)),:)=child2;
end
end