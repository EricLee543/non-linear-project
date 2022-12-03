function nsga2_main(inputArg1,inputArg2)
pop_size = 20; % 种群数量
V = 3; % （x1,x2...xv) 的维度，若翼形为16维描述符V=16
M = 2; % cons的数量
no_runs = 1; %
lb = [0.5 -2 -2]; %初始化的上下界
ub = [1 2 2];
x = lb+((ub-lb).*rand(pop_size, V)); %初始种群
for i =1:pop_size
[ff(i,:) err(i,:)] =feval("f_t", x(i,:));
end
error_norm=normalisation(err);
population_init=[x ff error_norm];
[population front]=NDS_CD_cons(population_init,V,M); 


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