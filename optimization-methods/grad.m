function [k,x,val] = grad(fun,gfun,x0,epsilon)
% 功能：最速下降法求解无约束优化问题：min f(x)
%输入：fun,gfun分别是目标函数及其梯度，x0是初始点，epsilon为容许误差
%输出：k是迭代次数，x，val分别是近似最优点和最优值


maxk=5000; %最大迭代次数
beta=0.5; sigma=0.4;
k=0;
while(k<maxk)
    gk=feval(gfun,x0); %计算梯度
    dk=-gk; %计算搜索方向
    if(norm(gk)<epsilon), break; end %检验终止原则
    m=0; mk=0;
    while(m<20) %Armijo搜索求步长
        if(feval(fun,x0+beta^m*dk)<=feval(fun,x0)+sigma*beta^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    x0=x0+beta^mk*dk;
    k=k+1;
end
x=x0;
val=feval(fun,x0);
    