function [k,x,val]=dfp(fun,gfun,x0,epsilon,N)
%功能：DFP算法求解无约束问题：min f(x)
%输入：fun，gfun分别是目标函数及其梯度，x0是初始点，
%   epsilon是容许误差，N是最大迭代次数
%输出：k是迭代次数，x,val分别是近似最优点和最优值
if nargin<5, N=1000; end
if nargin<4, epsilon=1.e-5; end

beta=0.55; sigma=0.4;
n=length(x0); Hk=eye(n); k=0;
while(k<N)
    gk=feval(gfun,x0); %计算梯度
    if(norm(gk)<epsilon), break; end %检验终止准则
    dk=-Hk*gk; %计算搜索方向
    m=0; mk=0;
    while(m<20) %用Armijo搜索求步长
        if(feval(fun,x0+beta^m*dk)<=feval(fun,x0)+sigma*beta^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    %DFP校正
    x=x0+beta^mk*dk;
    sk=x-x0; yk=feval(gfun,x)-gk;
    if(sk'*yk>0)
        Hk=Hk-(Hk*yk*yk'*Hk)/(yk'*Hk*yk)+(sk*sk')/(sk'*yk);
    end
    k=k+1;
    x0=x;
end
val=feval(fun,x0);