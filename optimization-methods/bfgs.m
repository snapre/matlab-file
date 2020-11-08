function [k,x,val]=bfgs(fun,gfun,x0,varargin)
%功能：BFGS算法求解无约束问题：min f(x)
%输入：fun，gfun分别是目标函数及其梯度，x0是初始点，varargin是输入的
%   可变参数变量,简单调用BFGS时可以忽略它,但若其他程序循环调用该程
%   序时将发挥重要的作用
%输出：k是迭代次数，x,val分别是近似最优点和最优值
N=1000; %给出最大迭代次数
epsilon=1.e-5; %给定容许误差

beta=0.55; sigma=0.4;
n=length(x0); Bk=eye(n);
k=0;
while(k<N)
    gk=feval(gfun,x0,varargin{:}); %计算梯度
    if(norm(gk)<epsilon), break; end %检验终止准则
    dk=-Bk\gk; %解方程组Bk*dk=-gk,计算搜索方向
    m=0; mk=0;
    while(m<20) %用Armijo搜索求步长
        newf=feval(fun,x0+beta^m*dk ,varargin{:});
        oldf=feval(fun,x0,varargin{:});
        if(newf<=oldf+sigma*beta^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    %BFGS校正
    x=x0+beta^mk*dk;
    sk=x-x0;
    yk=feval(gfun,x,varargin{:})-gk;
    if(yk'*sk>0)
        Bk=Bk-(Bk*sk*sk'*Bk)/(sk'*Bk*sk)+(yk*yk')/(yk'*sk);
    end
    k=k+1;
    x0=x;
end
val=feval(fun,x0,varargin{:});