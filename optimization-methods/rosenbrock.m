% 利用最速下降法求解Rosenbrock函数的局部极小值
% Meringue
% 2017/4/14
% ---------------------------
% ---------------------------
clc
clear all
close all

global xk dk iter

err = 1e-3; % 精度要求
iter = 0; % 迭代次数
iterMax = 1e5; % 最大迭代次数
x0 = [0,0]';% 初始点
[~,gk] = Rosenbrock(x0); % 求初始梯度


xk(:,1) = x0; 
while norm(gk) > err % 未达到精度要求
    if iter >= iterMax
        fprintf('达到最大迭代次数！');
        break
    end
    iter = iter+1;
    [~,gk] = Rosenbrock(xk(:,iter));
    dk = -gk; % 确定当前搜索方向

    % 确定搜索区间（2选1）
    % [lambdaMin,lambdaMax,~] = Trial(@Rosenbrock,xk(:,iter),dk,1,1,2,1e-3); % 进退搜索法
    lambdaMin = 0; lambdaMax = 10;% 步长一般很小，直接给出粗略区间

    % 精确直线搜索确定最优步长（2选1）
    % lambdak = P618(@Rosenbrock,xk(:,iter),dk,[lambdaMin,lambdaMax],1e-3); % 黄金分割法
    lambdak = fminbnd(@Rosenbrock_t,lambdaMin,lambdaMax); % 调用fminbnd函数
    xk(:,iter+1) = xk(:,iter)+lambdak*dk; % 搜索下一点
end

% 结果显示
fprintf('一共迭代了%d次',iter);
plot(xk(1,:),xk(2,:))
title('最速下降法求Rosenbrock函数的局部极小值')
xlabel('x_1'); ylabel('x_2');