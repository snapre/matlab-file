clc; clear
n = 256; m = 128;
A = randn(m,n);
u = sprandn(n,1,0.1);
b = A*u;
figure(1);subplot(2,1,1);plot(1:n,u);title('exact solu');

cvx_solver sdpt3
cvx_begin
    variable x(n)
    minimize (norm(x,1))
    subject to
        A*x == b
cvx_end
subplot(2,1,2);plot(1:n, x);title('l1 solu');
        
        
    
