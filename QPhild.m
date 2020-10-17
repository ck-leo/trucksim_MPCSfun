function x = QPhild(E,F,M,gama)
% Hildreth's algorithm
% solve the problem:
% min J = 1/2*(eta'*H*eta)+eta'*f;  (J = 1/2*x'*E*x+x'*F)
% A_cons*eta<=b                     (M*x<=gama)
% E=H
% F=f
% M=A_cons
% gamma=b
% eta=x
[n1,m1]=size(M);
x = -E\F;
kk=0;
for i=1:n1
    if(M(i,:)*x>gama(i))
        kk=kk+1;
    else
        kk = kk+0;
    end
end
if(kk==0) return; end
H = M*(E\M');
K = M*(E\F)+gama;
[n,m] = size(K);
x_ini = zeros(n,m);
lambda = x_ini;
a1 = 10;
for km = 1:38
% find the elements in the solution vector one by one
% km could be larger if the Lagranger multiplier has a slow
% convergence rate
    lambda_p = lambda;
    for i = 1:n
        w = H(i,:)*lambda - H(i,i)*lambda(i,1);
        w = w+K(i,1);
        la = -w/H(i,i);
        lambda(i,1) = max(0,la);
    end
    a1 = (lambda-lambda_p)'*(lambda-lambda_p);
    if(a1<10e-8) break; end
end
x = -E\F-E\M'*lambda;
end

