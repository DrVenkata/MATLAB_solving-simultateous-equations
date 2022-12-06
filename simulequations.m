clear all 
clc
% a is coefficient matrix
a=[10 -2 -1 -1;-1 10 -1 -1;-1 -1 10 -2;-1 -1 2 10]
% c is constant matrix
c=[3;15;27;-9]

% To find unknown variables  x1,x2,x3,x4
%Matrix inversion method
fprintf('Matrix inversion method\n')
Xiv=matrixinversion(a,c)

% Gauss jordan Elimination Method
fprintf('Gauss jordan Elimination Method\n')
Xgje=GaussJordanElimination(a,c)

%Gauss Seidal Iterative method
fprintf('Gauss Seidal Iterative method\n')
[Xgs,i]=GaussSeidalIterative(a,c)
fprintf('solution converged at iteration %d \n',i)

%Jacobi Iterative method\n
fprintf('Jacobi Iterative method\n')
[Xji,j]=JacobiIterative(a,c)
fprintf('solution converged at iteration %d \n',j)



function Xgje=GaussJordanElimination(a,c)
    % Augmented matrix 
    % i.e. coefficient matrix(a) and constant matrix(c)
            M=[a c];
    % Applying Gauss Jordan Elimination method of obtaining identity % matrix and variable matrix side by side i.e. R matrix 
    R=rref(M);
    % x1,x2,x3,x4 from Gauss Jordan Elimination method
    Xgje=R(:,end);
end


function Xiv=matrixinversion(a,c)
A=inv(a);
Xiv=A*c;
end

function [Xgs,i]=GaussSeidalIterative(a,c)
    % n is no of variables
    N=size(a);
    n=N(1);
    %Extract only the elements below the main diagonal of 'a' matrix
    L=tril(a,-1) ;
    %Extract only the elements above the main diagonal of 'a' matrix
    U=triu(a,1) ;
    %Extract only diagnol elements of 'a' matrix
    D=diag(diag(a));
    X=zeros(n,1);
    i=1;
    d=0;
    while i>=1
        % Gauss seidal iteration equation
        X=(inv(D+L))*(c-(U*X));
        % store each iteration value in Xi array
        Xi(:,i)=X; 
        % Convergence Criteria satisfies break the loop
        
            if a*X-c <= 0.00001
               d=1 ;
            end
       
        if d~=1 
           i=i+1;
        else            
            break;
        end
    end 
    Xi
    Xgs=Xi(:,i) ;
end

function [Xji,j]=JacobiIterative(a,c)
    % n is no of variables
    N=size(a);
    n=N(1);
    %Extract only the elements below the main diagonal of 'a' matrix
    L=tril(a,-1) ;
    %Extract only the elements above the main diagonal of 'a' matrix
    U=triu(a,1) ;
    %Extract only diagonal elements of 'a' matrix
    D=diag(diag(a));
    X=zeros(n,1); 
    syms d
    j=1;
    while j>=1
        % Jacobi iteration equation
        X=(inv(D))*(c-((L+U)*X));
        % store each iteration value in Xi array
        Xii(:,j)=X;
            % Convergence Criteria difference between two iteration
              
            if abs(a*X-c) <= 0.00001 
               d=1 ;
            end
        
        if d~=1 
           j=j+1;
        else            
            break;
        end
    end
    Xii
    Xji=Xii(:,end) ;
end

