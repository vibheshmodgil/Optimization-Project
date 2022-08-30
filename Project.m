
[a,b] = Bounding_Phase_method() ;
secant_method(a,b) ;

function z = secant_method(x1 , x2)
e = 0.1 ;
z = 1;
k = 0 ;
graph = [] ;
derivative = 1 ;
while abs(derivative) > e
     k = k + 1 ;
     f_dash_x2 = D(x2) ;
     f_dash_x1 = D(x1) ;
    z = x2 - (f_dash_x2/(f_dash_x2 - f_dash_x1)*(x2 - x1)) ;
    f_dash_z = D(z) ;
    graph(k) = Objective_Fun(z) ;
    if f_dash_z < 0
        x1 = z ;
        derivative =  f_dash_z ;
    elseif f_dash_z > 0
        x2 = z ;
         derivative =  f_dash_z ;
    end
end
 plot(1:k , graph) ;
 title("Convergence plot for Secant Method") ;
 xlabel("No. of iterations ") ;
 ylabel("Function Value ") ;
 fprintf("The no. of iterations done for secant method is : %d\n",k)
 fprintf("The solution to the given problem is : %f \n",z) ;
end
function [a,b] = Bounding_Phase_method()
x0 = input('Enter the initial approximation : '); % initial approximation
d = input('Enter the delta value : ') ; %delta value
x1 = x0 + d ;
y0 = Objective_Fun(x0) ;
y1 = Objective_Fun(x1) ;
y_1 = Objective_Fun(x0 - d) ;

if ((y1 > y0)&&...
        (y_1 < y0))
    d = -d ;
    y1 = y_1 ;
    x1 = x0 + d ;
elseif ((y1 < y0)&&...
        (y_1 > y0))
     d = d ;
else 
   error("Something wrong with input values \n%s ",...
       "Try different value of x0 and delta" ) ;
end
% we have decided in which direction to go till now 
x = 0 ;
n = 1 ;
fprintf("The range for iteration %d is : (%f , %f )\n",n,x0,x1) ;

while y0 > y1
 n = n + 1 ;
 x = x0 ;
 x0 = x1 ;
 x1 = x1 + (2^n)*d ;
 fprintf("The range for iteration %d is : (%f , %f )\n",n,x,x1) ;
 y0 = y1 ;
 y1 = Objective_Fun(x1) ;
 a = x ; b = x1 ;
end
fprintf("The no. of iterations done for Bounding Phase Method method is : %d\n",n) ;
fprintf("The of function evaluations are : %d\n", n+2 );
fprintf("The final solution is : (%f , %f )\n",x,x1) ;
end
    
function fun_val = Objective_Fun(x)
     fun_val = -1*((2*x-5)^4 - (x^2-1)^3)    ;           
    %fun_val = -1*(20*sin(x) - 15*x^2) ;
    %fun_val = -1*(4*x*sin(x)) ;
    % fun_val = 2*(x-3)^2 + exp(0.5*x^2) ;
end

 function d = D(a)
  b = 0.001 ;
  d = (Objective_Fun(a + b) - Objective_Fun(a))/b ;
 end
 