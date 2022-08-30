
% x1 = [5:0.1:16] ;
% x2 = [-19:0.1:4] ;
% [X1 , X2] = meshgrid(x1 ,x2) ;
% z = (X1-10).^3 + (X2-20).^3 ;
% contourf(X1,X2,z) ;
% colorbar ;
% hold on ;
% 
%  f = @(X1,X2) (X1-5).^2 + (X2-5).^2 -100 ;
%  fimplicit(f,[5 16 -19 4],'k','LineWidth',2)
% 
%   f = @(X1,X2) (X1-6).^2 + (X2-5).^2 -82.81 ;
%  fimplicit(f,[5 16 -19 4],'b','LineWidth',2)

 % ---------------------------------------Question 2--------------------
x1 = [0.9:0.01:2] ;
x2 = [3:0.01:5] ;
[X1 , X2] = meshgrid(x1 ,x2) ;
z = (((sin(2*pi.*X1)).^3).*sin(2*pi.*X2))./((X1.^3).*(X1+X2)) ;
contourf(X1,X2,z) ;
colorbar ;
hold on ;

  f = @(X1,X2) (X1).^2 - X2 +1 ;
  fimplicit(f,[1 2 3 5],'k','LineWidth',2)

  f = @(X1,X2) 1 + (X2-4).^2 - X1 ;
 fimplicit(f,[1 2 3 5],'b','LineWidth',2)

 xlabel('x1') ;
 ylabel('x2') ;
 title('contour plot') ;
  legend('((sin(2*pi*X1)).^3*sin(2*pi.*X2))/((X1^3)*(X1+X2))','X1^2 - X2 + 1 <= 0','1 - X1 + (X2-4)^2 <= 0')
r = 0.1 ;
n = input("Enter the no. of variables : ") ; % No. 0f variables
%-------code to get initial values of the variables-------------
x0 = zeros(n,1) ;
for i = 1:n 
    fprintf("Enter the value of x%d : ",i) ;
    x0(i) = input("") ;
end
fprintf("\n") 
fprintf("For iteration 1 : \n")
penalty = Function(x0,r,0)/r 
Function(x0,0,1)
fprintf("For iteration 2 : \n")
solution = project_phase_2(x0,r)      
penalty = Function(solution,r,0)/r 
Function(solution,0,1)
  plot([x0(1) solution(1)], [x0(2) solution(2)] ,'r--o','LineWidth',2)
k = 2 ; % iteration counter 

while true
k = k + 1 ; 
r = r*10 ;
fprintf("For iteration %d : \n",k)
solution_0  = solution ;
solution_1 = project_phase_2(solution_0,r) 
 plot([solution_0(1) solution_1(1)], [solution_0(2) solution_1(2)] ,'r--o','LineWidth',2)
if abs(Function(solution_0,r,1) - Function(solution_1,r/10,1)) < 0.001
    break
end
solution = solution_1 ;
penalty = Function(solution,r,0)/r 
Function(solution,0,1)
fprintf("\n")
end
Final_solution = solution_1
Function(solution_1,0,1)
% xlabel('x1') ;
%  ylabel('x2') ;
%  title('contour plot') ;
%  legend('(x1-10)^3 + (x2-20)^3','(x1-5)^2 + (x2-5)^2 - 100 >= 0','(x1-6)^2 + (x2-5)^2 - 82.81 <= 0')

 xlabel('x1') ;
 ylabel('x2') ;
 title('contour plot') ;
  legend('((sin(2*pi*X1)).^3*sin(2*pi.*X2))/((X1^3)*(X1+X2))','X1^2 - X2 + 1 <= 0','1 - X1 + (X2-4)^2 <= 0')


%----------------PROJECT PHASE 2 ------------------
function answer = project_phase_2(x0,r)
 % x0 = Initial Approximation 
 k = 1 ; % Iteration Counter 
 function_value = [] ; % For the iteration plot
 e = 0.001 ; % Accuracy reqired for activation of termination Condition
 m = 50 ; % Maxm no.  of iteration 
 project_1_function_eval = 0  ; % initilizing total no. of function evaluations in project phase 1 
 project_2_function_eval = 0 ; % initilizing total no. of function evaluations in project phase 2
 n_of_val = length(x0) ; % No. of variables

while true  
 gradient_0 = grad(x0,r) ;   % Gradient at initial Point
 project_2_function_eval =  project_2_function_eval +  n_of_val*2 ;
% N0. 0f variables = No. of first order partial differntials in gradient ;
% Each one requres 2 function evaluations
if magnitude(gradient_0) < e % Termination condition 1 
    answer = x0 ;
    break
end
  function_value(k) = Function(x0,r,1) ; % For the iteration plot
hessian = hess(x0,r) ; % Hessian matrix itinitial Point
project_2_function_eval =  project_2_function_eval +  4*n_of_val^2-n_of_val ; 
% n^2 second order partial derivative involved
% n diagonal elements needs 3n function evalutions 
% rest n^2-n elements needs 4*(n^2-n) function evaluations

%------checking that Hessian matrix of positive definite or not --------
% if k== 1
% eig_hessian = eig(hessian) ; % Eigen Values of Hessian Matrix 
% for i = 1 : length(x0) 
%  if eig_hessian(i) < 0 
%     error("For given initial conditions hessian matrix is not Positive definate") ;
%     break
%  end
% end
% end

direction = pinv(hessian)*gradient_0 ; % Direction for unidirectional Search 
[alfa , function_eval] = project_1(x0,direction,r) ; % Value of alfa by unidirectional Search using Project Phase 1
project_1_function_eval = project_1_function_eval + function_eval ; % Total no. of function evaluations in project phase 1 code
x1 = new_point(x0,alfa,direction) ;% New point after substituting the value
gradient_1 = grad(x1,r) ; % Gradient at new point
 project_2_function_eval =  project_2_function_eval +  n_of_val*2 ; % No. of function evaluations for this gradient
% Termination Condition 2 
 if (gradient_0*gradient_1' < e ) 
    answer = x1 ;
    break
end
% Termination Condition 3
if ((magnitude(x1 - x0))/(magnitude(x0))) < e 
answer = x1 ;
   break
end
if k > m % Termination Condition 4 For maxm no. of iterations 
    answer = x1 ;
    fprintf('Maxm no. of iterations exceeds our expectations \n') ;
    fprintf('Upto now : \n') ;
    break 
end
k = k + 1  ; % First iteration completes 2nd iteration starts 
x0 = x1 ;
end
% function_value(k+1) = Function(x1,r,1) ; % Above 2 expression are not cosidered if any above termination condition is used
% fprintf("The No. of iterations are : %d \n",k) 
% fprintf("The No. of function evaluations by project phase 1 code  : %d \n",project_1_function_eval)
% fprintf("The No. of function evaluations by project phase 2 code  : %d \n",project_2_function_eval)
% fprintf("The total No. of function evaluations by the code are : %d \n",project_1_function_eval+project_2_function_eval)
% fprintf("The minimum function value is : %f \n",Function(x1,r,1))

%--------------Ploting the iteration Plot -----------
% iteration = [0:k]  ;
% plot(iteration,function_value)
% title('Iteration plot') ;
% xlabel(' No. of iteration ') ;
% ylabel(' Function value ') ;
 end

% ------- Multivariable Function ------------
 function fun_val = Function(x,r,a)
%--------------Question 1----------------------------------
% fun_val_1 = (x(1)-10)^3 + (x(2)-20)^3  ;
% fun_val_2 = min([0 , ((x(1)-5)^2 + (x(2)-5)^2 -100)]);
% fun_val_3 = min([0 , -1*((x(1)-6)^2 + (x(2)-5)^2 -82.81)]) ;
% fun_val_4 = min([0 , x(1)-13]) ;
% fun_val_5 = min([0 , 20-x(1)]) ;
% fun_val_6 = min([0 , x(2)]) ;
% fun_val_7 = min([0 , 4-x(2)]) ;
% penalty = fun_val_2^2 + fun_val_3^2 + fun_val_4^2 + fun_val_5^2 + fun_val_6^2 + fun_val_7^2 ;
% fun_val = a*fun_val_1 + r*penalty ;
end

%------------------Question 2------------------------------
% fun_val_1 = -1*((((sin(2*pi*x(1)))^3)*sin(2*pi*x(2)))/((x(1)^3)*(x(1)+x(2)))) ;
% fun_val_2 = min([0 , -1*((x(1))^2 - x(2) +1) ]) ;
% fun_val_3 = min([0 , -1*(1 + (x(2)-4)^2 - x(1))]) ;
% fun_val_4 = min([0 , x(1)]) ;
% fun_val_5 = min([0 , 10-x(1)]) ;
% fun_val_6 = min([0 , x(2)]) ;
% fun_val_7 = min([0 , 10-x(2)]) ;
% penalty = fun_val_2^2 + fun_val_3^2 + fun_val_4^2 + fun_val_5^2 + fun_val_6^2 + fun_val_7^2 ;
% fun_val = a*fun_val_1 + r*penalty ;
% end

%----------------------Question 3-------------------------------------- 
% fun_val_1 = x(1) + x(2) + x(3) ;
% fun_val_2 = min([0 , -1*(-1 + 0.0025*(x(4)+x(6))) ]) ;
% fun_val_3 = min([0 , -1*(-1 + 0.0025*(-x(4)+x(5)+x(7)))]) ;
% fun_val_4 = min([0 , -1*(-1 + 0.01*(-x(6)+x(8)))]) ;
% fun_val_5 = min([0 , -1*(100*x(1) - (x(1)*x(6)) + 833.33252*x(4) - 83333.333)]) ;
% fun_val_6 = min([0 , -1*(x(2)*(x(4) - x(7)) -1250*x(4) + 1250*x(5))]) ;
% fun_val_7 = min([0 , -1*(x(3)*(x(5)-x(8)) - 2500*(x(5)) + 1250000)]) ;
% penalty = fun_val_2^2 + fun_val_3^2 + fun_val_4^2 + fun_val_5^2 + fun_val_6^2 + fun_val_7^2 ;
% fun_val = a*fun_val_1 + r*penalty ;
% end


%--------------For Penalty-----------------
% function p = penalty(x)
% fun_val_2 = min([0 , ((x(1)-5)^2 + (x(2)-5)^2 -100)]);
% fun_val_3 = min([0 , -1*((x(1)-6)^2 + (x(2)-5)^2 -82.81)]) ;
% fun_val_4 = min([0 , x(1)-13]) ;
% fun_val_5 = min([0 , 20-x(1)]) ;
% fun_val_6 = min([0 , x(2)]) ;
% fun_val_7 = min([0 , 4-x(2)]) ;
% p = fun_val_2^2 + fun_val_3^2 + fun_val_4^2 + fun_val_5^2 + fun_val_6^2 + fun_val_7^2 ;
% end


%-----Function for gradient ----------------
function gradient = grad(x,r) 
gradient = zeros(length(x),1) ;
h = 0.001 ;
for i = 1:length(x)
    y = x ;
 y(i) = y(i)+h ;
 a = Function(y,r,1) ;
 y(i) = y(i)-2*h ;
 b = Function(y,r,1) ;
 gradient(i) = (a - b)/(2*h) ;
end
end

%-------------Function for magnitude of a vector--------------
function m = magnitude(gradient)
magnitude_square = gradient.*gradient ;
magnitude_square_sum = sum(magnitude_square) ;
m = sqrt(magnitude_square_sum) ;
end

%------------Function for Hessian matrix -------------------
function hessian = hess(x,r)
l = length(x) ;
hessian = zeros(l,l) ;
h = 0.001 ;
for i = 1:l
for j = 1:l
if i == j  % For 
 y = x ;
 y(i) = y(i)+h ;
 a = Function(y,r,1) ;
 y(i) = y(i)-2*h ;
 b = Function(y,r,1) ;
 c = Function(x,r,1) ;
 hessian(i,j) =  (a+b-2*c)/(h^2) ;
else
 a = x ;
 b = x ;
 c = x ;
 d = x ;
 a(i) = a(i) + h ;
 a(j) = a(j) + h ;
 first_term = Function(a,r,1) ;
 b(i) = b(i)+ h ;
 b(j) = b(j)- h ;
 second_term = Function(b,r,1) ;
 c(i) = c(i) - h ;
 c(j) = c(j) + h ;
 third_term = Function(c,r,1) ;
 d(i) = d(i) - h ;
 d(j) = d(j) - h ;
 forth_term = Function(d,r,1) ;
hessian(i,j)= (first_term - second_term - third_term + forth_term)/(4*h^2) ;
end
end
end
end

%------Function to find new vector after finding the value of alfa-------
function new_points = new_point(x,alfa,direction)
l = length(x) ;
for i = 1:l
    x(i) = x(i) - alfa*direction(i) ;
end
new_points = x ;
end

%---Fuction creation for single variable optimization for gradient---
function fun_val = Objective_Fun(y,x0,direction,r)
l = length(x0);
for i = 1:l
x0(i) = x0(i) - y*direction(i) ;
end
fun_val = Function(x0,r,1) ;
end


%------------Function of single variable optimization----------
%------------Same as Project Phase 1 --------------------------
function [sol , function_eval] = project_1(x0 , direction,r)
[a,b,function_eval] = Bounding_Phase_method(-1 , 0.1 ,x0 ,direction,r) ;
[ sol , function_eval] = secant_method(a,b ,x0 ,direction,function_eval,r) ;


    function [ z , function_eval] = secant_method(y1 , y2 ,x0 ,direction,function_eval,r)
e = 0.01 ; % error allowed in derivative 
k = 0 ; % No. of iterations
% graph = [] ; % Enpty set to store value of function at each iteration
derivative = 1 ; % To initially initialize the while loop
while abs(derivative) > e
     k = k + 1 ;
     f_dash_x2 = D(y2 , x0 ,direction,r) ;
     f_dash_x1 = D(y1 ,x0 ,direction,r) ;
    z = y2 - (f_dash_x2/(f_dash_x2 - f_dash_x1)*(y2 - y1)) ; % Secant formula for next approximation
    f_dash_z = D(z ,x0 ,direction,r) ;
% graph(k) = Objective_Fun(z) ; % Stores the value of function in each iteration
    derivative =  f_dash_z ; % To check the while codition
    if f_dash_z < 0   % Conditions to ensure minima lies in x1 and x2  
        y1 = z ;
    elseif f_dash_z > 0
        y2 = z ;
    end
end
function_eval = function_eval + 2*(k+1) ;
%  plot(1:k , graph) ; % For convergence plot
%  title("Convergence plot for Secant Method") ;
%  xlabel("No. of iterations ") ;
%  ylabel("Function Value ") ;
%  fprintf("The no. of iterations done for secant method is : %d\n",k)
%  fprintf("The solution to the given problem is : %f \n",z) ;
end
    function [a,b,function_eval] = Bounding_Phase_method(y0 , d ,x0 ,direction,r)
% y0 = initial approximation
% d = delta value
y1 = y0 + d ;
f0 = Objective_Fun(y0 , x0 ,direction,r) ;
f1 = Objective_Fun(y1 , x0 ,direction,r) ;
f_1 = Objective_Fun((y0 - d), x0 , direction,r ) ;
if ((f1 > f0)&&(f_1 < f0)) % f(x0 + d) > f(x0) > f(xo - d)
    d = -d ;
    f1 = f_1 ;
    y1 = y0 + d ;
elseif ((f1 < f0)&&(f_1 > f0)) % f(x0 + d) < f(x0) < f(x0 - d)
     d = d ;
else   % f(x0 + d ) > f(x0) < f(x0 - d) solution is in this range only
    a = y0 - d ;
    b = y0 + d ;
    function_eval = 3 ;
    return
%    error("Something wrong with input values \n%s ",...
%        "Try different value of x0 and delta" ) ;
end
% we have decided in which direction to go till now 
y = 0 ;
n = 1 ; % For no. if iterations .
% fprintf("The range for iteration %d is : (%f , %f )\n",n,f0,x1) ;

while f0 > f1 % Works till f(x0) > f(x0 + d)
 n = n + 1 ;
 y = y0 ;
 y0 = y1 ;
 y1 = y1 + (2^n)*d ; % new approximation formula by bounding phase method
%  fprintf("The range for iteration %d is : (%f , %f )\n",n,x,x1) ; % Prints the range after each iteration
 f0 = f1 ;
 f1 = Objective_Fun(y1 ,x0 ,direction,r) ;
 a = y ; b = y1 ;
 function_eval = n+2  ;
end
% fprintf("The no. of iterations done for Bounding Phase Method method is : %d\n",n) ;
% fprintf("The of function evaluations are : %d\n", n+2 );
% fprintf("The final solution is : (%f , %f )\n",y,x1) ;
end

 function d = D(y0,x0 ,direction,r)
  h = 0.001 ;
  d = (Objective_Fun((y0 + h),x0 ,direction,r ) - Objective_Fun((y0 - h),x0 ,direction,r))/(2*h) ; % Forward Difference Method for derivative
 end
end
   
    