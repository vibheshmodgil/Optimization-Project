
n = input("Enter the no. of variables : ") ; % No. 0f variables
%-------code to get initial values of the variables-------------
x0 = zeros(n,1) ;
for i = 1:n 
    fprintf("Enter the value of x%d : ",i)
    x0(i) = input("") ;
end
fprintf("\n")
solution = project_phase_2(x0) 

%-----Generate Randon initial value and check bunch of them at once-----
% no_val = input(" No.of variable has to be found out : ") ;
% a = input("Enter the lower limit of range : "); % Lower value of the range 
% b = input("Enter the upper limit of range ");   % Higher value of the range
% for i = 1:10 % Checks the solution for 10 random initial values of variables
% x0 = (b-a).*rand(no_val,1) + a  % Generates random points between a and b 
% solution = project_phase_2(x0) 
% end

%----------------PROJECT PHASE 2 ------------------
function answer = project_phase_2(x0)
 % x0 = Initial Approximation 
 k = 1 ; % Iteration Counter 
 function_value = [] ; % For the iteration plot
 e = 0.001 ; % Accuracy reqired for activation of termination Condition
 m = 50 ; % Maxm no.  of iteration 
 project_1_function_eval = 0  ; % initilizing total no. of function evaluations in project phase 1 
 project_2_function_eval = 0 ; % initilizing total no. of function evaluations in project phase 2
 n_of_val = length(x0) ; % No. of variables

while true  
 gradient_0 = grad(x0) ;  % Gradient at initial Point
 project_2_function_eval =  project_2_function_eval +  n_of_val*2 ;
% N0. 0f variables = No. of first order partial differntials in gradient ;
% Each one requres 2 function evaluations
if magnitude(gradient_0) < e % Termination condition 1 
    answer = x0 ;
    break
end
  function_value(k) = Function(x0) ; % For the iteration plot
hessian = hess(x0) ; % Hessian matrix itinitial Point
project_2_function_eval =  project_2_function_eval +  4*n_of_val^2-1 ; 

%------checking that Hessian matrix of positive definite or not --------
% eig_hessian = eig(hessian) ; % Eigen Values of Hessian Matrix 
% for i = 1 : length(x0) 
%  if eig_hessian(i) < 0 
%     error("For given initial conditions hessian matrix is not Positive definate") ;
%     break
%  end
% end
direction = inv(hessian)*gradient_0 ;% Direction for unidirectional Search 
[alfa , function_eval] = project_1(x0,direction) ; % Value of alfa by unidirectional Search using Project Phase 1
project_1_function_eval = project_1_function_eval + function_eval ; % Total no. of function evaluations in project phase 1 code
x1 = new_point(x0,alfa,direction) ;% New point after substituting the value
gradient_1 = grad(x1) ; % Gradient at new point
 project_2_function_eval =  project_2_function_eval +  n_of_val*2 ; % No. of function evaluations for this gradient
if (gradient_0*gradient_1' < e ) % Termination Condition 2 
    answer = x1 ;
    break
end
if ((magnitude(x1 - x0))/(magnitude(x0))) < e % Termination Condition 3
answer = x1 ;
   break
end
k = k + 1  ;
x0 = x1 ;
end
function_value(k+1) = Function(x1) ;
fprintf("The No. of iterations are : %d \n",k) 
fprintf("The No. of function evaluations by project phase 1 code  : %d \n",project_1_function_eval)
fprintf("The No. of function evaluations by project phase 2 code  : %d \n",project_2_function_eval)
fprintf("The minimum function value is : %f \n",Function(x1))
iteration = [0:k]  ;
plot(iteration,function_value)
title('Iteration plot') ;
xlabel(' No. of iteration ') ;
ylabel(' Function value ') ;
end
%--------------Ploting the iteration Plot -----------

% ------- Multivariable Function ------------
function fun_val = Function(x)

%--------------QUESTION 1 SUM SQUARES FUNCTION------------------
% n_val = length(x) ;
% fun_val = 0 ;
% for i = 1:n_val
%   fun_val = fun_val + i*x(i)^2 ;
% end
% -------------QUESTION 2 ROSENBROCK FUNCTION --------------------- 
% n_val = length(x) ;
% fun_val = 0 ;
% for i = 1:n_val-1
%   fun_val = fun_val + 100*(x(i+1) - x(i)^2)^2 + (x(i) - 1)^2 ;
% end

% %------------QUESTION 3 DIXON PRICE FUNCTION ------------------------
% n_val = length(x) ;
% fun_val = (x(1) - 1)^2  ;
% for i = 2:n_val
%   fun_val = fun_val + i*(2*x(i)^2 - x(i-1))^2 ;
%  end
%-------------QUESTION 4 TRID FUNCTION -------------------------
% n_val = length(x) ;
% fun_val = 0  ;
% for i = 1:n_val
%   fun_val = fun_val + (x(i)-1)^2 ;
% end
% for i = 2:n_val
%   fun_val = fun_val - x(i)*x(i-1) ;
% end
%-----------QUESTION 5 ZAKHAROV FUNCTION ------------------------------
n_val = length(x) ;
fun_val = 0  ;
first_term = 0 ;
sum = 0 ;
for i = 1:n_val
  first_term = first_term + x(i)^2 ;
  sum = sum + 0.5*i*x(i) ;
end
fun_val = first_term + sum^2 + sum^4 ;

%--------------Himmelblau function----------------------------------
 %   fun_val = ((x(1))^2 + x(2) - 11)^2 + (x(1) + (x(2))^2 - 7)^2 ;
end

%-----Function for gradient ----------------
function gradient = grad(x) 
gradient = zeros(length(x),1) ;
h = 0.001 ;
for i = 1:length(x)
    y = x ;
 y(i) = y(i)+h ;
 a = Function(y) ;
 y(i) = y(i)-2*h ;
 b = Function(y) ;
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
function hessian = hess(x)
l = length(x) ;
hessian = zeros(l,l) ;
h = 0.001 ;
for i = 1:l
for j = 1:l
if i == j 
 y = x ;
 y(i) = y(i)+h ;
 a = Function(y) ;
 y(i) = y(i)-2*h ;
 b = Function(y) ;
 c = Function(x) ;
 hessian(i,j) =  (a+b-2*c)/(h^2) ;
else
 a = x ;
 b = x ;
 c = x ;
 d = x ;
 a(i) = a(i) + h ;
 a(j) = a(j) + h ;
 first_term = Function(a) ;
 b(i) = b(i)+ h ;
 b(j) = b(j)- h ;
 second_term = Function(b) ;
 c(i) = c(i) - h ;
 c(j) = c(j) + h ;
 third_term = Function(c) ;
 d(i) = d(i) - h ;
 d(j) = d(j) - h ;
 forth_term = Function(d) ;
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
function fun_val = Objective_Fun(y,x0,direction)
l = length(x0);
for i = 1:l
x0(i) = x0(i) - y*direction(i) ;
end
fun_val = Function(x0) ;
end


%------------Function of single variable optimization----------
%------------Same as Project Phase 1 --------------------------
function [sol , function_eval] = project_1(x0 , direction)
[a,b,function_eval] = Bounding_Phase_method(-2 , 0.1 ,x0 ,direction) ;
[ sol , function_eval] = secant_method(a,b ,x0 ,direction,function_eval) ;


    function [ z , function_eval] = secant_method(y1 , y2 ,x0 ,direction,function_eval)
e = 0.01 ; % error allowed in derivative 
k = 0 ; % No. of iterations
% graph = [] ; % Enpty set to store value of function at each iteration
derivative = 1 ; % To initially initialize the while loop
while abs(derivative) > e
     k = k + 1 ;
     f_dash_x2 = D(y2 , x0 ,direction) ;
     f_dash_x1 = D(y1 ,x0 ,direction) ;
    z = y2 - (f_dash_x2/(f_dash_x2 - f_dash_x1)*(y2 - y1)) ; % Secant formula for next approximation
    f_dash_z = D(z ,x0 ,direction) ;
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
    function [a,b,function_eval] = Bounding_Phase_method(y0 , d ,x0 ,direction)
% y0 = initial approximation
% d = delta value
y1 = y0 + d ;
f0 = Objective_Fun(y0 , x0 ,direction) ;
f1 = Objective_Fun(y1 , x0 ,direction) ;
f_1 = Objective_Fun((y0 - d), x0 , direction ) ;
if ((f1 > f0)&&(f_1 < f0)) % f(x0 + d) > f(x0) > f(xo - d)
    d = -d ;
    f1 = f_1 ;
    y1 = y0 + d ;
elseif ((f1 < f0)&&(f_1 > f0)) % f(x0 + d) < f(x0) < f(x0 - d)
     d = d ;
else   % f(x0 + d ) > f(x0) < f(x0 - d) solution is in this range only
   error("Something wrong with input values \n%s ",...
       "Try different value of x0 and delta" ) ;
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
 f1 = Objective_Fun(y1 ,x0 ,direction) ;
 a = y ; b = y1 ;
 function_eval = n+2  ;
end
% fprintf("The no. of iterations done for Bounding Phase Method method is : %d\n",n) ;
% fprintf("The of function evaluations are : %d\n", n+2 );
% fprintf("The final solution is : (%f , %f )\n",y,x1) ;
end

 function d = D(y0,x0 ,direction)
  h = 0.001 ;
  d = (Objective_Fun((y0 + h),x0 ,direction ) - Objective_Fun((y0 - h),x0 ,direction))/(2*h) ; % Forward Difference Method for derivative
 end
end
   
    