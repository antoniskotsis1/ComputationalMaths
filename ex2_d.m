%------------------------Exersice 2(D)------------------------------------------
%--calculates the actual solution and approaches the differintial equation 
%--with Euler's and Modified Euler's methods
%--with the initial values of exercise 1c
%----uncomment line 106 to create errors plot(may take a while)
%-------------------------------------------------------------------------------
clear; clc; close('all');

function [r1,r2] = findRoots(a,b,c)
    delta = (b^2)-(4*a*c);
    r1 = (-b+sqrt(delta))/2;
    r2 = (-b-sqrt(delta))/2;
endfunction

function printSomeValues(actual_z,z_Euler,Z_ModEuler)
  fprintf('\ti   Actual      Euler      ModifiedEuler\n')
  t=(0:0.001:30);
  for i=1:100:length(z_Euler)-1
    fprintf('%10d%+10.4f %+10.4f %+10.4f\n', i-1,actual_z(t)(i), z_Euler(i), Z_ModEuler(i));
  end
  
endfunction

function [t,z] = Euler(x0,k0,t_start,h,t_end,fun)
  
  t = (t_start:h:t_end); 
  z = zeros(size(t));
  z(1) = x0;
  k(1) = k0;

  for i = 1:1:length(t) - 1
      z(i + 1) = z(i) + h * k(i);
      k(i+1) = k(i) + h*fun(z(i),k(i));
      %fprintf('%10d%+10.4f%+10.2f%+15.2f\n', i, z(i + 1), t(i + 1), fun(z(i + 1),k(i+1)));
  end
endfunction 


function findMethodsError(actual,euler,modEuler,t)
  fprintf('\t i   eulerError     modEulerError    actual\n');
  for i=1:1:length(t)-1;
    
    euler_error(i+1) = actual(t)(i+1)-euler(i+1);
    mod_euler_error(i+1) = actual(t)(i+1)-modEuler(i+1);
    fprintf('%10d%+10.4f%+10.4f%+15.4f\n', i-1, euler_error(i+1), mod_euler_error(i+1),actual(t)(i));
  end
  figure('name','Error Plot');
  plot(t,euler_error,t,mod_euler_error)
  xlabel ('time')
  ylabel('error')
  legend('Euler','ModEuler')
  print -djpeg errors_plot_overtime.jpg
  title("Error's plot over time")
endfunction
 
function [tt,zz] = ModifiedEuler(x0,k0,t_start,h,t_end,fun,input)
    tt = (t_start:h:t_end);

    zz = zeros(size(tt));
    kk = zeros(size(tt));
    zz(1) = x0;
    kk(1) = k0;
    for j = 1:1:length(tt) - 1
      zz(j + 1) = zz(j) + h * kk(j);
      kk(j+1) = kk(j) + h*fun(zz(j)+h/2,kk(j)+(h/2)*fun(zz(j),kk(j)));
     % fprintf('%10d%+10.4f%+10.2f%+15.2f\n', j, zz(j+1), tt(j + 1), fun(zz(j)+h/2,kk(j)+(h/2)*fun(zz(j),kk(j))));
    end

endfunction
%----------initial values-------------------------------
am=3106;
Kpz=5;
Kdz=15+(am/1000);
z0 =am/1000;
k0 = 0;
Zdes = am/200;
Cz = 3+(am/5000);
M=1;
g=9.81;
t_start = 0;
t_end = 30;
h=0.001;
%------------differencial equation------------------------
f = @(z,k) -(Kdz+Cz)*k - Kpz*z+Kpz*Zdes;

%------------solving characteristic equation--------------
[r1,r2] = findRoots(1,(Kdz+Cz),Kpz);

%------solving 2x2 system to find the general solution-----
A = [1 1;r1 r2];
b = [z0-Zdes;0];
x = inv(A)*b;
c1=x(1);
c2=x(2);
%----------------general solution--------------------------
general_solution = @(t) 15.53 + c1*(e.^(r1*t)) + c2*(e.^(r2*t));

t = 0:0.001:30 ;
[euler_t,euler_z] = Euler(z0,k0,t_start,h,t_end,f);
[modEuler_t,modEuler_z] = ModifiedEuler(z0,k0,t_start,h,t_end,f);

printSomeValues(general_solution,euler_z,modEuler_z)

%-----uncomment to generate error's plot (this may take a while)----------------

%findMethodsError(general_solution,euler_z,modEuler_z,t)



%--------------------plotting----------------------------
figure('name','Solution Plot');
plot(t,general_solution(t),modEuler_t,modEuler_z,euler_t,euler_z)
xlabel('time')
ylabel('z(t)')
title("Z's Transform Plot over time")
legend('Analytical Solytion','Modified Euler','Euler')
%print -djpeg final_plot.jpeg




