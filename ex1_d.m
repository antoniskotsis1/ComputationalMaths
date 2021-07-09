%-------------------------Exercise 1(D)-----------------------------------------
%-------------Euler's and Modified Euler's Method-------------------------------
%-------------------------------------------------------------------------------

clear; clc; close('all');
function createPlot(x,y,axisName,mehtodName,color)
    figure('name',mehtodName);
    plot(x,y,color);
    xlabel ("t axis");
    ylabel ([axisName "(t) axis"]);
    title ([axisName " Transform Plot"]);
    legend(mehtodName);
    
endfunction

function printSomeValues(z_Euler,Z_ModEuler)
  fprintf('\ti    Euler     ModifiedEuler\n')
  for i=1:100:length(z_Euler)-1
    fprintf('%10d%+10.4f  %+10.4f\n', i-1, z_Euler(i), Z_ModEuler(i));
  end
  
endfunction


function return_value = f1 (z,k) 
  persistent am = 3106; 
  persistent M=1;
  persistent g=9.81;
  persistent Kpz = 5;
  persistent Kdz = 15 - (am/1000);
  persistent zdes = am/200
  persistent Cz = 3+(am/5000);
  fz = M*g + Kpz*(zdes - z) - Kdz*k;
  
  return_value = (fz-g*M-Cz*abs(k)*k)/M;
endfunction


function return_value = f2 (y,k) 
  
  persistent am = 3106
  persistent Kpy = 5
  persistent Kdy = 20
  persistent ydes = am/3000
  persistent M = 1
  persistent g=9.81
  persistent Iz = 0.08
  persistent Cy = 5 

  Tz = Kpy*(ydes - y) - Kdy*k;
  
  return_value = (Tz - Cy*abs(k)*k)/Iz;
endfunction

function [t,z] = Euler(x0,k0,t_start,h,t_end,fun)

  t = (t_start:h:t_end); 
 
  z = zeros(size(t));
  k = zeros(size(t));
  z(1) = x0;
  k(1) = k0;

  for i = 1:1:length(t) - 1
      z(i + 1) = z(i) + h * k(i);
      k(i+1) = k(i) + h*fun(z(i),k(i));
      
  end
endfunction 
 
function [tt,zz] = ModifiedEuler(x0,k0,t_start,h,t_end,fun)
    tt = (t_start:h:t_end);

    kk = zeros(size(tt));
    zz = zeros(size(tt));
    zz(1) = x0;
    kk(1) = k0;
    for j = 1:1:length(tt) - 1
      zz(j + 1) = zz(j) + h * kk(j);
      kk(j+1) = kk(j) + h*fun(zz(j)+h/2,kk(j)+(h/2)*fun(zz(j),kk(j)));

    end
endfunction
%------------------Constant Values,Initial Values,step and time-----------------
am = 3106;
z0 = am/1000;
y0 = 0;
k0 = 0;
t_start = 0;
t_end = 30;
h = 0.001; 

%--------------------------Methods for the equation (1)------------------------- 
[t,z] = Euler(z0,k0,t_start,h,t_end,@f1);

[t,zz] =ModifiedEuler(z0,k0,t_start,h,t_end,@f1);

%--------------------------Methods for the equation (2)-------------------------
[t,y] = Euler(y0,k0,t_start,h,t_end,@f2);

[t,yy] =ModifiedEuler(y0,k0,t_start,h,t_end,@f2);

printSomeValues(z,zz)
%------------------------------creating the plots-------------------------------
createPlot(t,z,'Z',"Euler's Method",'r')
%print -djpeg Z_Euler.jpg
createPlot(t,zz,'Z',"Modified Euler's Method",'b')
%print -djpeg Z_ModEuler.jpg
createPlot(t,y,"Y~(z's orientation)","Euler's Method",'g')
%print -djpeg Y_Euler.jpg
createPlot(t,yy,"Y~(z's orientation)","Modified Euler's Method",'k')
%print -djpeg Y_ModEuler.jpg