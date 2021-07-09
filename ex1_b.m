%------------------------Exercise 1(B)------------------------------------------
%-------------Euler's and Modified Euler's Method-------------------------------
%-------------------------------------------------------------------------------
clear; clc; close('all');


function createPlot(x,y,axisName,mehtodName,color,input)
    figure('name',mehtodName);
    plot(x,y,color);
    xlabel ("t axis");
    ylabel ([axisName "(t) axis"]);
    title ([axisName " Transform Plot for input " input]);
    legend(mehtodName);
endfunction

function printSomeValues(z_Euler,Z_ModEuler)
  fprintf('\ti    Euler     ModifiedEuler\n')
  for i=1:100:length(z_Euler)-1
    fprintf('%10d%+10.4f  %+10.4f\n', i-1, z_Euler(i), Z_ModEuler(i));
  end
  
endfunction

function return_value = f1 (z,k,fz) 
  global am;
  global M ;
  global g;
  persistent Cz = 3-(am/5000);
  return_value = (fz-g*M-Cz*abs(k)*k)/M;
endfunction

function return_value = f2 (y,k,Tz) 
  global am;
  persistent Iz = 0.08;
  persistent Cy = 5 + (am/5000);
  return_value = (Tz - Cy*abs(k)*k)/Iz;
endfunction


function [t,z] = Euler(x0,k0,t_start,h,t_end,fun,input)
  
  t = (t_start:h:t_end); 
  
  %fprinf('ModifiedEuler\n')
  %fprinf('i   z   t   dz/dt\n')
    
  z = zeros(size(t));
  z(1) = x0;
  k(1) = k0;

  for i = 1:1:length(t) - 1
      z(i + 1) = z(i) + h * k(i);
      k(i+1) = k(i) + h*fun(z(i),k(i),input);
      %fprintf('%10d%+10.4f%+10.2f%+15.2f\n', i, z(i + 1), t(i + 1), fun(z(i + 1),k(i+1),input));
  end
endfunction 
 
function [tt,zz] = ModifiedEuler(x0,k0,t_start,h,t_end,fun,input)
    tt = (t_start:h:t_end);
    
    %fprinf('ModifiedEuler\n')
    %fprinf('i   z   t   dz/dt\n')
    
    zz = zeros(size(tt));
    kk = zeros(size(tt));
    zz(1) = x0;
    kk(1) = k0;
    for j = 1:1:length(tt) - 1
      zz(j + 1) = zz(j) + h * kk(j);
      kk(j+1) = kk(j) + h*fun(zz(j)+h/2,kk(j)+(h/2)*fun(zz(j),kk(j),input),input);
      %fprintf('%10d%+10.4f%+10.2f%+15.2f\n', j, zz(j+1), tt(j + 1), fun(zz(j)+h/2,kk(j)+(h/2)*fun(zz(j),kk(j),input),input));
    end
endfunction

%------------------Constant Vars---------------------
global am = 3106;
global M=1;
global g=9.81;
%------Initial Values,step and time------------------
z0 = am/1000;
y0 = 0;
t_start = 0;
k0 = 0;
h = 0.01; 
t_end = 30; 

%-------------------input fz = Mg - am/1000 Tz = 0------------------------------

%-----------Methods for the equation (1)-------------- 
[t,z] = Euler(z0,k0,t_start,h,t_end,@f1,M*g + am/1000);
[tt,zz] =ModifiedEuler(z0,k0,t_start,h,t_end,@f1,fz=M*g + am/1000);

%-----------Methods for the equation (2)--------------
[t_,y_] = Euler(y0,k0,t_start,h,t_end,@f2,0);
[t_t,y_y] =ModifiedEuler(y0,k0,t_start,h,t_end,@f2,0);

%-----------input fz = Mg  Tz = -am/10000---------------------------

%-----------Methods for the equation (1)-------------- 
[t1,z1] = Euler(z0,k0,t_start,h,t_end,@f1,M*g);
[tt1,zz1] =ModifiedEuler(z0,k0,t_start,h,t_end,@f1,M*g);
printSomeValues(z1,zz1)
%-----------Methods for the equation (2)--------------
[t2,y2] = Euler(y0,k0,t_start,h,t_end,@f2,-am/10000);
[tt2,yy2] =ModifiedEuler(y0,k0,t_start,h,t_end,@f2,-am/10000);
printSomeValues(y2,yy2)
%-----------generating the plots-------------------------
createPlot(t,z,'Z',"Euler's Method",'r',"fz = M*g-am/1000")
%print -djpeg Z_plot_for_the_first_fz.jpg
createPlot(tt,zz,'Z',"Modified Euler's Method",'b',"fz = M*g-am/1000")
%print -djpeg Z_plot_for_the_first_fzME.jpg
createPlot(t_,y_,"Y~(z's orientation)","Euler's Method",'g',"Tz = 0")
%print -djpeg Y_plot_for_the_first_fzE.jpg
createPlot(t_t,y_y,"Y~(z's orientation)","Modified Euler's Method",'k',"Tz = 0")
%print -djpeg Y_plot_for_the_first_fzME.jpg
createPlot(t1,z1,'Z',"Euler's Method",'m',"fz = M*g")
%print -djpeg Z_plot_for_the_second_fzE.jpg
createPlot(tt1,zz1,'Z',"Modified Euler's Method",'r',"fz = M*g")
%print -djpeg Z_plot_for_the_second_fzME.jpg
createPlot(t2,y2,"Y~(z's orientation)","Euler's Method",'b',"Tz = -am/10000")
%print -djpeg Y_plot_for_the_second_fzE.jpg
createPlot(tt2,yy2,"Y~(z's orientation)","Modified Euler's Method",'k',"Tz = -am/10000")
%print -djpeg Y_plot_for_the_second_fzME.jpg
