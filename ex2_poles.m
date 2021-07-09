%-------------------------Exercise 2(B)-----------------------------------------
%-------------caclculate the change of poles by altering Kpz and Kdz------------
%-------------------------------------------------------------------------------

clear; clc; close('all');
function createPlot(a,Name,input)
    figure('name',Name);
    plot(real(a), imag(a),'x');
    xlabel ("Re(Z)");
    ylabel (["Im(Z)"]);
    title ([" PZ map for variable " input]);
    legend(Name);
endfunction


function [s1,s2] = calculateRoots(a,b)
    delta = (b^2) - (4*a);
    s1 = (-b + sqrt(delta))/(2*a);
    s2 = (-b - sqrt(delta))/(2*a);
endfunction

function result = findPoles(Kpz,Kdz,Cz,M,result,ch_var,end_value)
  start_value = 0.1;
  while(start_value<=end_value)
    if strcmp(ch_var,'Kpz')
      [s1,s2] = calculateRoots(M/start_value,(Kdz+Cz)/start_value);
    else
       [s1,s2] = calculateRoots(M/Kpz,(start_value+Cz)/Kpz);
    endif
    if s1==s2
      display(Kpz)
    endif
    result=[result;complex(s1)];
    result=[result;complex(s2)];
    start_value = start_value+0.1;
  endwhile
endfunction
%--------initial values--------------------
Kdz = 15;
Kpz = 5;
am = 3106;
Cz = 3+ am/5000;
M=1;
end_value = 1000;
%--------------constant Kdz variable Kpz---------------------------------
variableKpzPoles = findPoles(Kpz,Kdz,Cz,M,[],'Kpz',end_value);
createPlot(variableKpzPoles,'poles','Kpz')

%--------------constant Kpz variable Kdz---------------------------------
variableKdzPoles = findPoles(Kpz,Kdz,Cz,M,[],'Kdz',end_value);
createPlot(variableKdzPoles,'poles','Kdz')
