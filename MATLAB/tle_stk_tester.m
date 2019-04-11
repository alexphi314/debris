% Alex Philpott
r = [-1761336.083,5782317.269,-3609765.529]./1000; %km
v = [-2316.498,-4430.769,-5686.826]./1000; %km/s
num = 23224;
[line1,line2] = tle_from_stk(r,v,num);

fprintf('Got tle:\n%s\n%s\n',line1,line2);