function lab2()

% REMINDER

[t,q] = Sinfn(1,1,4000);

figure;
plot(t,q);
title('Simple Sin Wave');


[t,q] = Sinfn(1,1,20);
figure;
stem(t,q);
title('Sampled Sin Wave');

% -------------------------------------------------------------EXERCISE UN

N=10;

x = Step(4,N);

figure;
stem(x);
title('4-Delayed Step Signal');

% RESPONSE OF THE GIVEN NON-CAUSAL SYSTEM TO THE DELAYED STEP SIGNAL
for i=2:N-1
   y(i) = x(i)/2 + x(i+1)/2
    
end

figure;
stem (y);
title('RESPONSE TO SYSTEM OF EX: 1.1');

% CAUSALITY PROPERTY DOES NOT STAND FOR THIS SYSTEM SINCE THE FIRST NON-ZERO
% VALUE OF THE RESPONSE APPEARS BEFORE THAT OF THE SIGNAL, MEANING THAT 
% THE SYSTEM RELIES ON THE FUTURE VALUES, MAKING IT ANTI-CAUSAL.
% 
% RESPONSE OF THE PROPOSED CAUSAL SIGNAL TO THE DELAYED STEP SIGNAL

for i=2:N-1
   y(i) = x(i)/2 + x(i-1)/2
    
end

figure;
stem (y);
title('RESPONSE TO SYSTEM OF CAUSAL SYSTEM');

% WE CAN SEE THAT THE FIRST VALUE OF THE RESPONSE OF THE CAUSAL SYSTEM
% STARTS FROM THE 5TH INDEX, WHICH CORRESPONDS TO THE FIRST NON-ZERO VALUE
% OF OUR DELAYED STEP SIGNAL

 
% ----------------------------------------------------------EXERCISE DEUX
% 
% 2.1
% 
% figure;
% d = prim(x);
% stem(d);
% title('4-DELAYED STEP RESPONSE TO PRIM FN');
% 
%     function r = prim(f)
%         w=0;
%         for e=1:length(f)
%             
%             f(e) = f(e) + w;
%             w=f(e)+w;
%         
%         end
%         figure;
%         stem(f);
%         r=f;
%     
%     end
% 
% THE PRIMITIVE OPERATOR IS NOT SABLE AS THE VALUES ARE CONTINUOUSLY
% INCREASING GOING TOWARDS INFINITY
% 
% 
% 
% 2.2
% 
% x= zeros(1,N);
% x(1) = 1;
% 
% figure;
% d = prim(x);
% stem(d);
% title('IMPULSE RESPONSE OF PRIM');
% disp(d);
% 
% THE IMPULSE RESPONSE OF PRIM OPERATOR IS ALSO UNSTABLE.


% 2.3

x= zeros(1,N);
x(2) = 1;

y = zeros(1,N);

for n=2:N

    y(n) = x(n) + 2*(y(n-1));

end

figure;
stem(y);
title('REPONSE TO SYSTEM IN EX: 2.3');

% NOT STABLE.

% 2.4


b = zeros(1,N);


for n=2:N

    b(n) = x(n) + (1/3)*(b(n-1));

end

figure;
stem(b);
title('REPONSE TO SYSTEM IN EX: 2.4');
%STABLE, AS THE OUTPUT HAS CONVERGED TO ONE.

% ----------------------------------------------------------EXERCISE TROIS

% 3.1
xa = [0 0 0 0 1 2 3 4 5 0 0 0 0 0 0 0 0 0 0]; 
xb = [0 0 0 0 0 0 0 0 0 4 3 2 1 0 0 0 0 0 0];

j = zeros(1,length(xa));
m = zeros(1,length(xb));

for v = 2:length(xa)-1
j(v) = 3*xa(v-1) - 2*xa(v) + xa(v+1);
end


for v = 2:length(xa)-1
m(v) = 3*xb(v-1) - 2*xb(v) + xb(v+1);
end


figure;
stem(j);
title('REPONSE OF xa TO SYSTEM IN EX: 3.1');

figure;
stem(m);
title('REPONSE OF xb TO SYSTEM IN EX: 3.1');

xc=xa+xb;

g = zeros(1,N);

for v = 2:length(xa)-1
g(v) = 3*xc(v-1) - 2*xc(v) + xc(v+1)
end


figure;
stem(g);
title('REPONSE OF xc TO SYSTEM IN EX: 3.1');

%AS WE CAN SEE FROM THE GRAPHS, Y(XA+XB) = Y(XA) + Y(XB). 
%THIS PROVES THE PROPERTY OF LINEAR SYSTEM.


xal = [0 0 0 0 0 0 1 2 3 4 5 0 0 0 0 0 0 0 0]; %xa delayed by 2 indices
xbl = [0 0 0 0 0 0 0 0 0 0 0 4 3 2 1 0 0 0 0]; %xb delayed by 2 indices


jl = zeros(1,length(xa));
ml = zeros(1,length(xa));

%calculating response of xal signal
for v = 2:length(xa)-1
jl(v) = 3*xal(v-1) - 2*xal(v) + xal(v+1);
end

%calculating response of xbl signal

for v = 2:length(xa)-1
ml(v) = 3*xbl(v-1) - 2*xbl(v) + xbl(v+1);
end

xcl = xal+xbl;

gl = zeros(1,length(xa));

%calculating response of xal + xbl signal
for v = 2:length(xa)-1
gl(v) = 3*xcl(v-1) - 2*xcl(v) + xcl(v+1)
end

figure;
stem(jl);
title('REPONSE OF 2-DELAYED XA');

figure;
stem(ml);
title('REPONSE OF 2 -DELAYED XB');

figure;
stem(gl);
title('REPONSE OF SUM OF 2-DELAYED XA AND 2-DELAYED XB');

if gl == (ml + jl)
   disp(' !first EQUAL SIGNALS! ');
else
    disp(' !first unnnnnnnnnEQUAL SIGNALS! ');
end


%AS WE CAN SEE FROM THE GRAPHS, Y(X(A-2)+X(B-2)) = Y(X(A-2)) + Y(X(B-2)). 
%THIS PROVES THE PROPERTY OF TIME INVARIANT SYSTEM.

% 3.3


y21 = zeros(1,length(xa));
for v = 2:length(xa)-1
y21(v) = (xa(v-1))^2 % PROPOSED NON-LTI SYSTEM
end


y22 = zeros(1,length(xa));
for v = 2:length(xa)-1
y22(v) = (xb(v-1))^2 % PROPOSED NON-LTI SYSTEM
end


y23 = zeros(1,length(xa));
for v = 2:length(xa)-1
y23(v) = (xc(v-1))^2 % PROPOSED NON-LTI SYSTEM
end

%figure;
%stem();
%title('REPONSE OF xc TO NON-LTI SYSTEM IN EX: 3.1');

if(y23 == (y21 + y22))
   disp(' !EQUAL SIGNALS! ');
else
    disp(' !unnnnnnnnnEQUAL SIGNALS! ');
end

% CLEARLY THIS SYSTEM IS LINEAR. MISSION FAILED. 

%-------------------------------------------------------------------LA FIN




















end