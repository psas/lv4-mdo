clear
clc
%  Ian Zabel
%% Provided Parameters
%  Objective Function: Minimize f(x1,x2)=(1-x1)^2+(2-x2)^2

Y = 5;     % Gamma:   Y > 1
B = 0.8;   % Beta:    0 < B < 1
E = 0.002; % Epsilon: Termination Parameter
a = 3;     % Scale Factor a
N = 2;     % Number of Variables N

d1 = a*(sqrt(N+1)+N-1)/(N+1); % Delta 1
d2 = a*(sqrt(N+1)-1)/(N+1);   % Delta 2

f=[];x1n=[];x2n=[];c=0; % Empty Matrices and Counter Value

%% Guesses
xg1 = [1 5];                   % Initial Guess
xg2 = [xg1(1)+d1,xg1(2)-d1];   % Second Guess
xg3 = [xg1(1)+d2,xg1(2)+d2];   % Third Guess
x1 = [xg1(1) xg2(1) xg3(1)];   % All Guess Values, x1
x2 = [xg1(2) xg2(2) xg3(2)];   % All Guess Values, x2
fg1 = (1-x1(1))^2+(2-x2(1))^2; % Function Guess 1
fg2 = (1-x1(2))^2+(2-x2(2))^2; % Function Guess 2
fg3 = (1-x1(3))^2+(2-x2(3))^2; % Function Guess 3
f   = [fg1 fg2 fg3];           % All Function Values
fc  = (1-(sum(x1)/N))^2 + (2-(sum(x2)/N))^2;
fb  = min(f);

%% Loop
while c < 1000
    x1b=[];x2b=[];
    f = [fg1 fg2 fg3];
for i = 1:3 % Populate x = [best, second best, worst]
    if f(i) == min(f)
      x1b(1) = x1(i); % Best x1
      x2b(1) = x2(i); % Best x2   
    elseif f(i) == max(f)
      x1b(3) = x1(i); % Worst x1
      x2b(3) = x2(i); % Worst x2   
    else
      x1b(2) = x1(i); % Second Best x1
      x2b(2) = x2(i); % Second Best x2
    end
end
    xw = [x1b(3) x2b(3)];                   % Worst Point
    xs = [x1b(2) x2b(2)];                   % Second Best Point
    xb = [x1b(1) x2b(1)];                   % Best Point
    xc = [(xs(1)+xb(1))/N (xs(2)+xb(2))/N]; % Center Point
    xr = [(2*xc(1)-xw(1)) (2*xc(2)-xw(2))]; % Reflection Point
    fw = (1-xw(1))^2+(2-xw(2))^2;           % Worst Function Point
    fs = (1-xs(1))^2+(2-xs(2))^2;           % Second Best Function Point
    fb = (1-xb(1))^2+(2-xb(2))^2;           % Best Function Point
    fc = (1-xc(1))^2+(2-xc(2))^2;           % Center Function Point
    fr = (1-xr(1))^2+(2-xr(2))^2;           % Reflection Function Point
    if fr < fb
      xn = (1+Y)*xc-Y*xw; % Expansion
    elseif fr >= fw
      xn = (1-B)*xc+B*xw; % Contraction
    elseif fs < fr && fr < fw
      xn = (1-B)*xc+B*xw; % Contraction
    end
    fn = (1-xn(1))^2+(2-xn(2))^2; % New Minimum Function
    x1 = [xb(1) xs(1) xn(1)];     % New x1 Matrix
    x2 = [xb(2) xs(2) xn(2)];     % New x2 Matrix
    fg3=fn;fg1=fb;fg2=fs;         % Repopulate f Matrix
    c=c+1;                        % Iterate Counter
    if sqrt((fg3-fc)^2)/(N+1) <= E
      break % Terminate Condition Check
    end
end
fprintf('Iterations: %3.0f \n\n',c)
fprintf('Minimum Points:\n')
fprintf('x1: %3.3f\n',xn(1))
fprintf('x2: %3.3f\n\n',xn(2))
fprintf('f(x) at Minimum: %3.3f \n',fn)
