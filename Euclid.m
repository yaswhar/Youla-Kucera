function [Ns, Xs, Ms, Ys] = Euclid(G, alpha)
% Implement Euclid's algorithm that is needed in Youla-Kucera Theorem.

Ns = []; Xs = []; Ms = []; Ys = [];
syms s lambda;

% Check whether G is asymptotically stable or not
if isstable(G)
    Ns = G;
    Xs = sym2tf(0*s+0);
    Ms = sym2tf(0*s+1);
    Ys = sym2tf(0*s+1);
    disp("G is asymptotically stable");
    return;
end

% Transform s to lambda
[Num, Den] = tfdata(G, 'v');
G_syms = poly2sym(Num, s) / poly2sym(Den, s);
G_lambda = subs(G_syms, s, (1-alpha*lambda)/lambda);
[N, M] = numden(G_lambda);
n = sym2poly(N);
m = sym2poly(M);

% Check whether degree(n(lambda)) is greater that degree(m(lambda))
if length(n)-1 < length(m)-1
    temp = n;
    n = m;
    m = temp;
end

% Initialization
num = 1000; % maximum number of steps
q = cell(num, 1);
r = cell(num, 1);

% Step 1
[q{1}, r{1}] = deconv(n, m);
r{1} = nonzeros(r{1})';
showSteps(1, n, m, q{1}, r{1});

if length(r{1}) > 1
    % Step 2
    [q{2}, r{2}] = deconv(m, r{1});
    r{2} = nonzeros(r{2})';
    showSteps(2, m, r{1}, q{2}, r{2});
    
    if length(r{2}) > 1
        % Step i
        i = 2;
        while length(r{i}) > 1
            [q{i+1}, r{i+1}] = deconv(r{i-1}, r{i});
            r{i+1} = nonzeros(r{i+1})';
            i = i + 1;
        end
        showSteps(i, r{i-2}, r{i-1}, q{i}, r{i});
        
        % Solve equations as demonstrated in page 11 of "22-23 IMC" handout
        qsym = sym('q', [1, i]);
        A = sym(eye(i));
        A(2, 1) = qsym(2);
        for j = 3:i
            A(j, j-2) = -1;
            A(j, j-1) = qsym(j);
        end
        B = sym(zeros(i, 1));
        syms nsym msym;
        B(1) = nsym - msym * qsym(1);
        B(2) = msym;
        Roots = inv(A) * B;
        [temp1, temp2] = coeffs(Roots(i), [nsym, msym]);
        xsym = temp1(1);
        X = xsym;
        ysym = temp1(2);
        Y = ysym;
        for j = 1:i
            X = subs(X, qsym(j), poly2sym(q{j}, lambda));
            Y = subs(Y, qsym(j), poly2sym(q{j}, lambda));
        end
        rlast = poly2sym(r{i}, lambda);
        [N, X, M, Y] = remainder(rlast, N, X, M, Y);
    else
        q1 = poly2sym(q{1}, lambda);
        q2 = poly2sym(q{2}, lambda);
        r2 = poly2sym(r{2}, lambda);
        X = -q2;
        Y = (1+q1*q2);
        [N, X, M, Y] = remainder(r2, N, X, M, Y);
    end
else
    q1 = poly2sym(q{1}, lambda);
    r1 = poly2sym(r{1}, lambda);
    X = 1;
    Y = -q1;
    [N, X, M, Y] = remainder(r1, N, X, M, Y);
end

% Results
Ns = subs(N, lambda, 1/(s+alpha));
Ns = sym2tf(Ns);
Xs = subs(X, lambda, 1/(s+alpha));
Xs = sym2tf(Xs);
Ms = subs(M, lambda, 1/(s+alpha));
Ms = sym2tf(Ms);
Ys = subs(Y, lambda, 1/(s+alpha));
Ys = sym2tf(Ys);
end

function [N, X, M, Y] = remainder(rlast, N, X, M, Y)
% Check the sign of the remainder

if rlast > 0
    X = X / rlast;
    Y = Y / rlast;
elseif rlast < 0
    rlast = abs(rlast);
    N = -N;
    M = -M;
    X = X / rlast;
    Y = Y / rlast;
else
    X = 0;
    Y = 1;
end
end

function showSteps(numStep, a, b, c, d)
% Show steps in the required fromat

aname = inputname(2);   
bname = inputname(3);
if isequal(aname, 'n') && isequal(bname, 'm')
    cname = 'q1';
    dname = 'r1';
elseif isequal(aname, 'm')
    bname = 'r1';
    cname = 'q2';
    dname = 'r2';
else
    aname = strjoin(["r", num2str(numStep-2)], "");
    bname = strjoin(["r", num2str(numStep-1)], "");
    cname = strjoin(["q", num2str(numStep)], "");
    dname = strjoin(["r", num2str(numStep)], "");
end
temp = strjoin([num2str(numStep), ": ", aname, " = ", bname, " * ", cname, " + ", dname], "");
disp(temp);

astr = makeSteps(a);
bstr = makeSteps(b);
cstr = makeSteps(c);
dstr = makeSteps(d);
temp = strjoin([astr, " <strong>=</strong> ", bstr, " <strong>*</strong> ", cstr, " <strong>+</strong> ", dstr], "");
disp(temp);
fprintf("\n");
end

function str = makeSteps(var)
% Produce the required format of each step

lambda = char(955);
str = strings;
j = 1;

for i = length(var):-1:1
    formatSpec1 = "(%s^%d)";
    if i > 1
        str1 = compose(formatSpec1, lambda, i-1);
    else
        str1 = num2str(1);
    end
    temp = arrayfun(@string, poly2sym(var(j)));
    formatSpec2 = "%s*%s";
    str2 = compose(formatSpec2, temp, str1);
    j = j + 1;
    str = [str, str2, " + "];
end

str = str(1:end-1);
str = ["<strong>[</strong>", str, "<strong>]</strong>"];
str = strjoin(str, "");
end

function Gtf = sym2tf(Gsym)
% Transform symbolic expressions to transfer functions
[symNum, symDen] = numden(Gsym);
TFnum = sym2poly(symNum);  
TFden = sym2poly(symDen);
Gtf = tf(TFnum, TFden);
end
