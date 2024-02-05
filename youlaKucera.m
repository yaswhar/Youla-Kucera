function [K,Q,cd,cden]=youlaKucera(Xs,Ms,Ns,Ys,Q,key)
% Inputs: Xs, Ms, Ns, Ys from Euclid Fucntion in TF fromat, Q in Symbolic
% Format, and key as the number of integrators we want K to have 
%
% Outputs: K, Q, cd *The coefficients of denominator which should be zero
% to generate the integrators*, cden *The coefficients of the denominator*
%
% Currently, I couldn't automate it for the key above 2. I will try to find
% a solution for that too.
%
% Course: Foundamentals of Automatic Control Design. 28-255 , Semester:
% 1402-1
% Sharif University of Technology, Department of Mechanical Engineering, Tehran, Iran.
% Prepared by: Yashar Zafari, S.N.: 99106209
syms s
Ns = tf2sym(Ns);
Xs = tf2sym(Xs);
Ms = tf2sym(Ms);
Ys = tf2sym(Ys);
K=simplifyFraction((Xs+Ms*Q)/(Ys-Ns*Q));
pretty(K)
[~,den]=numden(K);
[cden,~]=coeffs(den,s,"All");
disp('Coeffs to be zero:')
cd=cden(end-key+1:end)
eq=cd==0;
switch key
    case 0
        return
    case 1
        svars=symvar(cd);
        sol=solve(eq,svars(1));
        K=simplifyFraction(subs(K,svars(1),sol));
        Q=simplifyFraction(subs(Q,svars(1),sol));
    case 2
        svars=symvar(cd);
        com=intersect(symvar(cd(1)),symvar(cd(2)));
        sol=struct2cell(solve(eq,setdiff(svars,com)));
        K=simplifyFraction(subs(K,setdiff(svars,com),[sol{1} sol{2}]));
        Q=simplifyFraction(subs(Q,setdiff(svars,com),[sol{1} sol{2}]));
end
disp('K=')
pretty(K)
disp('Q=')
pretty(Q)
end