%% Yashar Zafari - 99106209
% Bahman 1402 - Feb. 2024
% Youla-Kucera Solving function
% Test case script for example on Page 14 of 22-23 IMC Lesson
clear;clc;
s=tf('s');
G=1/((s-1)*(s-2));
[Ns,Xs,Ms,Ys]=Euclid(G,1);
syms a b c s
Q=(a*s+b)/(s+c); % Format of Q should be specified, **SYMBOLIC*
% Ramp following key=2
[K,Q,cd,~]=youlaKucera(Xs,Ms,Ns,Ys,Q,2);
%% Substitute c with the desired value
K=subs(K,c,2);
Q=subs(Q,c,2);
%% Transforming them into Transfer Function
K=sym2tf(K);
Q=sym2tf(Q);
% For using the transfer function, keep it in mind that
% we changed 's' from transfer function symbol to symbolic value in order
% to perform the calculation, thus we need to define it again:
s=tf('s');
%% Checking the performance
step(feedback(K*G,1)/s,1/s);
axis auto