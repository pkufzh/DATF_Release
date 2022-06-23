%% 程序初始化
clear all
clf
clc

%% 翼型斜率计算
X1=0.898634206312060;
X2=0.901105671280223;
Y1=0.019093445326927;
Y2=0.0186601365338074;
Xp=0.9;
Yp=interp1([X1,X2],[Y1,Y2],Xp,'linear');
Lp=0.30;
beta=atan((X2-X1)/(Y1-Y2));
Xt=Xp+Lp*cos(beta)
Yt=Yp+Lp*sin(beta)
