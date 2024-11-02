clc; close all; clear;


[qun, Qun] = trussglobal3D(4,15,6);

disp('qun (mm):')
disp(qun*1000)

disp('Qun (kN):')
disp(Qun)

diary('Print Out.txt')