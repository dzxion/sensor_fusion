% TestScript.m
%
% This script tests the quaternion library functions to ensure that each
% function output is consistent.
%
% Date          Author          
% 2021/11/09    Deng zhengxiong 

%% Start of script

close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal

global ToRad ToDeg

% q = [cos(pi/8),        0,        0,sin(pi/8);
%      cos(pi/8),        0,sin(pi/8),        0;
%      cos(pi/8),sin(pi/8),        0,        0];

q = [cos(pi/8),        0,0,        sin(pi/8)];
%% Quaternion to ZYX Euler angles
euler = quatern2euler(q) * ToDeg;

%% Quaternion to rotation matrix
acc_w = [0 0 -1]';
R = quatern2rotMat(q);
acc_w = R * acc_w

%% End of script