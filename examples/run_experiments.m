% Script to reproduce the numerical experiments of [CCZ22, Sec. 4]
%
% List of experiments:
% (1) Validation of the code
% (2) Evolutionary advection--diffusion--reaction equation
% (3) Allen--Cahn equation
% (4) Brusselator equation
%
% [CCZ22] M. Caliari, F. Cassini, and F. Zivcovich,
%         A mu-mode approach for exponential integrators: actions of
%         phi-functions of Kronecker sums, Submitted, 2022

clc
clear all
close all

code_validation;
example_adr;
example_allencahn;
example_brusselator;
