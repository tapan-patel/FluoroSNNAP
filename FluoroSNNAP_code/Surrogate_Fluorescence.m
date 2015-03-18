function F = Surrogate_Fluorescence(spks,frames,fps)
% Given spike times (array), return simulated fluorescence trace
% [Ca2+]t - [Ca2+]_t-1 = - deltaT/tau_Ca* [Ca2+]_t-1 + A*n
% where [Ca2+]t is the calcium concentration at time t, and [Ca2+]t-1 is
% the calcium concentration at time t-1. deltaT is time-step of image
% acqusition and tau_Ca is decay constant, 1s. A is the increase in Ca2+
% after a single spike, set to 50uM and n is the number of spikes at time
% t.

% Input: spks = array of spiking frames
% frames = total number of frames in image stack
% fps = imaging frame rate
% For example, a neuron fired at frames 10, 40, 70 and 130. total 200
% frames and frame rate was 5
% F = Surrogate_Fluorescence([10,40,70,130],200,5);

% This function is part of the bigger, FluoroSNNAP, package. Please see
% http://www.seas.upenn.edu/~molneuro/fluorosnnap.html for more details
% Author: Tapan P Patel, PhD, tapan.p.patel@gmail.com

deltaT = 1/fps;
tau = 1;
A = 50;
Ca = zeros(frames,1);
for i=2:frames
    Ca(i) = -deltaT/tau*Ca(i-1) + Ca(i-1);
    if(nnz(i==spks))
        Ca(i) = Ca(i) + A;
    end
end

% Fluorecence = [Ca2+]/([Ca2+]+Kd) + mu
% where Kd for fluorophore is 300 and mu is a normally distributed gaussian
% noise with standard deviation of 0.001 and mean 0
F = Ca./(Ca+300)+normrnd(0,0.001,[frames,1]);