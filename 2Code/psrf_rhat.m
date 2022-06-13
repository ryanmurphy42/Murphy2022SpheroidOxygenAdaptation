function [Rhat,Vh,W,B] = function_psrf_rhat(varargin)
% Potential scale reduction factor (also referred to as Rhat)
% based on Gelman--Bayesian Data Analysis--3rd edition, CRC Press

% assumes all chains the same size and length of chain is an even number.

% 1) split each chain in half
chain_tmp = varargin{1};
chain_tmp_size = size(chain_tmp);
N =  chain_tmp_size(1)/2;% length of chain
M = length(varargin)*2; % number of chains
P = chain_tmp_size(2);% number of parameters

chains_matrix = zeros(N,P,M);

chain_counter = 0;
for i=1:length(varargin)
    chain_tmp = varargin{i};
    chain_counter = chain_counter +1;
    chains_matrix(1:N,1:P,chain_counter) = chain_tmp(1:N,:);
    chain_counter = chain_counter +1;
    chains_matrix(1:N,1:P,chain_counter) = chain_tmp((N+1):end,:);
end


if N>2

    % 2) Calculate B: between-sequence variance and calculate W: within-sequence variance
    B = zeros(1,P);
    W = zeros(1,P);
    for p_loop=1:P
        phi_dot_j = zeros(1,M); % mean per chain
        for b_loop=1:M
            phi_dot_j(b_loop) = (1./N).*sum(chains_matrix(:,p_loop,b_loop));
            s_j_squared(b_loop) = (sum((chains_matrix(:,p_loop,b_loop)-phi_dot_j(b_loop)).^2));
        end
        phi_dotdot = (1./M).*sum(phi_dot_j); % mean of mean per chain
        B(p_loop) = (N./(M-1)).*sum((phi_dot_j - phi_dotdot).^2); % between-sequence variance
        W(p_loop) = (1./(N-1)).*(1./M)*sum(s_j_squared);
    end

    % 3) calculate Rhat
    var_phi_y = (N-1)/N*W + B./N;
    Rhat = sqrt(var_phi_y./W);

else
    Rhat = zeros(length(varargin),1); % not calculated as only 1 chain.
end

end