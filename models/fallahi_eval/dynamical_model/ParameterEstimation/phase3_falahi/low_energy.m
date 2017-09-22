clear param
j=1;
minenergy = 1;
while j <= 100
    [c,i] = min(energy_chain(1,1:750));
    param(:,j) = params_chain(1,:,i);
    index(1,j) = i;
    energy(1,j) = energy_chain(1,i);
    energy_chain(1,i) = 100000;
    j = j+1;
end

% param_best = param(:,1:10);
% param_worst = param(:,90:1000);