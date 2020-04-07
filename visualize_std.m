% Opens up the file soyou can see the standard deviation
% Detailed explanation of this function.
function std_fig = visualize_std(amplifier_data)
std_fig = figure;
hold on
for i = 1:2:16
    
    plot(amplifier_data(i,:)+i*600)
    for j = 1:6
        plot(-1*j*std(amplifier_data(i,:)) * ones(length(amplifier_data(i,:)),1) + i*600)
    end
end
end

