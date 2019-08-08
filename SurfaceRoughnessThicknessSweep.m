% This program is meant to sweep through an array of thickness for the
% surface roughness and plot them.


T = 0.05:.10:1;% thickness array

att_constant_dB = -20*log10(att_constant);
% correct the attentuation constant with the roughness coefficient


Corrected_attenuations = [];


for idx = 1:numel(T)

    rms = T(idx);
    rms = rms.*1e-6
    D = 1./sqrt(s*pi*uo*ur*freq); % skin depth
    x = (rms./D).^2;
    Ksr = 1+(2*atan(1.4*x))./pi; % copper surface roughness correction factor
    
    Corrected_attenuations(:,idx) = att_constant_dB .* Ksr
    
    
end

f = freq./1e9

N = numel(T);

M = [1 2 3 4 5 6 7 8 9 10];
LegendString = cell(1,numel(T));


figure
hold all
box on
title('Corrected Attentuation Constant')
ylabel('dB')
xlabel('Frequency(GHz)')
xlim([1,20])
legendCell = cellstr(num2str(thickness', 'M=%-d'));
% legend(legendCell)
for idx2 = 1:N
    plot(f,Corrected_attenuations(:,idx2))
    LegendString{idx2} = sprintf('n = %i',T(idx2));
end


% x = 1:10;
% y = rand(1,10);
% 
% n=[2 4 6 8 10];

%// Initialize the cell containing the text. For each "n" there is a cell.
% LegendString = cell(1,numel(n));

%// Plot every curve and create the corresponding legend text in the loop.
% hold all
% for k = 1:numel(n)
% 
%     plot(x,n(k)*y)
%     LegendString{k} = sprintf('n = %i',n(k));
% end

%// Display the legend
legend(LegendString)

