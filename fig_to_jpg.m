figs = openfig('C2H6_C3H8_C4H10_1.5_260.fig');
for K = 1 : length(figs)
   filename = 'C2H6_C3H8_C4H10_1.5_260.jpg';
   saveas(figs(K), filename);
end