clear all;
close all;
clc;
%%
AA = ReadData('../../../Dropbox\Applied Statistics\Project 2\Data/Tjaereborg_20090901_1000_000035.tim',174);
%%
figure (2)
hold on
set(0,'defaultaxescolororder',winter(6))
for i = 113:118;
    h = plot(AA(:,i));
    
end