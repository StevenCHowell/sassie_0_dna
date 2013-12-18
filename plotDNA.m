% $Id: plotDNA.m,v 1.2 2013-12-18 14:53:00 schowell Exp $
% script to plot the results of calculating the Re and Rg of DNA
%
% data columns: moves  rg/lp      re/lp           a       r
clc;clc;
% close all;

d.files = {'131216_100lpA',...
    '131216_200lpA',...
    '131216_20lpA',...
    '131216_40lpA',...
    '131216_60lpA',...
    '131216_6lpA',...
    '131216_80lpA',...
    '131217_2lpB',...
    '131217_2lpA',...   
    '131217_1lpA',...
    '131217_2lpC',...
    '131217_4lpA'};
%     '131216_400lpA',...
d.name = {'lp100',...
    'lp200a',...
    'lp020a',...
    'lp040a',...
    'lp060a',...
    'lp006a',...
    'lp080a',...
    'lp002b',...          
    'lp002a',...
    'lp001a',...
    'lp002c',...
    'lp004a'};
%     'lp400a',...
d.L = [100,...
    200,...
    020,...
    040,...
    060,...
    006,...
    080,...
    002,...          
    002,...
    001,...
    002,...
    004];
%     400,...
d.path = './validate/';
d.ext = '_dnaMoves.o';

rg = loadxy('extractedRg.csv');
[d.L,j] = sort(d.L);
d.files = d.files(j);
d.name = d.name(j);
lag = 100;
for i = 1:length(d.name)
    d.(d.name{i}).cmp = spline(rg(:,1),rg(:,2),d.L(i));
    d.(d.name{i}).all = loadxy([d.path,d.files{i},d.ext]);
    d.(d.name{i}).rg = d.(d.name{i}).all(:,[1,2]);
    d.(d.name{i}).re = d.(d.name{i}).all(:,[1,3]);
    d.(d.name{i}).a = d.(d.name{i}).all(:,4);
    d.(d.name{i}).iter = sum(sum(d.(d.name{i}).all(:,[4,5])));
    d.(d.name{i}).accep = sum(d.(d.name{i}).a)/d.(d.name{i}).iter;
    d.(d.name{i}).mav = d.(d.name{i}).rg;
    n = length(d.(d.name{i}).rg);
    if lag > n
        d.(d.name{i}).mav(:,2) = tsmovavg(d.(d.name{i}).rg(:,2),'s',ceil(n/4),1);
    else
        d.(d.name{i}).mav(:,2) = tsmovavg(d.(d.name{i}).rg(:,2),'s',lag,1);
    end
    mySubplot;
    hold all;
%     xyplot(rg)
    xyplot(d.(d.name{i}).rg,'s')
    plot(d.(d.name{i}).rg(:,1),d.(d.name{i}).cmp*ones(length(d.(d.name{i}).rg),1),'-','linewidth',2);
    xyplot(d.(d.name{i}).mav,'.','linewidth',2)
    title(['Rg of',d.name{i},'  (a=',num2str(d.(d.name{i}).accep),')']);
    xlabel('iterations')
    ylabel('Rg/lp')
end
legend('Rg','actual Rg','Moving Av','location','NorthEastOutside')
legend boxoff
