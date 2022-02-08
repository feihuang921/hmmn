year = [1893:1908];
 F_data= importdata('p_f.txt');
 M_data= importdata('p_m.txt');
 
 
F_100=[381,428, 418, 436, 412, 502, 496, 482, 478, 529, 595, 581, 610, 639, 586, 657];
%F_100_point=[389.2, 406.5,435.7,441.7,477.8, 493.6,503.3,488.6,539.9,550.2,560.4,631.4,601.1,595.6,635,780];
%F_100_lower=[240.6, 249.8, 267.8, 271.5, 295.2,309.3,303.3, 300.3, 333.7, 338.2, 344.5, 393.1,371.8,366.3,392.9,488.8];
%F_100_upper=[465.6,489.8,525.1,532.2,572, 594.8,606.5,588.6,646,662.8,674.8,748.9,718.4,717,758.8,917.5];

 F_100_point=F_data(1,1:16);
 F_100_lower=[];
 for i=1:length(year)
     F_100_lower=[F_100_lower, CI_f(1, 1, i)];
 end
 
 F_100_upper=[];
 for i=1:length(year)
     F_100_upper=[F_100_upper, CI_f(1, 2, i)];
 end
 
 
F_105=[18,28,31, 17,22, 22,31,28,20,32,35,28,44,52,38,37];
%F_105_point=[28, 25.9, 28.2, 28.2, 36.7,32.0,32.6,29.8,39.3, 33,30.5,50,38.5,29.7,38.9,73.4];
%F_105_lower=[19.4,17.7,19.3,19.3,25.3,21.9,22.3,20.4,27.3,22.7,21.1,35.2,27,20.7,27.4, 51.8];
%F_105_upper=[29.6,27.7,30.2,30.1,39,34.2,34.9,31.7,41.7,35.1,32.3,52.6,40.6,31.4,40.9, 76.9];

 F_105_point=F_data(2,1:16);
 F_105_lower=[];
 for i=1:length(year)
     F_105_lower=[F_105_lower, CI_f(2, 1, i)];
 end
 
 F_105_upper=[];
 for i=1:length(year)
     F_105_upper=[F_105_upper, CI_f(2, 2, i)];
 end
 
M_100=[123,132, 108,116,103,96,104,111,108,110,103,119,112,99, 128, 99];
%M_100_point=[194.1,129.0,196.5,115.4, 108,104.4,115.7,107.7,123.9,110.5,98.8,143.6,113.3,122.1,116, 108.6];
%M_100_lower=[128.9,85.4,130.8,71.1,66.5,64.4,73.4,67.7,79.4,68.7,60.8,93.2,69.8,77.6,71.5, 68.3];
%M_100_upper=[212.1, 141.4,214.4,138.6,129.7, 125.4,134,126.2, 141.6, 131.3, 118.7, 161.6,136.1,141, 139.4,  127.1];

 M_100_point=M_data(1,1:16);
 M_100_lower=[];
 for i=1:length(year)
     M_100_lower=[M_100_lower, CI_m(1, 1, i)];
 end
 
 M_100_upper=[];
 for i=1:length(year)
     M_100_upper=[M_100_upper, CI_m(1, 2, i)];
 end

M_105=[5,5,5, 5, 2 ,3,4, 4, 3, 3,10,6,2, 6, 11,9];
%M_105_point=[24.4,10.2,24,3.4,3.4,2.9,5.3,3.8,6.9,3.5,3.1,8.6,3.6,5,3.3, 3.7];
%M_105_lower=[18.6,8,18.3,2.5,2.5,2.1,4.1,2.9,5.3,2.6,2.2, 6.7,2.6,3.9,2.4,2.8];
%M_105_upper=[25,10.4,24.5,3.5,3.6,3,5.5,3.9,7.1,3.6,3.2,8.8,3.7, 5.1,3.4,3.8];

 M_105_point=M_data(2,1:16);
 M_105_lower=[];
 for i=1:length(year)
     M_105_lower=[M_105_lower, CI_m(2, 1, i)];
 end
 
 M_105_upper=[];
 for i=1:length(year)
     M_105_upper=[M_105_upper, CI_m(2, 2, i)];
 end

 errorbar(year, F_100_point, F_100_point-F_100_lower, F_100_upper-F_100_point, 'o')
 xlim([1892 1909])
 ylim([0 1000])
 %ylabel('Observed and expected number alive, 100+, Female')
 %xlabel('Cohort')
 hold on;
 scatter(year, F_100, 'red');
 line(year, F_100, 'color','red');
 hold off;
% %  
%  errorbar(year, F_105_point, F_105_point-F_105_lower, F_105_upper-F_105_point, 'o')
%  xlim([1892 1909])
%  ylim([0 150])
% %ylabel('Observed and expected number alive, 105+, Female')
% %xlabel('Cohort')
%  hold on;
%  scatter(year, F_105, 'red');
%  line(year, F_105, 'color','red');
%  hold off;
% % % %  
%  errorbar(year, M_100_point, M_100_point-M_100_lower, M_100_upper-M_100_point, 'o')
%  xlim([1892 1909])
%  ylim([0 500])
%  %ylabel('Observed and expected number alive, 100+, Male')
%  %xlabel('Cohort')
%  hold on;
%  scatter(year, M_100, 'red');
%  line(year, M_100, 'color','red');
%  hold off;
%  
%  errorbar(year, M_105_point, M_105_point-M_105_lower, M_105_upper-M_105_point, 'o')
%  xlim([1892 1909])
%  ylim([0 100])
% %  ylabel('Observed and expected number alive, 105+, Male')
% %  xlabel('Cohort')
%  hold on;
%  scatter(year, M_105, 'red');
%  line(year, M_105, 'color','red');
%  hold off;