% Regression Analysis
% variable = Range,Vel,theta
% functional form
% ttp = A*Range + B*Vel +C*cos(theta)+D*Range*Vel+E*Vel*cos(theta)*F*Range*cos(theta)+G*R*Vel*cos(theta)
% Matlab Lasso Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%ACO LAT/LON
HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018

% icListen LAT/LON
icListen_lat = 22.739153;
icListen_lon = -158.0061254;

day = 27:30 ;               %  Edit
start_hour = 3;             % Edit
end_hour = 14;              % EDIT
hydrophone = "icListen";    % EDIT
% extract tx rx information
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(day,start_hour,end_hour,hydrophone);

% Travel time perterbation
ttp = (act_arrival - real(est_arrival))*3600*24*1000;     % ms
%% Calculate Azimuth and Heading relative to the ACO
azmth = ones(length(tx_lat),1);
% 0-360
for i=1:length(tx_lat)
    azmth(i) = azimuth(icListen_lat,icListen_lon,tx_lat(i),tx_lon(i));
end
theta = tx_heading';
theta = tx_heading' - azmth;

for i = 1:length(theta)
   if theta(i) < 0 
       theta(i) = theta(i)+360;
   end
    
end

%% Regression analysis
% 1. Prep predictor matrix
% truncate data
t_start = datenum('2018-10-27 03:00','yyyy-mm-dd HH:MM');
t_end = datenum('2018-10-30 15:30','yyyy-mm-dd HH:MM');
keep_ind = find((tx_t>= t_start)&(tx_t<=t_end));
predictorlen = length(keep_ind);

% create independent terms vales
pre_Range = range(keep_ind)';
pre_Vel= tx_xvel(keep_ind)';
pre_theta= theta(keep_ind);
pre_azmth = cos(azmth(keep_ind)/180*pi);

% create cross term's values
pre_RV = pre_Range.*pre_Vel;
pre_Vthe = pre_Vel.*cos(pre_theta/180*pi);
pre_Rthe = pre_Range.*cos(pre_theta/180*pi);
pre_RVT = pre_Range.*pre_Vel.*cos(pre_theta/180*pi);
pre_theta= cos(theta(keep_ind)/180*pi);
% predictor matrix
X = [pre_Range pre_Vel pre_theta pre_azmth pre_Vthe pre_RV  pre_Rthe pre_RVT];

% response vector
y = ttp(keep_ind)';

%2. Regression
B = lasso(X,y);
[~,lambda_num] = size(B);
[~,num_term] = size(X);
% visualize
f1 = figure(1);
f1.Units = 'normalized';
f1.Position = [0 0.8 0.5 0.7];
clf
for k = 1:num_term
    f1;
    plot([1:lambda_num]./100,B(k,:),'LineWidth',3);
    hold on
    text(10,max(B(k,:))+0.01,string(k))
end

title({'Lasso Regression','TTP = \beta_{1}*Range+\beta_{2}*Vel+\beta_{3}*cos(\Theta)+\beta_{4}*Azimuth+\beta_{5}*Vel*cos(\Theta)+\beta_{6}*Range*Vel+\beta_{7}*Range*cos(\Theta)+\beta_{8}*R*Vel*cos(\Theta)'},'Fontsize',14);
grid on
xlabel('Lasso Penalty Coefficient (\lambda)','fontsize',14)
ylabel('Coefficient','fontsize',14)
legend({'1 Range','2 Vel','3 Relative Heading cos(\theta)','4 Azimuth','5 Vcos(\theta)','6 RV','7 Rcos(\theta)','8 RVcos(\theta)'},'FontSize',14,'Location','northeast')
set(gca,'fontsize',14)
% Test the model with various lambda by test data and plot Sum Sqaured Error
% prep test data
test_ind =  1:length(ttp);
tes_len = length(test_ind);
% create independent terms vales
test_Range = range(test_ind)';
test_Vel= tx_xvel(test_ind)';
test_theta= theta(test_ind);
test_azmth = azmth(test_ind);
% create cross term's values
test_RV = test_Range.*test_Vel;
test_Vthe = test_Vel.*cos(test_theta/180*pi);
test_Rthe = test_Range.*cos(test_theta/180*pi);
test_RVT = test_Range.*test_Vel.*cos(test_theta/180*pi);

% test data matrix
T = [test_Range test_Vel test_theta test_azmth test_RV test_Vthe test_Rthe test_RVT];

% measurement y0
y0 = ttp(test_ind)';

% Cal SSE
SSE = zeros(1,lambda_num);

% test all possible models
for jj = 1:lambda_num
    Coeff = B(:,jj);
    sq_err = 0;
    for ii = 1:tes_len
        predicted_y = Coeff'*T(ii,:)';
        sq_err = sq_err+ (predicted_y-y0(ii))^2;
    end
    SSE(jj) = sq_err/tes_len;
end
MSE = min(SSE); 
min_lambda = find(SSE == MSE);
% min_lambda = 69;
figure(2)
clf
plot([1:lambda_num]/100,SSE)
hold on
scatter(min_lambda/100,MSE,'r')
grid on
xlabel('\lambda')
ylabel('Sum Squared Error (TTP)')
text(min_lambda/100-0.1,MSE+5,{['\lambda = ' num2str(min_lambda/100)], ['MSE = ' num2str(MSE)]},'EdgeColor','k')
best_coeff = B(:,min_lambda);
best_coeff = best_coeff(find(abs(best_coeff) >=0.01 ));
best_f = sprintf('TTP = %.4f Range+%.4f Vel*cos(\\Theta)+%.4fRange*cos(\\Theta)',best_coeff);
title({'Test Models', best_f},'Fontsize',12)
set(gca,'Yscale','log')
ylim([1e-1 1e3])
%% SSE for the whole dataset

% Cal SSE
sq_err = 0;

if length(best_coeff) < 3
    best_coeff = [best_coeff;0;0];
end
for jj = 1: length(ttp)
%         predicted_y = best_coeff'*[Range(jj);Vel(jj)*cos(theta(jj)/180*pi);Range(jj)*cos(theta(jj)/180*pi)];
        predicted_y = best_coeff'*[range(jj);tx_xvel(jj)*cos(theta(jj)/180*pi);range(jj)*cos(theta(jj)/180*pi)];
        sq_err = sq_err+(predicted_y-ttp(jj))^2;

end
MSE_all = sq_err/length(ttp)







