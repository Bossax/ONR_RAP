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
hydrophone = "HEM";    % EDIT
% extract tx rx information
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,Vel,Range,act_arrival,est_arrival] = tx_rx_extraction_Oct(day,start_hour,end_hour,hydrophone);

% Travel time perterbation
ttp = (act_arrival - est_arrival)*3600*24*1000;     % ms
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
% 27 03:00 - 29 12:30 (CW square)
% truncate data
t_start = datenum('2018-10-27 03:00','yyyy-mm-dd HH:MM');
t_end = datenum('2018-10-30 15:30','yyyy-mm-dd HH:MM');
keep_ind = find((tx_t>= t_start)&(tx_t<=t_end));
predictorlen = length(keep_ind);

% create independent terms vales
pre_Range = Range(keep_ind)';
pre_Vel= Vel(keep_ind)';
pre_theta= cos(theta(keep_ind)/180*pi);

max_R = max(pre_Range);
max_V = max(pre_Vel);
max_theta = max(pre_theta);
% create cross term's values
pre_RV = pre_Range.*pre_Vel;
pre_Vthe = pre_Vel.*pre_theta;
pre_Rthe = pre_Range.*pre_theta;
pre_RVT = pre_Range.*pre_Vel.*pre_theta;

max_RV = max(pre_RV);
max_Vthe = max(pre_Vthe);
max_Rthe = max(pre_Rthe);
max_RVT = max(pre_RVT);
% predictor matrix
X = [pre_Range pre_Vel pre_theta pre_RV pre_Vthe pre_Rthe pre_RVT];
dimless_X = X./(ones(size(X,1),(size(X,2))).*[max_R max_V max_theta max_RV max_Vthe max_Rthe max_RVT]);
% response vector
y = ttp(keep_ind)';

%% Dimensionless equation
dimless_R = pre_Range./max_R;
dimless_Vtheta = pre_Vthe./max_Vthe;

a = 0.1:.1:18;
b =  0.1:.1:40;
ss = [];
for ii = 1:length(a)
    p_ttp = [];
    for jj =1:length(b)
        p_ttp = a(ii)*dimless_R + b(jj)*dimless_Vtheta;
        ss(end+1) = sum((p_ttp -ttp').^2);
    end
end
[MSE,I] =min(ss);
figure(2)
clf
plot(ss)
a_best = a(floor(I/length(b)))
b_est = b(I-length(b)*floor(I/length(b)))


%% 2. Regression
%dimless_ttp = ttp.*Vel./Range;
B = lasso(dimless_X,y);
[~,lambda_num] = size(B);
[~,num_term] = size(X);
% visualize
f1 = figure(1);
f1.Units = 'normalized';
f1.Position = [0 0.8 0.5 0.7];
clf
for k = 1:num_term
    f1;
    plot([1:lambda_num]./100,B(k,:));
    hold on
    text(10,max(B(k,:))+0.01,string(k))
end

title({'Lasso Regression','TTP = \beta_{1}*Range+\beta_{2}*Vel+\beta_{3}*cos(\Theta)+\beta_{4}*Range*Vel+\beta_{5}*Vel*cos(\Theta)+\beta_{6}*Range*cos(\Theta)+\beta_{7}*R*Vel*cos(\Theta)'},'Fontsize',12);
grid on
xlabel('Lasso Penalty Coefficient (\lambda)')
ylabel('Coefficient')
legend('1 Range','2 Vel','3 cos(\theta)','4 RV','5 Vcos(\theta)','6 Rcos(\theta)','7 RVcos(\theta)')

% Test the model with various lambda by test data and plot Sum Sqaured Error
% prep test data
test_ind =  1:length(ttp);
tes_len = length(test_ind);
% create independent terms vales
test_Range = Range(test_ind)';
test_Vel= Vel(test_ind)';
test_theta= theta(test_ind);
% create cross term's values
test_RV = test_Range.*test_Vel;
test_Vthe = test_Vel.*cos(test_theta/180*pi);
test_Rthe = test_Range.*cos(test_theta/180*pi);
test_RVT = test_Range.*test_Vel.*cos(test_theta/180*pi);

% test data matrix
T = [test_Range test_Vel test_theta test_RV test_Vthe test_Rthe test_RVT];

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
best_coeff = best_coeff(find(best_coeff >=0.01 ));
best_f = sprintf('TTP = %.4f Range+%.4f Vel*cos(\\Theta)+%.4fRange*cos(\\Theta)',best_coeff);
title({'Test Models', best_f},'Fontsize',12)

%% Residula sum squares of each term
a1 = 0.4867;
a2 = 0.3992;
a3 = 0.0614;

ttp1 = a1*pre_Range;
ttp2 = a2*pre_Vthe;
ttp3 = a3*pre_Rthe;
sq_err1 = sum((ttp'-ttp1).^2);
sq_err2 = sum((ttp'-ttp2).^2);
sq_err3 = sum((ttp'-ttp3).^2);

sq = [sq_err1 sq_err2 sq_err3];

bar(sq)
set(gca,'XTicklabel',{' 0.4867Range' ,'0.3992V*cos(\Theta)' ,'0.0614R*cos(\Theta)'})
grid on
title('Function =  0.4867Range + 0.3992V*cos(\Theta) + 0.0614R*cos(\Theta)');
ylabel('Residual Sum Square (ms^{2})')
set(gca,'fontsize',12)

