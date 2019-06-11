function [c]=cssprofile( z, lat);
%function [c]=cssprofile( z);
% canonical sound speed profile
%
% z (meters)
% c (meters/second)
%
% if abs(lat)<67.5 temperate ocean (default)
% else polar ocean
%

gammaa=0.0113/1000;
za=1000;
h=1000;
ca=1500;

if nargin==1 lat=0; end

if abs(lat)<67.5
  %temperate ocean
  eta=2*(za-z)./h;
  c=ca*(1+h*gammaa*(exp(eta) - eta - 1)/2); 
else
  %polar ocean
  c=ca*(1+gammaa*z);
end

