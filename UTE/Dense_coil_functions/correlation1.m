function [ r ] = correlation1( vec_a, vec_b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here




vec_abs_a(:,1)=squeeze(abs(vec_a));

vec_abs_b(:,1)=squeeze(abs(vec_b));

n=size(vec_abs_a,1);

haut=n*sum(vec_abs_a(:,1).*vec_abs_b(:,1))-(sum(vec_abs_a(:,1))*sum(vec_abs_b(:,1)));

bas1=sqrt(n*sum(vec_abs_a(:,1).*vec_abs_a(:,1))- (sum(vec_abs_a(:,1))*sum(vec_abs_a(:,1))));

bas2=sqrt(n*sum(vec_abs_b(:,1).*vec_abs_b(:,1))- (sum(vec_abs_b(:,1))*sum(vec_abs_b(:,1))));

r=haut/(bas1*bas2);

msg_str=sprintf('haut: %f   bas1: %f   bas2: %f   r: %f \n' , haut, bas1, bas2, r);
%disp(msg_str);



end

