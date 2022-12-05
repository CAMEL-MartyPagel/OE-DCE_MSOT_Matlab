function [ places ] = choose_timeproj( M,over_fact,wl_num )  
% places = choose_timeproj(M,over_fact,wl_num)
% ------------------------------------------------------------------------
% returns locations of relevant times and projections
%
% Parameters:
% - M:          Matrix
% - over_fact:  Threshold
% - wl_num:     ???
    MM=max(abs(M),[],2);
    ele_num=over_fact*wl_num;
    M_s=sort(MM,'descend');
    thresh=M_s(ele_num);
    places=find(MM>thresh)';
end