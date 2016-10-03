function [OK] = test_mfun_rainread()


%% Test read
[T1,step1] = CD1_mfun_rainread('test1_mfun_rainread.ixx');
[T2,step2] = CD1_mfun_rainread('test2_mfun_rainread.ixx');





%% checks:
OK = 1;
if(step1~=60 || step2~=300)
    OK = 0;
end


end

