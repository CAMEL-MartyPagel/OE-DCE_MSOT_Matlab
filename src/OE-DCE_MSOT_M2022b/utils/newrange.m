function myImgNorm = newrange(myImg,myRange)


%% Normalize the Image:
%myRange = getrangefromclass(myImg(1));
newMax = myRange(2);
newMin = myRange(1);

myImgNorm = (myImg - min(myImg(:)))*(newMax - newMin)/(max(myImg(:)) - min(myImg(:))) + newMin;

end