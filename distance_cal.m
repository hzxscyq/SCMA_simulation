function [dis]=distance_cal(C1,C2,C3)
len = length(C1);
count=1;
for i=1:len
    for j=1:len
        for k=1:len
            d1=real(C1(i)-C2(j))^2+imag(C1(i)-C2(j))^2;
            d2=real(C1(i)-C3(k))^2+imag(C1(i)-C3(k))^2;
            d3=real(C2(j)-C3(k))^2+imag(C2(j)-C3(k))^2;
            d_min(count)=min([d1,d2,d3])
            count=count+1;
        end
    end
end
dis=min(d_min);

end




