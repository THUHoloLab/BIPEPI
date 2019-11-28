function [out]=Weightfunction(w)
w=abs(w);
c=-0.5;
if (w<=1)&&(w>=0)
	out=1-(c+3)*w^2+(c+2)*w^3;%weight function
    else if (w>1)&&(w<=2)
        out=c*w^3-5*c*w^2+(8*c)*w-4*c;%weight function
    else
        out=0;
    end
end
end