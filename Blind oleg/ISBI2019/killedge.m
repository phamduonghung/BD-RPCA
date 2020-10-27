function [g,h] = killedge(g,h)

% KILLEDGE Taper edges of complex-valued images.
%   KILLEDGE(g,h) extends the functionality of 'edgetaper' to the case of
%   complex-valued images 'g' and PSFs 'h'.
%
%   Written by O. Michailovich, 07/2018. 

if ~isreal(h)
    h=abs(h);
end

if isreal(g)
    g=edgetaper(g,h);
else
    gr=edgetaper(real(g),h);
    gi=edgetaper(imag(g),h);
    g=complex(gr,gi);
end

end