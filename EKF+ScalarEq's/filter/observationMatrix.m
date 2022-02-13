function H = observationMatrix(attitude, Borb)
Jh = jacobiH(attitude, Borb);
H = Jh;
end

