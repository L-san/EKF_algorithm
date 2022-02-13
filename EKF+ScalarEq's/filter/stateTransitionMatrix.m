function F = stateTransitionMatrix(attitude,dt)
Jf = jacobiF(attitude,dt);
F = Jf;
end