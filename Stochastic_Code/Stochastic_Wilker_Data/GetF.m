function F = GetF(k, l, fitness_m)

if (k+l) == 0
    F = 0; return;
end

fitness_w = 1;
F = (k/(k+l))*fitness_m + (l/(k+l))*fitness_w;
