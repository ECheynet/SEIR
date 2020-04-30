function r = rNumber(infected)
    shifted = infected(2: length(infected));
    infected = infected(1: length(infected) - 1);
    r = infected .\ shifted;
end