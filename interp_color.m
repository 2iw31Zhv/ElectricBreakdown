function v = interp_color(x, xmin, xmax)
xs = linspace(xmin, xmax, 101);
colors = parula(101);
v = interp1(xs, colors, x);

end

