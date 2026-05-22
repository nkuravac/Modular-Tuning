function G = color_gradient(rgbcolor1,rgbcolor2,n)
    G = [linspace(rgbcolor1(1),rgbcolor2(1),n)',linspace(rgbcolor1(2),rgbcolor2(2),n)',linspace(rgbcolor1(3),rgbcolor2(3),n)'];
end
