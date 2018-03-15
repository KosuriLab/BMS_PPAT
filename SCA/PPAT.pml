delete all
load 1H1T.pdb, main
hide all
bg_color white
show cartoon, (chain E)
color white

set_color color1, [0.000,1.000,0.400]
create sector1, (resi 8,10,11,12,14,16,17,19,23,41,50,72,88,90,94,98,105,114,127,128,146) & (chain A)
color color1, sector1
show spheres, sector1
show surface, sector1

set_color color2, [1.000,0.000,0.000]
create sector2, (resi 4,9,20,36,38,54,67,76,119,126,129,132,134) & (chain A)
color color2, sector2
show spheres, sector2
show surface, sector2

set_color color3, [0.200,0.000,1.000]
create sector3, (resi 7,34,43,69,89,91,92,93,95,96,99,100,118,131) & (chain A)
color color3, sector3
show spheres, sector3
show surface, sector3

set_color color4, [1.000,0.900,0.000]
create sector4, (resi 44,45,49,102,116,123,133,135,138,149,156) & (chain A)
color color4, sector4
show spheres, sector4
show surface, sector4

set_color color5, [1.000,0.000,0.600]
create sector5, (resi 6,53,74,86,87,117,153) & (chain A)
color color5, sector5
show spheres, sector5
show surface, sector5

set_color color6, [0.000,1.000,1.000]
create sector6, (resi 15,28,71,73,75,97,101,104,113,125,140,144) & (chain A)
color color6, sector6
show spheres, sector6
show surface, sector6

set_color color_ic1, [1.000,0.000,0.000]
create ic_1, (resi 8,10,11,12,14,16,17,19,23,41,50,72,88,90,94,98,105,114,127,128,146) & (chain A)
color color_ic1, ic_1
show spheres, ic_1
show surface, ic_1

set_color color_ic2, [1.000,0.857,0.000]
create ic_2, (resi 4,9,20,36,38,54,67,76,119,126,129,132,134) & (chain A)
color color_ic2, ic_2
show spheres, ic_2
show surface, ic_2

set_color color_ic3, [0.286,1.000,0.000]
create ic_3, (resi 7,34,43,69,89,91,92,93,95,96,99,100,118,131) & (chain A)
color color_ic3, ic_3
show spheres, ic_3
show surface, ic_3

set_color color_ic4, [0.000,1.000,0.571]
create ic_4, (resi 44,45,49,102,116,123,133,135,138,149,156) & (chain A)
color color_ic4, ic_4
show spheres, ic_4
show surface, ic_4

set_color color_ic5, [0.000,0.571,1.000]
create ic_5, (resi 6,53,74,86,87,117,153) & (chain A)
color color_ic5, ic_5
show spheres, ic_5
show surface, ic_5

set_color color_ic6, [0.286,0.000,1.000]
create ic_6, (resi 15,28,71,73,75,97,101,104,113,125,140,144) & (chain A)
color color_ic6, ic_6
show spheres, ic_6
show surface, ic_6

set_color color_ic7, [1.000,0.000,0.857]
create ic_7, (resi 29,63) & (chain A)
color color_ic7, ic_7
show spheres, ic_7
show surface, ic_7

zoom
set transparency, 0.4
ray
png PPAT
