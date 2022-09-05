#!/usr/bin/env python3
import sys
import numpy as np

##### Choose #####
# type = "x11"
type = "cairolatex"

no_deco = False

timelabel = True

xaxislabel = True
xticslabels = True
yaxislabel = True
yticslabels = True

colorbox = True
cbaxislabel = True
cbticslabels = True

plot_cross_sections = True
plot_cross_sections_number = 30
plot_cross_sections_every_dframe = 5 # 5, 10, 20

plot_curvature = True
plot_normal_velocity = True
plot_velocity_vs_curvature = True
set_curvature_as_one_over_r = False

frame = 70

###### initialisation ######

# Add a convenience function 'writeb' to 'stdout' to write strings as encoded binaries with out.writeb("string") rather than with out.write("string".encode()).
if sys.version_info >= (3,):
    out = sys.stdout.buffer # alias
    
    def writeb(self, s):
        import codecs
        self.write(codecs.utf_8_encode(s)[0]) # codecs.latin_1_encode or codecs.utf_8_encode
    from types import MethodType
    out.writeb = MethodType( writeb, out)
else:
    sys.exit(0, "Please use python 3")

# read metadata:

metadata = open("datadir/metadata.dat", "r").read().replace(r'"""', '"').replace("false", "0").replace("true", "1")
exec(metadata)
out.writeb(metadata) # includes setting parameters to gnuplot


##### min/max and multiplot layout parameters #####
out.writeb(r"""
x_min=-dx_coarse/2.
x_max=Lx+dx_coarse/2.
y_min=-dy_coarse/2.
y_max=Ly+dx_coarse/2.

u_min = 0
u_max = 1
u_incr = 0.1

curv_min = -0.15/dx
curv_max = -curv_min
curv_zero = 0

v_min = 0 # Allen-Cahn: [-0.35, 0.15], otherwise [0, 1] or [0, 10] etc.
v_max = 1

# # outward motion:
# curv_min = 0
# curv_max = R
# v_min = 0
# v_max = 0.4


# multiplot layout: 5 plots, 3 in top row (p3, p4, p5), 2 in bottom row (p1, p2).

lm = 23. # [mm] space between plots and canvas edge. Left and bottom larger for labels
rm = 25.
bm = 16.
tm = 10.

# width and height of all plots:
h1 = 70.
# w1 = 1.5*h1 # golden aspect ratio
w1 = 1.618*h1 # golden aspect ratio

h2 = h1
w2 = w1

h3 = 70.
w3 = h3 # square aspect ratio

h4 = h3
w4 = w3

h5 = h3
w5 = w3

# spacings between plots
h34 = 30. # horizontal spacing p3 <-> p4 and p4 <-> p5 in top row
h12 = w3 + w4 + w5 + 2*h34 - w1 - w2 # Adjust spacing p1 <-> p2 in bottom row to same total length as top row. We want w1 + w2 + h12 = w3 + w4 + w5 + 2*h34
v12 = 30. # vertical spacing row 1 <-> row 2

# bottom left corner of all plots:
l1 = lm
b1 = bm

l2 = l1 + w1 + h12
b2 = b1

l3 = l1
b3 = b1 + h1 + v12

l4 = l3 + w3 + h34
b4 = b3

l5 = l4 + w3 + h34
b5 = b3

W = l1+w1+h12+w2+rm
H = b1+h1+v12+h3+tm

########################################################################################
""")


if (plot_curvature or plot_velocity_vs_curvature) and (not calculate_curvature):
    import sys
    sys.stderr.write("Curvature was not calculated. Can't plot curvature")
    plot_curvature = False
    plot_velocity_vs_curvature = False

if (plot_normal_velocity or plot_velocity_vs_curvature) and (not calculate_normal_velocity):
    import sys
    sys.stderr.write("Normal velocity was not calculated. Can't plot normal velocity")
    plot_normal_velocity = False
    plot_velocity_vs_curvature = False

def reset_axes_decorations():
    global no_deco, timelabel, xaxislabel, xticslabels, yaxislabel, yticslabels, colorbox,cbaxislabel, cbticslabels
    if no_deco:
         timelabel = False
         xaxislabel = False
         xticslabels = False
         yaxislabel = False
         yticslabels = False
         colorbox = False
         cbaxislabel = False
         cbticslabels = False        
    if not timelabel:
        out.writeb(r"""
        # unset label 1
        """)        
    if not xaxislabel:
        out.writeb(r"""
        unset xlabel
        """)
    if not yaxislabel:
        out.writeb(r"""
        unset ylabel
        """)
    if not colorbox:
        out.writeb(r"""
        unset colorbox
        """)
    if not cbaxislabel:
        out.writeb(r"""
        unset cblabel
        """)
    if not xticslabels:
        out.writeb(r"""
        set format x ''
        """)
    elif type=="cairolatex":
        out.writeb(r"""
        set format xy '\footnotesize $%g$'
        """)
    if not yticslabels:
        out.writeb(r"""
        set format y ''
        """)
    if not cbticslabels:
        out.writeb(r"""
        set format cb ''
        """)
    elif type == "cairolatex":
        out.writeb(r"""
        set format cb '\footnotesize $%g$'
        """)

# When making several plots, we need to redefine axes etc. We do this with these functions
def reset_plot_u(): # p1
    out.writeb(r"""
set lmargin at screen l3/W
set rmargin at screen (l3+w3)/W
set bmargin at screen b3/H
set tmargin at screen (b3+h3)/H

set xrange [x_min : x_max]
set yrange [y_min : y_max]

set xlabel '\small $x$' offset 0,0.6 # '\small $x$ axis [$\mm$]' offset 0,0.3
set ylabel '\small $y$' offset 2.2,0 # epslatex: 3,0
set title '\small density $u$' offset 0,-0.3
set xtics out offset 0,0.3
set ytics out offset 0.6,0
    
set palette model RGB functions 1.1*gray**0.25, gray**0.75, 0 # black to gold *******BEST for color
# set palette defined (0 "black", 1 "#fefefe")
set cblabel ''
# set cblabel '\small $u$' offset -5,6.5 rotate by -0
set cbrange [u_min : u_max]
set zrange [u_min : u_max]
""")

def reset_plot_curvature(): # p2
    out.writeb(r"""
set lmargin at screen l4/W
set rmargin at screen (l4+w4)/W
set bmargin at screen b4/H
set tmargin at screen (b4+h4)/H

set xrange [x_min : x_max]
set yrange [y_min : y_max]

set xlabel '\small $x$'
set ylabel ''

set format y '' # remove tic labels
    
set title '\small curvature $\kappa$'
 
set cblabel ''
# set cblabel '\small $\kappa$'
set cbrange [curv_min : curv_max]
set zrange [curv_min : curv_max]

# unset colorbox
""")

def reset_plot_normal_velocity():
    out.writeb(r"""
set lmargin at screen l5/W
set rmargin at screen (l5+w5)/W
set bmargin at screen b5/H
set tmargin at screen (b5+h5)/H

set xrange [x_min : x_max]
set yrange [y_min : y_max]

set xlabel '\small $x$'
set ylabel ''
set title '\small velocity $v$'

set format y '' # remove tic labels

set palette model RGB functions 1.1*gray**0.25, gray**0.75, 0 # black to gold *******BEST for color

set cblabel ''
# set cblabel '\small $v$'
# set cbrange [*:*]
# set zrange [*:*]
set cbrange [v_min : v_max]
set zrange [v_min : v_max]
""")    


def reset_plot_cross_sections():
    out.writeb(r"""
set lmargin at screen l1/W
set rmargin at screen (l1+w1)/W
set bmargin at screen b1/H
set tmargin at screen (b1+h1)/H

set xrange [-0.05:1.05*Lx]
set yrange [u_min-0.05:u_max+0.05]

set xlabel '\small $x$'
set ylabel '\small $u(x,L_y/2,t)$'
set title '\small density profiles'

""")
    
def reset_plot_velocity_vs_curvature():
    out.writeb(r"""
set lmargin at screen l2/W
set rmargin at screen (l2+w2)/W
set bmargin at screen b2/H
set tmargin at screen (b2+h2)/H

# set xrange [-15:-0.5] # Allen-Cahn
set yrange [v_min:v_max]
set xrange [-60:20]
# set yrange [0:10]
# set xrange [*:*]
# set yrange [*:*]
# set xrange [-60:20]
# set yrange [0:0.2]
# set xrange [-60:20]
# set yrange [0:10]
# set xrange [-60:20]
# set yrange [0:1]
# set xrange [-250:20]
# set yrange [0:3]
# set yrange [0:10]
# set xrange [-0.1:1] # outward, D=0.005
# set yrange [0.10:0.15] # outward, D=0.005
# set xrange [0.1:1] # outward, D=0.02
# set yrange [0:0.3] # outward, D=0.02
# set xrange [curv_min:curv_max]
# set yrange [v_min:v_max]
    
set xlabel '\small $\kappa$'
set ylabel '\small $v$' offset 3,0
set title ''
# set title '\small growth law $v(\kappa)$'

set format xy '\footnotesize $%g$'

set palette model RGB functions 1.1*gray**0.25, gray**0.75, 0 # black to gold *******BEST for color
# set palette defined (0 "black", 1 "#fefefe")
# set cblabel ''
set cblabel '\small $u$' offset -7.1,7 rotate by -0
set cbrange [u_min : u_max]
set zrange [u_min : u_max]    

""")


##### gray out the region outside of the domain #####
reset_plot_u()
outside = open("datadir/outside.dat", "rb").read()

# contour at boundary of the domain
out.writeb(r"""
set contour base
set view map
set cntrparam levels discrete 1.0
set cntrlabel onecolor
unset surface
set format '%g'# reset otherwise with latex, points are outputted wrapped in latex commands!
set table $boundary
splot \
    '-' binary array=(""" + str(Nx_coarse) + "," + str(Ny_coarse) + r""") format='%float32' transpose origin=(x_min+0.5*dx_coarse,y_min+0.5*dy_coarse,0) dx=dx_coarse dy=dy_coarse notitle
""")
out.write(outside)
out.writeb(r"""
unset table
""")

# plot substrate as a separate png
out.writeb(r"""
set term pngcairo transparent size Nx_coarse, Ny_coarse

set margins 0,0,0,0

set border 0
unset xtics
unset ytics    
unset xlabel
unset ylabel
unset colorbox

set palette defined (0 "white", 1 "#eeeeee")

set output 'outside.png'

# set xrange [-1:1]
# set yrange [-1:1]

set xrange [*:*]
set yrange [*:*]

# see https://stackoverflow.com/questions/25321213/transparency-for-specific-values-in-matrix-using-gnuplot-while-preserving-the-pa
plot \
    '-' binary array=(""" + str(Nx_coarse) + "," + str(Ny_coarse) + r""") transpose format='%float32' using ($1 == 0 ? NaN : $1) with image notitle,\
    $boundary using ($1/dx_coarse):($2/dy_coarse) with lines lt 1 lw 1.5 lc rgb 'black' notitle
""")

# ----- Output initial data -----
out.write(outside)
out.writeb(r"""
unset output
""")

    
################



out.writeb(r"""
# Define a more sophisticated palette for the curvature map: red for negative values, green for zero, blue for positive values. We define three functions Red(gray), Green(gray), and Blue(gray), where gray is between 0 and 1 and corresponds to the value to plot via a linear conversion, see curv(gray). Here we define three functions Red(curv), Green(curv), Blue(curv) such that if curv < 0, we interpolate from min_RGB to zero_RGB, and if curv > 0, we interpolate from zero_RGB to max_RGB, with a nonlinearity provided by gamma. Note that for gnuplot, word(min_RGB, n) takes the nth word in the string min_RGB.

# Definition of colors to interpolate:
RGB_min=" 1.0 0.1 0.0 " # 'red'
RGB_zero=" .81640625 1.0 0.7578125 " # 'flash green' 208, 255, 194
RGB_max=" 0. .28 1. " # 'blue'        
RGB_min_sat=" 100./255 0.0 0.0 " # 'dark red'
RGB_max_sat=" 0.0 0.0 100./255 " # 'dark blue'

# Definition of nonlinearity (gamma correction)
gamma_min2zero = 1
gamma_zero2max = 1./gamma_min2zero

# Convenience:
Red_min = word(RGB_min, 1) + 0. # add 0. to convert to float
Red_zero = word(RGB_zero, 1) + 0.
Red_max = word(RGB_max, 1) + 0.
Red_min_sat = word(RGB_min_sat, 1) + 0.
Red_max_sat = word(RGB_max_sat, 1) + 0.
Green_min = word(RGB_min, 2) + 0. # add 0. to convert to float
Green_zero = word(RGB_zero, 2) + 0.
Green_max = word(RGB_max, 2) + 0.
Green_min_sat = word(RGB_min_sat, 2) + 0.
Green_max_sat = word(RGB_max_sat, 2) + 0.
Blue_min = word(RGB_min, 3) + 0. # add 0. to convert to float
Blue_zero = word(RGB_zero, 3) + 0.
Blue_max = word(RGB_max, 3) + 0.
Blue_min_sat = word(RGB_min_sat, 3) + 0.
Blue_max_sat = word(RGB_max_sat, 3) + 0.

# Helper functions for curvature
gray(curv) = (curv-curv_min)/(curv_max-curv_min)
curv(gray) = curv_min + gray*(curv_max - curv_min) # gnuplot assumes gray is between 0 and 1.

# Helper functions for curvature
gray(curv) = (curv-curv_min)/(curv_max-curv_min)
curv(gray) = curv_min + gray*(curv_max - curv_min) # gnuplot assumes gray is between 0 and 1.

# Helper function for palette
rgb(val, val_min, val_max, rgb_min, rgb_max, gamma) = rgb_min + (rgb_max - rgb_min)*((val-val_min)/(val_max-val_min))**(1./gamma)

# Functions to use when defining the palette: Red(gray), Green(gray), and Blue(gray)
Red(gray) = curv(gray) < 0 ? (curv(gray) < curv_min? Red_min_sat : rgb(curv(gray), curv_min, curv_zero, Red_min, Red_zero, gamma_min2zero)) : (curv(gray) > curv_max? Red_max_sat : rgb(curv(gray), curv_zero, curv_max, Red_zero, Red_max, gamma_zero2max))

Green(gray) = curv(gray) < 0 ? (curv(gray) < curv_min? Green_min_sat : rgb(curv(gray), curv_min, curv_zero, Green_min, Green_zero, gamma_min2zero)) : (curv(gray) > curv_max? Green_max_sat : rgb(curv(gray), curv_zero, curv_max, Green_zero, Green_max, gamma_zero2max))

Blue(gray) = curv(gray) < 0 ? (curv(gray) < curv_min? Blue_min_sat : rgb(curv(gray), curv_min, curv_zero, Blue_min, Blue_zero, gamma_min2zero)) : (curv(gray) > curv_max? Blue_max_sat : rgb(curv(gray), curv_zero, curv_max, Blue_zero, Blue_max, gamma_zero2max))

set xrange [x_min : x_max]
set yrange [y_min : y_max]

# set palette rgbformulae 22,13,10 # matlab-like
# set palette gray
set palette defined (0 "black", 1 "#fefefe")

set cbrange [u_min : u_max]
set zrange [u_min : u_max]

set tics out
set cbtics in

# set size ratio -1

k1 = 2.4048 # first zero of the J0 Bessel function
set key spacing 2
""")

reset_axes_decorations()



###### plotting ######
out.writeb("# FRAME %d" % frame) # not necessary, but helps to debug
u = open("datadir/u_%d.dat" % frame, "rb").read()
if plot_curvature and calculate_curvature:
    kappa = open("datadir/curvature_%d.dat" % frame, "rb").read()        
if plot_normal_velocity and calculate_normal_velocity:
    v = open("datadir/v_%d.dat" % frame, "rb").read()
        
        
# ===== contour lines =====
# preparation: calculate contour levels and save into gnuplot data table $contourdata (by calling splot) for later plotting

reset_plot_u()
# contour at u_c
out.writeb(r"""
set contour base
set view map
set cntrparam levels discrete """ + "%f" % u_c + r"""
set cntrlabel onecolor
unset surface
set format '%g'# reset otherwise with latex, points are outputted wrapped in latex commands!
set table $zerocontour
# set table "zerocontour-".sprintf("%.3d",""" + "%d" % frame + r""").".dat" # for debugging
splot \
    '-' binary array=(""" + str(Nx_coarse) + "," + str(Ny_coarse) + r""") format='%float32' transpose origin=(x_min+0.5*dx_coarse,y_min+0.5*dy_coarse,0) dx=dx_coarse dy=dy_coarse every 1:1
""")# offset by half a grid size in 'origin' unless wanting to show cell-centered values
out.write(u)

# contour at u_initial
out.writeb(r"""
set cntrparam levels discrete u_initial
set table $boundary
# set table "boundary-".sprintf("%.3d",""" + "%d" % frame + r""").".dat" # for debugging
splot \
    '-' binary array=(""" + str(Nx_coarse) + "," + str(Ny_coarse) + r""") format='%float32' transpose origin=(x_min+0.5*dx_coarse,y_min+0.5*dy_coarse,0) dx=dx_coarse dy=dy_coarse every 1:1
""")
out.write(u)
out.writeb(r"""
unset table
""")
    
# plot other contours:
out.writeb(r"""
set cntrparam levels incremental u_min, u_incr, u_max
set table $contourdata
# set table "contourdata-".sprintf("%.3d",""" + "%d" % frame + r""").".dat" # for debugging
splot \
    '-' binary array=(""" + str(Nx_coarse) + "," + str(Ny_coarse) + r""") format='%float32' transpose origin=(x_min+0.5*dx_coarse,y_min+0.5*dy_coarse,0) dx=dx_coarse dy=dy_coarse every 1:1
""")
out.write(u)
out.writeb(r"""
unset table
""")


# ===== Plot cell density u + contours into png =====
# plot using pngcairo first to rasterise the data, then use the created png data with epslatex or pdflatex
reset_axes_decorations()
reset_plot_u()

out.writeb(r"""
# set term pngcairo transparent size w3/10. cm, h3/10. cm
# set term pngcairo transparent size (Nx_coarse-1)/2, (Ny_coarse-1)/2
set term pngcairo transparent size Nx_coarse, Ny_coarse

set margins 0,0,0,0
    
set border 0
unset xtics
unset ytics    
unset xlabel
unset ylabel
unset colorbox
    
set output 'udata.png'
""")
    
out.writeb(r"""
plot \
    '-' binary array=(""" + str(Nx_coarse) + "," + str(Ny_coarse) + r""") format='%float32' transpose origin=(0.5*dx_coarse, 0.5*dy_coarse) dx=dx_coarse dy=dy_coarse every 1:1 title '' with image,\
    $contourdata with lines lt 1 lw 1 lc rgb 'white' dt 2 title '',\
    $boundary with lines lt 1 lw 1.5 lc rgb 'black' title '',\
    $zerocontour with lines lt 1 lw 1 lc rgb 'red' title ''
""")# offset by half a grid size in 'origin' unless wanting to show cell-centered values


# ----- Output all data in order -----
out.write(u)    
out.writeb(r"""
unset output
""")

    
# ===== Plot curvature map into png =====
reset_plot_curvature()
out.writeb(r"""
# set term pngcairo transparent size w3/10. cm, h3/10. cm
# set term pngcairo transparent size (Nx_coarse-1)/2, (Ny_coarse-1)/2
set term pngcairo transparent size Nx_coarse, Ny_coarse

set margins 0,0,0,0
    
set border 0
unset xtics
unset ytics    
unset xlabel
unset ylabel
unset colorbox
    
set output 'kappadata.png'
""")
    
out.writeb(r"""
plot \
    '-' binary array=(""" + str(Nx_coarse) + "," + str(Ny_coarse) + r""") format='%float32' transpose origin=(x_min+0.5*dx_coarse,y_min+0.5*dy_coarse) dx=dx_coarse dy=dy_coarse every 1:1 using (255*Red(gray($1))):(255*Green(gray($1))):(255*Blue(gray($1))):($1 == NaN? 0 : 255) title '' with rgbalpha
""")# offset by half a grid size in 'origin' unless wanting to show cell-centered values

# ----- Output all data in order -----    
out.write(kappa)

out.writeb(r"""
unset output
""")
            
# ===== Plot normal velocity map into png =====
reset_plot_normal_velocity()
out.writeb(r"""
# set term pngcairo transparent size w3/10. cm, h3/10. cm
# set term pngcairo transparent size (Nx_coarse-1)/2, (Ny_coarse-1)/2
set term pngcairo transparent size Nx_coarse, Ny_coarse

set margins 0,0,0,0
    
set border 0
unset xtics
unset ytics    
unset xlabel
unset ylabel
unset colorbox
    
set output 'vdata.png'
""")
out.writeb(r"""
plot \
    '-' binary array=(""" + str(Nx_coarse) + "," + str(Ny_coarse) + r""") format='%float32' transpose origin=(x_min+0.5*dx_coarse,y_min+0.5*dy_coarse) dx=dx_coarse dy=dy_coarse every 1:1 title '' with image
""")# offset by half a grid size in 'origin' unless wanting to show cell-centered values

# ----- Output all data in order -----    
out.write(v)
out.writeb(r"""
unset output
""")

# ====== Plot curvature vs velocity points into png =====
reset_plot_velocity_vs_curvature()
out.writeb(r"""
resolution=1
set term pngcairo transparent size resolution*w2/10. cm, resolution*h2/10. cm

set margins 0,0,0,0
    
set border 0
unset xtics
unset ytics    
unset xlabel
unset ylabel
unset colorbox
    
set output 'vkudata.png'
""")

u_np = np.fromfile("datadir/u_%d.dat" % frame, dtype=np.float32, count=Nx_coarse*Ny_coarse).reshape(Nx_coarse,Ny_coarse)

v_np = np.fromfile("datadir/v_%d.dat" % frame, dtype=np.float32, count=Nx_coarse*Ny_coarse).reshape(Nx_coarse,Ny_coarse)

if set_curvature_as_one_over_r: # better precision for radially symmetric solutions
    kappa_np = np.fromfunction(lambda i,j: -1./np.sqrt(((i-Nx_coarse//2)*dx_coarse)**2 + ((j-Ny_coarse//2)*dx_coarse)**2), (Nx_coarse, Ny_coarse), dtype=np.float32)
else:
    kappa_np = np.fromfile("datadir/curvature_%d.dat" % frame, dtype=np.float32, count=Nx_coarse*Ny_coarse).reshape(Nx_coarse,Ny_coarse)



# only include points satisfying some criteria. 
def include_point(i,j):
    # uncomment the inclusion criterion that applies

    # # inclusion criterion: subdiagaonal in 1st quadrant
    # return i*dx_coarse>=Lx/2. and j*dy_coarse>=Ly/2. and j<=i*Ly/Lx

    # inclusion criterion: points in 1st quadrant below the 60 degrees line
    import math
    return i*dx_coarse>=Lx/2. and j*dx_coarse>=Lx/2. and j<=i*math.sqrt(3)*Ly/Lx

    # # all included
    # return True

# 1) count how many, to allocate the right array size
N_elem = 0
for i in range(Nx_coarse):
    for j in range(Ny_coarse):
        if include_point(i,j):
            N_elem = N_elem + 1

# 2) allocate array and fill in with data
v_k_u = np.zeros((N_elem, 3), dtype=np.float32)
N_elem = 0
for i in range(Nx_coarse):
    for j in range(Ny_coarse):
        if include_point(i,j):
            v_k_u[N_elem, :] = [kappa_np[i,j], v_np[i,j], u_np[i,j]]
            N_elem = N_elem + 1

out.writeb(r"""
plot \
    '-' binary record=""" + str(N_elem) + r""" format='%float32%float32%float32' using 1:2:3 title '' with points pt 1 lc palette
""")

out.write(v_k_u.tobytes())
out.writeb(r"""
unset output
""")



################################# multiplot #################################
if type == "x11":
    out.writeb(r"""
set term x11 size 2*W, 2*H #font "default" 10 #"LMRoman10-Regular" 10
""")

elif type == "wxt":
    out.writeb(r"""
set term wxt persist size 2*W, 2*H #font "default" 10 #"LMRoman10-Regular" 10
""")

elif type == "epslatex":
    out.writeb(r"""
set term epslatex standalone color level3 size 0.5*W/10. cm, 0.5*H/10. cm font "default, 10" colortext header '\usepackage{amsmath}\usepackage{amssymb}\usepackage[utopia,greekuppercase=italicized]{mathdesign}\usepackage{newtxtext}\edef\partial{\mathchar\number\partial\noexpand\!}\usepackage[T1]{fontenc}\usepackage[expproduct=cdot]{siunitx}\newcommand{\um}{\ensuremath{\micro\metre}}'

set format xy '\tiny %g' # latex command may offset the ylabel a lot. Use offset in set ylabel.
""")

elif type == "pdfcairo":
    out.writeb(r"""
set term pdfcairo enhanced size 0.5*W/10. cm, 0.5*H/10. cm
""")

elif type == "cairolatex":
    out.writeb(r"""
set term cairolatex pdf standalone color size 0.5*W/10. cm, 0.5*H/10. cm colortext font "default, 10" header '\usepackage{amsmath}\usepackage{amssymb}\usepackage[utopia,greekuppercase=italicized]{mathdesign}\usepackage{newtxtext}\edef\partial{\mathchar\number\partial\noexpand\!}\usepackage[T1]{fontenc}\usepackage[expproduct=cdot]{siunitx}\newcommand{\um}{\ensuremath{\micro\metre}}'

set format xy '\tiny $%g$' # latex command may offset the ylabel a lot. Use offset in set ylabel.
""")
    
# ===== output =====
if type == "epslatex" or type == "cairolatex":
    out.writeb(r"""
set output "frame-".sprintf("%.4d",""" + "%d" % frame + r""").".tex"
""")
        
elif type == "pdfcairo":
    out.writeb(r"""
set output "frame-".sprintf("%.4d",""" + "%d" % frame + r""").".pdf"
""")

# ===== Start multiplot and time label =====
out.writeb(r"""
set multiplot title sprintf('\footnotesize$t=%5.2f$', """ + "%d" % frame + r"""*dt_frame) offset screen -0.46,0.015

set colorbox
set border
""")


# ===== Plot density u from png =====
reset_plot_u()
out.writeb(r"""
set format xy '\footnotesize $%g$'
    
plot 'udata.png' binary filetype=png with rgbalpha axes x2y2 title '',\
    'outside.png' binary filetype=png with rgbalpha axes x2y2 notitle,\
    1/0 with lines lc palette title '' # force showing colorbox
""")
    
# ===== Plot curvature map from png =====
reset_plot_curvature()
out.writeb(r"""
plot 'kappadata.png' binary filetype=png with rgbalpha axes x2y2 title '',\
    1/0 with lines lc palette title '' # force showing colorbox
""")
    
# ===== Plot normal velocity map from png =====
reset_plot_normal_velocity()
out.writeb(r"""
plot 'vdata.png' binary filetype=png with rgbalpha axes x2y2 title '',\
    1/0 with lines lc palette title '' # force showing colorbox
""")
    
# ====== Plot curvature vs velocity points from png =====
reset_plot_velocity_vs_curvature()
out.writeb(r"""
plot 'vkudata.png' binary filetype=png with rgbalpha axes x2y2 title '',\
    1/0 with lines lc palette title ''
""")
        
    
# ====== Plot cross-sections =====
reset_plot_cross_sections()

out.writeb(r"""
set arrow 1 from 0.0,u_c to 1.0*Lx,u_c nohead dt 3 lw 3 lc rgb 'red' front
""")
    
out.writeb(r"""
plot \
""")

# frame_end = Nt_frame
frame_end = frame
for f in range(0,frame_end,plot_cross_sections_every_dframe):
    out.writeb(r"""'-' binary record=""" + "%d" % Nx_coarse + r""" format='%float32' using ($0*dx_coarse):1 with lines lw 1 lc rgb 'black' title '',\
""")
out.writeb(r"""'-' binary record=""" + "%d" % Nx_coarse + r""" format='%float32' using ($0*dx_coarse):1 with lines lw 2 lc rgb 'red' title '',\
""")    
out.writeb(r"""
""")
for f in range(0,frame_end,plot_cross_sections_every_dframe):
    u = np.fromfile('datadir/u_%d.dat' %  f, dtype=np.float32)
    u.shape = (Nx_coarse,Ny_coarse)
    out.write(u[:,Ny_coarse//2].tobytes())

u = np.fromfile('datadir/u_%d.dat' % frame, dtype=np.float32)
u.shape = (Nx_coarse,Ny_coarse)
out.write(u[:,Ny_coarse//2].tobytes())

out.writeb(r"""
unset arrow 1
#refresh
    """)
    

# ===== Frame end (unnecessary, but helps to debug) =====
out.writeb(r"""
unset multiplot
unset output
""")


