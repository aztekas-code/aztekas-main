#!/usr/bin/gnuplot -c

#####################################
# Filame as argument
# You can add your own arguments from 
# ARG3 TO ARG9
file_name=ARG1 
output_name=ARG2

#Do not erase these lines
len=strlen(output_name)
out_name=substr(output_name,1,len-4)
ext=substr(output_name,len-2,len)

if ( file_name eq '' || output_name eq '' ){
###################################################
# This file is run in your terminal as:
# $ gnuplot -c plot_1D.gp file.dat
print "Gnuplot script for simple 1D aztekas graph"
print ""
print "Execute this file in your terminal as:"
print "  $ ./plot_1D.gnu file.dat output.ext"
print "where"
print "  file.dat: your aztekas output"
print "  output.ext: name of the output file with"
print "              'ext' the extension (png,pdf)"
print ""

set out  

}else{
##################################################
# Writting info on the screen
# $ gnuplot -c plot_1D.gp file.dat
print "Gnuplot script for simple 1D aztekas graph"
print "File name to plot: ", file_name
print "Output name of the plot: ", output_name

#######################################
# This gnuplot graph is for 1D aztekas
# ASCII output in which the columns are
# that are needad in @ at line:
# > plot filename u @:@
# are
# 1 : X1
# 2 : Density
# 3 : Pressure
# 4 : Velocity in X1
#
# This file can also be used for 2D 
# aztekas ASCII output, if you consider
# 1 : X1
# 2 : X2
# 3 : Density
# 4 : Pressure
# 5 : Velocity in X1
# 6 : Velocity in X2
#######################################

if ( ext eq 'pdf' ){

##############################################
## If you want to export the graph using LaTex
## uncomment the folowing lines and comment
## the png ones.
## To generate the plot run
## $ ./gtex2pdf plot_1D
## that generate plot_1D_out.pdf, or
## $ ./gtex2eps plot_1D
## that generate plot_1D_out.eps 
set terminal epslatex colour input "" 12
set output out_name.".tex"

##################################
# Modify labels, titles, ticks,
# scale, range... whatever you
# want (Check $ man gnuplot)
# For LaTeX output, all math that
# you want is written between $ $
# and symbols with double \\, e.g.
# set ylabel "$\\rho$"
set xlabel "$x$"
set ylabel "$u$"
}

if ( ext eq 'png' ){
################################
## Standard output of a png file
set terminal png
set output out_name.".png"

##################################
# Modify labels, titles, ticks,
# scale, range... whatever you
# want (Check $ man gnuplot)
set xlabel "x"
set ylabel "u"
}

########################################
# Plot your file with all the parameters
# legends and style you want 
# (Check $ man gnuplot)
plot file_name u 1:2 notitle

set out

if ( ext eq 'pdf' ){
   system("$AZTEKAS_PATH/Graphics/gtex2pdf ".out_name)
}

}
