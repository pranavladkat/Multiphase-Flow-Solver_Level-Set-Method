#!/bin/bash

# cd to folder with images saved by screenshot function
cd Re5Bo1

# convert eps to png
for f in `ls *.eps`; do
     convert -density 200 $f -flatten ${f%.*}.png;
done

#build movie
mencoder "mf://*.png" -mf fps=10 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4

#remove images
rm *.eps
rm *.png
