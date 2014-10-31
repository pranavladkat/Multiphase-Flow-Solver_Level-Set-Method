#!/bin/bash

cd Re5Bo1

for f in `ls *.eps`; do
     convert -density 200 $f -flatten ${f%.*}.png;
done

mencoder "mf://*.png" -mf fps=10 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4

rm *.eps
#rm *.png
