svg="$1.svg"
#echo $svg

pdf="$1.pdf"
#echo $pdf

if [ "$svg" -nt "$pdf" ]; then
    /Applications/Inkscape.app/Contents/Resources/bin/inkscape -z -D --file="$svg" --export-pdf="$pdf" --export-latex --export-area-drawing
fi