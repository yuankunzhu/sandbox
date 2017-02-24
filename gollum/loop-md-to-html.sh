FILES="*.md"
for f in $FILES
do
    base=`basename $f ".md"`
    pandoc -s -f markdown $f > "$base".html
done