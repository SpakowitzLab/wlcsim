#/bin/bash



# compress files with python
temp="$1/binaryData/r"
for file in "$temp"*.npz
do
    base=`basename "$file"`
    newbase=`echo "$base" | cut -d'.' -f1`
    python dataTools/expand_r.py $file > "$1/$newbase";
    rm "$file"
done
