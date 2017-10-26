#/bin/bash


# Error if there are already npz files
#ls ~/dataTools/r*.npz 2>/dev/null  | wc -l > tempxxx
#read lines < tempxxx
#rm tempxxx
#if [ $lines -gt 0 ]
#then
#    echo "there are already $lines npz files"
#    exit
#fi


# compress files with python
for file in $1/r* 
do base=`basename "$file"`;
   dir=`dirname "$file"`;
   mkdir -p "$dir/binaryData"
   python dataTools/condense_r.py $file "$dir/binaryData/$base";
   rm "$file"
done
