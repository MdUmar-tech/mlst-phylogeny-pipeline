cd "/Users/mdumar/Desktop/md10/MLST/mlst_extracted_genes"

for d in */; do
    newname=$(echo "$d" | sed 's/Strain//g' | tr -s ' ')
    newname=$(echo "$newname" | sed 's/ *$//')  # Remove trailing space
    mv "$d" "$newname"
done
