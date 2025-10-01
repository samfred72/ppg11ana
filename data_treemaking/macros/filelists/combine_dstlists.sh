paste <(ls dst_calo_*) <(ls dst_calofitting*) | while read file1 file2; do
  length1=`ls $file1 | wc -l`
  length2=`ls $file2 | wc -l`
  if [ $length1 -ne $length2 ]; then
    continue
  fi
  paste <(cat $file1) <(cat $file2) | while read dst1 dst2; do 
    echo "$dst1" "$dst2" >> combined_$file1 
  done
done
