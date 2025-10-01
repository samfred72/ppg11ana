rm combined_dstlist.list
touch combined_dstlist.list
paste <(ls dst_calo_cluster.list) <(ls dst_global.list) <(ls g4hits.list) <(ls dst_mbd_epd.list) | while read file1 file2 file3 file4; do
  length1=`ls $file1 | wc -l`
  length2=`ls $file2 | wc -l`
  length3=`ls $file3 | wc -l`
  length4=`ls $file4 | wc -l`
  if [ $length1 -ne $length2 ]; then
    continue
  fi
  if [ $length1 -ne $length3 ]; then
    continue
  fi
  if [ $length1 -ne $length4 ]; then
    continue
  fi
  if [ $length2 -ne $length3 ]; then
    continue
  fi
  if [ $length2 -ne $length4 ]; then
    continue
  fi
  if [ $length3 -ne $length4 ]; then
    continue
  fi
  paste <(cat $file1) <(cat $file2) <(cat $file3) <(cat $file4) | while read dst1 dst2 dst3 dst4; do 
    echo "$dst1 $dst2 $dst3 $dst4" >> combined_dstlist.list
  done
done
