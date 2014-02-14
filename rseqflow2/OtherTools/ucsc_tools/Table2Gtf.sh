cat $1 |grep -v 'name' >temp1.txt
cut temp1.txt -f2- >temp2.txt
genePredToGtf file temp2.txt $2
rm temp1.txt -f
rm temp2.txt -f

