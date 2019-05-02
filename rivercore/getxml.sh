#/bin/sh -f

xmldir='/home/dai.56//data4/EarthDEM/alaska_2018oct26/qb_wv_alaska_metadata/';
dir='/fs/byo/howat-data5/pgc_deliv/chunli/sag/deliv3/ortho_imagery_non-max_ona/';
#find /*/ArcticDEM/region*/strips/2m/ -type f -name '*meta.txt' > Reglist
find $dir -type f -name '*.tif' > Reglist
for line in `cat Reglist`; do
#infile=$dir/$line
infile=$line
fbname0=$(basename $infile)
fbname=${fbname0/_u16ns3413.tif/.xml}
fbname2=${fbname0/_u16ns3413.tif/_u16ns3413.xml}
xmlfile=$xmldir/$fbname
echo $line  $xmlfile 
cp $xmlfile $dir/$fbname2
#rm $dir/$fbname
done

