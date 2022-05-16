cut -f 1 ../test/samples.lst | while read line; do echo "cd $line && make && cd -";done > out/1.standard/standard.sh
cd out/1.standard | sh standard.sh
cd out/1.standard && make && cd -
cd out/2.identification && make && cd -
cd out/3.filter && make && cd -
cd out/4.ds_score && make && cd -
