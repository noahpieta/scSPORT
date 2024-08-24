for dir in ~/DATA/project/scSHAPE_Map/RiboSwitches_scSHAPE/riboswitch_in_h9_transcriptome/done/R*
do
    s=${dir#*done/}
    echo $s

    mkdir $s
    cd $s
    cp $dir/*.mutrate.txt.gz .
    cd ../

done

