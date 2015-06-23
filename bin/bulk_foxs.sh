#!/bin/bash -f

QMAX=0.2
PROFILESIZE=50

while [[ $# > 0 ]]
do
key="$1"
case $key in
    -q|--max_q)
    QMAX="$2"
    shift # past argument
    ;;
    -s|--profile_size)
    PROFILESIZE="$2"
    shift # past argument
    ;;
    *) files+=("$1") ;;
esac
shift # past argument or value
done

for file in "${files[@]}" 
do
    if [[ $file == *".pdb" ]]
    then
        if [ ! -f $file.dat ]
            then
            echo
            # echo running: foxs -q $QMAX -s $PROFILESIZE $file 
            eval \foxs -q $QMAX -s $PROFILESIZE $file
            echo
        else
            echo $file.dat already calculated
        fi
    fi
done


