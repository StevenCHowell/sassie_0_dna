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

SUFFIX=_cap.pdb
for file in "${files[@]}" 
do
    if [[ $file == *".pdb" ]]
    then
        CAPPDB=${file%.pdb}$SUFFIX
        if [ ! -f $CAPPDB.dat ]
        then
            echo
            echo creating $CAPPDB then running:   foxs -q $QMAX -s $PROFILESIZE $CAPPDB 
            eval \awk "'/ P / || / CA /'" $file > $CAPPDB
            eval \foxs -q $QMAX -s $PROFILESIZE $CAPPDB
            echo
        else
            echo $CAPPDB.dat already calculated
        fi
    fi
done


