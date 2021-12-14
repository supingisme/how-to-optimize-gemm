#!/bin/ksh

if [ ! -d  "./result" ]; then
    mkdir result
fi

for file in `ls test*`
do
    NEW=${file%.*}
    echo $NEW
    echo "version = '${NEW}';" > ./result/$NEW.m
    on -C7 ./${file}  >> ./result/$NEW.m
    cat ./result/$NEW.m
done
