#!/bin/bash

if [ ! -d  "./bin" ]; then
    mkdir bin
fi

for file in `ls MMult*`
do
    NEW=${file%.*}
    echo $NEW
    make NEW=${file%.*}  
    mv  test_$NEW.x ./bin/
done
