#!/bin/bash

for i in int double ; do
    sed -e "s|@@TYPE@@|${i}|g" jkmarray_source.c >jkmarray_${i}.c
    indent jkmarray_${i}.c
    sed -e "s|@@TYPE@@|${i}|g" jkmarray_source.h >jkmarray_${i}.h
    indent jkmarray_${i}.h

    sed -e "s|@@TYPE@@|${i}|g" jmarray_source.c >jmarray_${i}.c
    indent jmarray_${i}.c
    sed -e "s|@@TYPE@@|${i}|g" jmarray_source.h >jmarray_${i}.h
    indent jmarray_${i}.h
done
