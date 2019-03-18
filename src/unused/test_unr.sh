#!/bin/bash

number_a=`grep 'Number of alpha electrons' $1 | awk '{print $6}'`
number_b=`grep 'Number of beta electrons' $1 | awk '{print $6}'`

t=`echo "$((number_a - number_b))"`

if [ $t -eq 0 ] ; then
echo "restricted $number_a $number_b $t"
else
echo "unrestricted $number_a $number_b $t"
fi


