#!/bin/bash
declare -a num
declare -a message

for file in $(find "$1" -name "*.log")
do

    DoesItMeetInitial="false"
    IsItMessage="false"

    i=0
    j=0

    while read -r line
    do
        current=$line
        if [[ $current == *"== inital =="* ]]; then
            DoesItMeetInitial="true"
        fi

        if [ "$DoesItMeetInitial" == "true" ] ; then
            if [[ $current == *"======="* ]]; then
	        message[i]=$current
                ((i++))
            elif [[ "$current" == *"Number of"* ]];then
                temp_num=`echo $current | cut -d":" -f2 | tr -d " "`
                num[j]=`expr ${num[j]} + $temp_num`
                ((j++))
	        temp_message=`echo $current | cut -d":" -f1`
                message[i]=$temp_message
                ((i++))
            fi
        fi

    done < $file

    if [ "$DoesItMeetIntial" == "false" ]; then
        echo "cannot find initial!"
	exit 1
    fi

done

k=0
if [ "$DoesItMeetInitial" == "true" ]; then
    ELEMENTS=${#message[@]}
    for ((i = 0; i < $ELEMENTS; i++)); do
        if [[ "${message[${i}]}" == *"===="* ]]; then
            echo ${message[${i}]}
    elif [[ "${message[${i}]}" == *"Number of"* ]]; then
	    echo ${message[${i}]}": "${num[k]}
	    ((k++))
    fi
    done
fi

