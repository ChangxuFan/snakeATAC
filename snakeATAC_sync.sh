#!/bin/bash
if [ $1 == "ris" ]
then
	echo "ris"
	target=fanc@compute1-client-1.ris.wustl.edu:/home/fanc/software/snakeATAC/
else
	echo "htcf"
	target=fanc@htcf.wustl.edu:/home/fanc/software/snakeATAC/
fi

rsync -a --include-from="/bar/cfan/software/snakeATAC/include.txt" --exclude="*" \
--delete --delete-excluded \
~/software/snakeATAC/ $target