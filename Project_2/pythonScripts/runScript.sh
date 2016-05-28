FILES="../dataFiles/smaug/localenergiesN6*"
for f in $FILES
do
	python blocking.py $f 30000
done
