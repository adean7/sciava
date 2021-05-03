
current_dir=`pwd`

cd /Users/adam/Documents/projects/sciava/

python setup.py install >> /dev/null 2>&1

cd $current_dir

python test.py

