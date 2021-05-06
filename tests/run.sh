
current_dir=$(pwd) ;

cd /Users/adam/Documents/projects/sciava/ || exit ;

python setup.py install >> /dev/null 2>&1 ;

cd "$current_dir" || exit ;

python test.py ;
