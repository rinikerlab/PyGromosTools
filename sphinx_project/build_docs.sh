#clean up:
make clean;

#make doku
##configurations
sphinx-apidoc -o _source ../../pygromos;

cp ../../examples/ex*ipynb ./Examples;
cp -r ../../examples/developer_examples ./Examples;
cp -r ../../examples/example_files ./Examples;
cp ../../examples/t*ipynb ./Tutorials;
cp -r ../../examples/example_files ./Tutorials;

python conf.py;

##execute making docu
make html;

cp -rf _build/html/* ../
