#clean up:
make clean;

#make doku
##configurations
sphinx-apidoc -o _source ../../pygromos

cp ../../examples/ex*ipynb ./Examples
cp ../../examples/t*ipynb ./Tutorials

python conf.py

##execute making docu
make html

cp -r _build/html/* ../
