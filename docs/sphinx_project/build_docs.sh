#clean up:
make clean;

#make doku
##configurations
sphinx-apidoc -o _source ../../pygromos

cp ../../examples/*ipynb ./Examples

python conf.py

##execute making docu
make html
make latex

ln ../index.html /_build/html/index.html