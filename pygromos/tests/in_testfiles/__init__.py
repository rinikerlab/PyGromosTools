import os
in_test_file_path = os.path.dirname("/home/runner/work/PyGromosTools/pygromos/tests/in_testfiles/__init__.py") #__file__)

print("Input Dir: ",in_test_file_path)
print("Input Dir ABS: ", os.path.abspath(in_test_file_path))
print("Im here: ", os.getcwd())
 
