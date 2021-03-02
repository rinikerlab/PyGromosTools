import os
import importlib
import unittest


if __name__ == "__main__":
    root_dir = os.path.dirname(__file__)
    test_files = []

    #gather all test_files
    for dir in os.walk(root_dir):
        test_files.extend([dir[0]+"/"+path for path in dir[2] if( path.startswith("dev") and path.endswith(".py") and not "test_run_all_tests" in path)])

    #print("\n".join(test_files))
    #get module import paths
    modules = []
    for test_module in test_files:
        module_name = "pygromos" + test_module.replace(os.path.dirname(root_dir), "").replace("/", ".").replace(".py", "")
        modules.append(module_name)

    #LOAD TESTS
    print("LOAD TESTS")
    suite = unittest.TestSuite()
    test_loader = unittest.TestLoader()
    first = True

    for test_module in modules:
        imported_test_module = __import__(test_module, globals(), locals(), ['suite'])
        if(first):
            suite = test_loader.loadTestsFromModule(imported_test_module)
            first = False
        else:
            tmp = test_loader.loadTestsFromModule(imported_test_module)
            suite.addTest(tmp)


    print("RUN TESTS")
    try:
        print("TEST SUIT TESTS: ", suite.countTestCases())
        test_runner = unittest.TextTestRunner(verbosity=5)
        test_runner.run(suite)
        exit(0)
    except Exception as err:
        print("Test did not finish successfully!\n\t"+"\n\t".join(err.args))
        exit(1)
