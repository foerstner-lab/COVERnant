import unittest

all_tests = unittest.TestLoader().discover("./tests")
unittest.TextTestRunner(verbosity=2).run(all_tests)
