from pygromos.analysis.error_estimate import error_estimator
import numpy as np
import unittest


class test_ee(unittest.TestCase):
    error_estimate_class = error_estimator
    test_array = np.arange(10000)

    def test_constructor(self):
        print(self.error_estimate_class(self.test_array))

    def test_error_estimate(self):
        # Create array
        error_estimate = self.error_estimate_class(self.test_array).calculate_error_estimate()
        np.testing.assert_almost_equal(desired=866.5903464697677, actual=error_estimate, decimal=2)
