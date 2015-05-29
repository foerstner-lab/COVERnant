import sys
sys.path.append(".")
import unittest
from covernantlib.extract import CoverageExtractor


class TestCoverageExtractor(unittest.TestCase):

    def setUp(self):
        self.coverage_extractor = CoverageExtractor()

if __name__ == "__main__":
        unittest.main()
