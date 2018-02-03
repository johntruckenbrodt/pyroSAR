import unittest
import pyroSAR.ancillary as anc

class TestListFunctions(unittest.TestCase):

    def test_dissolve_with_lists(self):
        self.assertEqual(anc.dissolve([[1, 2], [3, 4]]), [1, 2, 3, 4])

    def test_dissolve_3_list_recursion(self):
        self.assertEqual(anc.dissolve([[[1]]]), [1])

    def test_dissolve_with_tuples(self):
        self.assertEqual(anc.dissolve(((1, 2,), (3, 4))), [1, 2, 3, 4])

    def test_dissolve_tuple_recursion(self):
        self.assertEqual(anc.dissolve(((1, 2), (1, 2))), [1, 2, 1, 2])

    def test_union_(self):
        self.assertEqual(anc.union([1], [1]), [1])


if __name__ == "__main__":
    unittest.main()