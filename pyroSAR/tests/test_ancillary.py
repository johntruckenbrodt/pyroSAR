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

    def test_dictmerge(self):
        self.assertEqual(anc.dictmerge({'a': 1, 'b': 2}, {'c': 3, 'd': 4}), {'a': 1, 'b': 2, 'c': 3, 'd': 4})

    def test_parse_literal(self):
        self.assertEqual(anc.parse_literal(['1', '2.2', 'a']), [1, 2.2, 'a'])

    def test_parse_literal_error(self):
        with self.assertRaises(IOError):
            anc.parse_literal(1)

    def test_seconds(self):
        self.assertEqual(anc.seconds('test_20151212T234411'), 3658952651.0)


if __name__ == "__main__":
    unittest.main()
