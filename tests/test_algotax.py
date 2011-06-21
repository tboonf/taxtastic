from Bio import Phylo
from StringIO import StringIO
import unittest
from taxtastic import algotax

def treestring(s):
    return Phylo.read(StringIO(s), 'newick')

class ColoredTreeTestMixin(object):
    @classmethod
    def setup_class(cls):
        cls.parsed_tree = treestring(cls.tree)
        cls.colors = {n: n.name for n in cls.parsed_tree.get_terminals()}
        cls.metadata = algotax.color_clades(cls.parsed_tree, cls.colors)

    @classmethod
    def setUpClass(cls):
        super(ColoredTreeTestMixin, cls).setUpClass()
        cls.setup_class()

class CladeColorTestMixin(ColoredTreeTestMixin):
    def test_parents(self):
        self.assertEqual(
            len(self.metadata.parents),
            sum((1 for _ in self.parsed_tree.find_clades())))
        for child, parent in self.metadata.parents.iteritems():
            if parent is None:
                self.assertIs(child, self.parsed_tree.root)
            else:
                self.assertIn(child, parent.clades)

    def test_cut_colors(self):
        numbering = dict(enumerate(
            self.parsed_tree.find_clades(order='postorder')))
        rev_numbering = {v: k for k, v in numbering.iteritems()}
        for color, nodes in self.cut_colors.iteritems():
            if color is None:
                for node in nodes:
                    cut_colors = self.metadata.cut_colors[numbering[node]]
                    self.assert_(not cut_colors)
                continue

            for node in numbering.itervalues():
                if node is self.parsed_tree.root:
                    continue
                self.assertEqual(
                    rev_numbering[node] in nodes,
                    color in self.metadata.cut_colors[node])

class CladeColorTest1(CladeColorTestMixin, unittest.TestCase):
    tree = '((A,A),(B,B))'
    cut_colors = {
        'A': {0, 1},
        'B': {3, 4},
        None: {2, 5},
    }

class CladeColorTest2(CladeColorTestMixin, unittest.TestCase):
    tree = '((A,B),((A,B),A))'
    cut_colors = {
        'A': {0, 2, 7, 5, 3, 6},
        'B': {1, 2, 7, 5, 4},
    }

class CladeColorTest3(CladeColorTestMixin, unittest.TestCase):
    tree = '(((A,B),(C,D)),(A,B),(C,A))'
    cut_colors = {
        'A': {0, 2, 6, 9, 7, 12, 11},
        'B': {1, 2, 6, 9, 8},
        'C': {3, 5, 6, 12, 10},
        'D': set(),
        None: {4},
    }

class AlgotaxWalkTestMixin(ColoredTreeTestMixin):
    @classmethod
    def setup_class(cls):
        super(AlgotaxWalkTestMixin, cls).setup_class()
        cls.nodeset = algotax.walk(cls.parsed_tree.root, cls.metadata)

    setUpClass = setup_class

    def test_walk(self):
        self.assertEqual(len(self.nodeset), self.convex_tree_size)

class AlgotaxWalkTest1(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '((A,A),(B,B))'
    convex_tree_size = 4

class AlgotaxWalkTest2(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '((A,B),(A,B))'
    convex_tree_size = 3

class AlgotaxWalkTest3(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '(((A,B),(A,B)),(C,C))'
    convex_tree_size = 5

class AlgotaxWalkTest4(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '(((A,B),B),(A,A))'
    convex_tree_size = 4

class AlgotaxWalkTest5(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '(A,(A,(B,C)))'
    convex_tree_size = 4

def subrk_min_of_map(subrk_map):
    def subrk_min(terminals):
        names = {n.name for n in terminals}
        return subrk_map.get(''.join(sorted(names)), 0)
    return subrk_min

def node_number(tree, node):
    return next(e
        for e, n in enumerate(tree.find_clades(order='postorder'))
        if n is node)

class RerootingTests(unittest.TestCase):
    def test_single_root(self):
        tree = treestring('(A,(B,C))')
        subrk_map = {
            'A': 1,
            'B': 1,
            'C': 2,
        }
        first_root, other_roots = algotax.reroot(
            tree.root, subrk_min_of_map(subrk_map))
        self.assertFalse(other_roots)
        self.assertEqual(node_number(tree, first_root), 3)

    def test_two_roots(self):
        tree = treestring('(A,B,(C,D))')
        subrk_map = dict.fromkeys('ABCD', 1)
        first_root, other_roots = algotax.reroot(
            tree.root, subrk_min_of_map(subrk_map))
        self.assertTrue(other_roots)
        self.assertIn(node_number(tree, first_root), {4, 5})
        self.assertEqual(len(other_roots), 1)
        other_root, = other_roots
        self.assertIn(node_number(tree, other_root), {4, 5})
        self.assertNotEqual(first_root, other_root)

    def test_three_roots(self):
        tree = treestring('(A,B,(C,(D,(E,(F,G)))))')
        subrk_map = dict.fromkeys('ABCDEFG', 1)
        first_root, other_roots = algotax.reroot(
            tree.root, subrk_min_of_map(subrk_map), False)
        self.assertTrue(other_roots)
        self.assertIn(node_number(tree, first_root), {8, 9, 10})
        self.assertEqual(len(other_roots), 2)
        first_other_root, second_other_root = other_roots
        self.assertIn(node_number(tree, first_other_root), {8, 9, 10})
        self.assertIn(node_number(tree, second_other_root), {8, 9, 10})
        self.assertNotEqual(first_root, first_other_root)
        self.assertNotEqual(first_root, second_other_root)
