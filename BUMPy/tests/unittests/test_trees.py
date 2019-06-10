import unittest
from bumpy import Tree_node, Shape_tree


class test_shape_node(unittest.TestCase):

    def test_constructor_works(self):
        node = Tree_node("name")
        self.assertEqual(node.name, "name")
        self.assertIsNone(node.parent)
        self.assertEqual(node.children, [])

        node2 = Tree_node("name2", node)
        self.assertEqual(node2.parent, node)

    def test_add_children(self):
        node = Tree_node("name")
        node2 = Tree_node("name2")
        node3 = Tree_node("name3")
        node.add_child(node2)
        node.add_child(node3)
        self.assertTrue(node2 in node.children)
        self.assertTrue(node3 in node.children)
        self.assertEqual(node, node2.parent)

    def test_has_child(self):
        node = Tree_node("name")
        node2 = Tree_node("name2")
        node3 = Tree_node("name3")
        node.add_child(node2)
        self.assertTrue(node.has_child(node2))
        self.assertFalse(node.has_child(node3))


class test_shape_tree(unittest.TestCase):

    def test_constructor_works(self):
        tree = Shape_tree("base_shape")
        self.assertEqual(tree.root.name, "base_shape")
        self.assertEqual(len(tree.nodes), 1)

    def test_add_node(self):
        tree = Shape_tree("base_shape")
        tree.add_node(tree.root, "child_node")
        self.assertEqual(len(tree.root.children), 1)
        # make sure base is unchanged
        self.assertEqual(tree.root.name, "base_shape")
        self.assertEqual(tree.root.children[0].name, "child_node")


class test_shape_info(unittest.TestCase):

    def test_register_subshape(self):
        pass

    def test_get_shape_string(self):
        pass
