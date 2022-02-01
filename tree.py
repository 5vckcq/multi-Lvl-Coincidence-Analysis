#!/usr/bin/env python3

# teilweise uebernommen von https://developpaper.com/question/how-does-python-generate-trees-recursively/

class Tree:

    def __init__(self, name):
        self.name = name
        self.children = {}

    def __iter__(self):
        return iter(self.children)

    def __str__(self):
        return self.name

    def __repr__(self):
        return 'Tree("{}")'.format(self.name)

    def add_child(self, child):
        self.children[child] = child
        return child
    
    
    
    def dfs(self, include_self=True):
        if include_self:
            yield self
        for child in self.children:
            yield child
            yield from child.dfs(False)

    def bfs(self, include_self=True):
        if include_self:
            yield self
        trees = list(self.children.keys())
        while True:
            if not trees:
                break
            tree = trees.pop(0)
            yield tree
            trees += list(tree.children.keys())
