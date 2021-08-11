'''
This is a solved problem from Rosalind website on bioinformatics problems.
Title: Introduction to Pattern Matching
ID: TRIE

Given: A list of at most 100 DNA strings of length at most 100 bp, none of which is a prefix of another.
Return: The adjacency list corresponding to the trie T for these patterns, in the following format. If T has n nodes,
first label the root with 1 and then label the remaining nodes with the integers 2 through n in any order you like. Each
edge of the adjacency list of T will be encoded by a triple containing the integer representing the edge's parent node,
followed by the integer representing the edge's child node, and finally the symbol labeling the edge.
'''


from pptree import *


class Tree:

    def __init__(self, data, edge='', name='root'):
        self.data = data
        self.name = f'-{edge}- {name} ({data})'
        self.edge = edge
        self.children = []

    def set_traits(self, data = None, edge = None, name = None):
        if data is not None:
            self.data = data
        if edge is not None:
            self.edge = edge
        if name is not None:
            self.name = name

    def set_children(self, *children):
        for child in children:
            if type(child) != list:
                self.children.append(Tree(data = child, name='child'))
                ch = self.get_children(-1)
                ch.set_traits(name = f'-- child[{self.children.index(ch)}] ({child})')
            else:
                self.children.append(Tree(data = child[0], edge = child[1], name='child'))
                ch = self.get_children(-1)
                ch.set_traits(name = f'-{child[1]}- child[{self.children.index(ch)}] ({child[0]})')

    def get_traits(self, trait=None):
        if trait is not None:
            if trait.lower() == 'data':
                return self.data
            if trait.lower() == 'edge':
                return self.edge
            if trait.lower() == 'name':
                return self.name
        return None

    def get_children(self, ind=None, data=None, edge=None):
        if ind is not None:
            return self.children[ind]

        if edge is not None and data is not None:
            children_list = []
            for child in self.children:
                if edge == child.get_traits('edge') and data == child.get_traits('data'):
                    children_list.append(child)
            if len(children_list) == 1:
                return children_list[0]
            return children_list

        if edge is not None:
            children_list = []
            for child in self.children:
                if edge == child.get_traits('edge'):
                    children_list.append(child)
            if len(children_list) == 1:
                return children_list[0]
            return children_list

        if data is not None:
            children_list = []
            for child in self.children:
                if data == child.get_traits('data'):
                    children_list.append(child)
            if len(children_list) == 1:
                return children_list[0]
            return children_list
        return self.children

    def connect(self, tree):
        self.children.append(tree)

    def print_tree(self):
        return print_tree(self)


def build_trie(text, tree=Tree('')):

    """    Input: a text of DNA sequences lines (text)
    Output: a trie containing all of the DNA bases of the given DNA sequences lines with their edges labeled in
    ascending order of 1 starting from the number 2"""

    trie = tree
    pointer = trie
    overlap = True
    edge = 2
    for i in text.replace('\n', '$'):
        if i != '$':
            overlap = False  # there is no overlap until found otherwise
            for child in pointer.get_children():
                # if the character is equal to the data of the child, go to that child to check for further overlap
                if i == child.get_traits('data'):
                    pointer = child
                    overlap = True
                    break
            # if overlap was not found in the children, add the rest of the characters in the line to the tree
            if not overlap:
                pointer.set_children(i)
                pointer = pointer.get_children(ind=-1)
                pointer.set_traits(edge=edge)
                # changing the name for visualisation reasons when printing the tree
                pointer.set_traits(name='-' + str(edge) + pointer.get_traits('name')[1:])
                edge += 1
        else:
            # stat from the beginning for the next line
            pointer = trie
            overlap = True
    return trie


def pattern_matching(trie, count=0, edge=0):

    """    Input: a trie (trie)
    Output: a string (pattern) with n lines for each of the n nodes in the trie (except the root) where each line have
    the edge of the parent node, followed by the edge of the node itself and the node data (the DNA base). If the parent
    is the root of the trie, the parent edge would be represented as 1"""

    pattern = ''
    pointer = trie
    edge = pointer.get_traits('edge')
    for ind in range(len(pointer.get_children())):
        if edge == '':
            pattern += '1 '
        else:
            pattern += str(edge) + ' '
        pattern += str(pointer.get_children(ind=ind).get_traits('edge')) + ' '
        pattern += pointer.get_children(ind=ind).get_traits('data') + '\n'
        pattern += pattern_matching(pointer.get_children(ind=ind))

    return pattern


with open('...//rosalind_trie.txt', 'r') as file:
    f = file.read()

trie = build_trie(f)
print(pattern_matching(trie))
#t.print_tree()
