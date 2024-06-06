import numpy as np

class Tree_node:
  def __init__(self, pos_x, pos_y, pos_z):
    self.x = pos_x
    self.y = pos_y
    self.z = pos_z
    self.pos = np.array([self.x, self.y, self.z])

  def __repr__(self):
    return self.pos.__repr__()


class Attraction_point:
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z
    self.pos = np.array([self.x, self.y, self.z])

  def __repr__(self):
    return self.pos.__repr__()
   

class Tree:
  def __init__(self, root):
    self.root = root
    self.nodes = [self.root]
    self.transition_map = {}

  def add_child(self, parent, child):

    if child in self.nodes:
      raise ValueError

    if parent in self.nodes:
      self.transition_map[child] = parent
      self.nodes.append(child)
    else:
      raise ValueError

  def is_leaf(self, node):
    if node not in self.nodes:
      raise ValueError
    if node not in self.transition_map.values():
      return True
    else:
      return False

  def num_children(self, node):
    if node not in self.nodes:
      raise ValueError

    return list(self.transition_map.values()).count(node)

  def get_children(self, parent):
    if parent not in self.nodes:
      raise ValueError

    return [child for child in self.transition_map if self.transition_map[child] == parent]

  def get_level(self, node):
    if node not in self.nodes:
      raise ValueError

    if node == self.root:
      return 0

    x = self.transition_map[node]
    level = 1
    while x != self.root:
      x = self.transition_map[x]

      if self.num_children(x) > 1:
        level += 1

    return level
