## Software Pre Requisites ##
Requires: [Shapely](https://pypi.python.org/pypi/Shapely)

On Ubuntu:
```
#!shell
apt-get install libgeosX.X
```
```
#!shell
apt-get install python-shapley
```

On Mac OS:
```
#!shell
sudo port install geos
```
```
#!shell
sudo port install py-shapely
```



# Quad Tree Exercise #

The quad tree object oriented implementation we propose, relies on two objects, ``Node`` and ``QuadTree``.
In the following we assume ``python`` as programming language. Of course the same ideas apply independently of the language.

##Node##

Node represents a single node on a quad tree. It is characterized by the knowledge of several informations:

```
#!python
class Node:
	def __init__(self,p_in):
		self.polygon = p_in
		self.population = []
		self.kind = ''
```
Let's first focus on the type of information and then we discuss the meaning of a function ``__init__()``.

* ``polygon`` the spatial representation of the tree node.
* ``population`` is the list of particles (or targets) contained in the ``Node``.
* ``kind`` it is a string that defines if the node is ``"root"``, ``"node"``, or ``"leaf"``.

The strange function ``__init__()``  is the definition for the **Constructor**. 
At this level we introduce two concepts: the *Constructor* itself and its ``python``
syntax. Focus first on concepts. When defining (or *declaring*) an object inside a main file
we need to *construct* it. Then **Constructors** are a particular kind of functions that 
*define*, *declare*, *construct* (call this action the way you like, the technical expression 
is **declare**) an object. The result in the main application will be:

```
#!python
node = Node(polygon)
```
Now that we know what is a constructor, it is straight to understand the syntax. 
``__init__()`` is simply the function that initializes the object. The other strange word is 
``self``. You can read it as "the object itself". If you use it as a prefix before a variable, 
for example ``self.population``, it makes the variable **public**, meaning that now that variable can 
be accessed bot from inside the object, and from outside. In a Fortran style we are declaring module 
variables, being the module ``Node``. When ``self`` is used as first argument of a function, defines
the function as **member function**, meaning a function that is member of the object.

This a good point to introduce *which* member functions are useful to a node. In particular we focused on 
two member functions. The first one is devoted to fill the list of particles contained in the 
``population`` public variable. It is designed to get a list of particles, and define the subset of these 
particles that populate the node. The *interface* is:
```
#!python
	def evaluate_population(self,target_in):
            ...
            return
```
The second member function subdivides the node into for equal children.  No surprise it is called ``subdivide``.  In this case
as well it is pretty easy to understand that the function simply returns the list of the four children:
```
#!python
      def subdivide_into_poligons(self):
             ...
	     return children 
```
In this function name we specify ``subdivide _into_polygons`` because we designed this function to produce nodes where only the 
``polygon`` member variable is specified. In a sense, the ``Node`` only handles local information. It is left to the ``QuadTree`` to 
correctly locate the node along the tree and complete with the necessary information.

## QuadTree ##

The second object is ``QuadTree``. In the suggested design we *encapsulate* the local information and functions in the object 
``Node``, while  ``QuadTree`` is in charge of information and methods that involve the whole tree. When constructed, the ``QuadTree`` object, 
builds the node corresponding to the whole domain, and gets the list of targets (particles to be located). The constructed ``node`` spans over the whole unit box, and its ``population`` is the whole set of particles. By definition this is a ``"root"`` node, and we append it to the public variable ``node_list`` which is the *goal* of our algorithm.
```
#!python
class QuadTree:
	node_tree = []
	def __init__(self,target_in):
		node = Node(Polygon([(0,0),(1,0),(1,1),(0,1)]))
		node.population = target_in
		node.kind = 'root'
		self.node_tree.append(node)
		return
```
The first thing that  *tree* naturally does, it to *bloom*, we obtain a list of four children. As we plan to do some more manipulations on the ``children`` list, we declare it to be public with the ``self.`` prefix.
```
#!python
	def bloom_children(self,parent):

		self.children = parent.subdivide_into_poligons()
```
Now for every ``kid`` in this list ``QuadTree`` completes the set of information (the missing information are: ``population`` and ``kind`` check the node definition). Using the member function  ``evaluate_population``
we define the ``kid.population`` being a subset of ``parent.population``. Now if the number of particles inside the ``kid`` is larger than one, the ``kid`` is defined to be a ``"node"``. If the number of particles inside the ``kid`` is equal to one, the ``kid`` is a ``"leaf"``. 
```
#!python
		for kid in self.children:

			kid.evaluate_population(parent.population)
			n_targets =  len(kid.population)

			if n_targets > 1:
				kid.kind = 'node'
			if n_targets == 1:
				kid.kind = 'leaf'
```


It is a natural selection matter that not all of the children shall survive. Then we need a nasty member function, called ``prune_children``. 
In this function only ``"leaf"`` and ``"node"`` nodes will  be appended to the tree. When the node is a ``"node"`` we ask the tree itself to 
``refine``.
```
#!python
	def prune_children(self):

		for node in self.children:
			if node.kind== 'node':
			    self.node_tree.append(node)
				self.refine(node)
			if node.kind== 'leef':
				self.node_tree.append(node)

		return
```
If the node is a ``node`` we will ``refine`` it. At this point the reader should have realized that this function
refine is recursive, and repeats the processes of *blooming* and *pruning*:
```
#!python
	def refine(self,node):

		self.bloom_children(node)
		self.prune_children()

		return
```

##The Main File##
As usual in object oriented practice, the main file becomes no more than a few lines:
```
#!python
def main():

	x = np.array([.1,.101,.11,.9,.91])
	y = np.array([.9,.9,.7,.1,.09])

	coords_list = []
	for ix, iy in zip(x,y):
		coords = (ix,iy)
		coords_list.append(coords)

	quad_tree = QuadTree(coords_list)
	quad_tree.refine(quad_tree.get_root_node())

	quad_tree.plot()
	
	return

if __name__ == '__main__':
	main()
```
The reader may notice, that two elementary methods are not discussed here, ``get_root_node`` and ``plot``. It is opinion of the author that no detailed discussion is needed about these two methods.