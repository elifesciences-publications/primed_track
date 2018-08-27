# -*- coding: utf-8 -*-

# Calculate total distance between lineage trees.
#
# Code for the paper:
#
# Welling et al. "High fidelity lineage tracing in mouse pre-implantation 
# embryos using primed conversion of photoconvertible proteins".
#
# This Imaris XTension required pIceImarisConnector to run.
# See: https://github.com/aarpon/pIceImarisConnector
#
# <CustomTools>
#  <Menu>
#   <Submenu name="BSSE">
#    <Submenu name="Lineage Tracer">
#     <Item name="Calculate total distance between lineage trees" icon="Python">
#      <Command>PythonXT::bsseLineageTracerTotalDistanceBetweenLineageTrees(%i)</Command>
#     </Item>
#    </Submenu>
#   </Submenu>
#  </Menu>
# </CustomTools>
#
# Aaron Ponti (BSSE) 2017, 2018

from pIceImarisConnector import pIceImarisConnector as pIce
import numpy as np
import treelib
import zss
import scipy.spatial
from Tkinter import *
from tkSimpleDialog import Dialog

# SET THE PARAMETERS FOR THE TREE DISTANCE CALCULATION. 
PARAMETERS = dict()

# Condense trees for distance calculation? A condensed tree is a tree that
# preserves only the branching events. All "linear stretches" are condensed;
# however, the information on how many nodes were present between two branch
# points is not discared. It is stored in the node .data property and is used
# to calculate the distance between nodes using the 'use_diff_num_children'
# strategy below.
PARAMETERS['condense'] = True

# Distance between two nodes: if True, return distance between two nodes as
# the difference the original number of children before tree condensation;
# if False, return 1 if the original numbers of children before tree
# condensation differ, 0 otherwise.
PARAMETERS['use_diff_num_children'] = False


class ListSingleSelectDialog(Dialog):
    """
    List dialog (single selection).
    """

    def __init__(self, title, prompt, choice_list, parent=None):

        if not parent:
            import Tkinter
            parent = Tkinter._default_root
            parent.withdraw()

        self.__prompt = prompt
        self.__choices = choice_list
        self.__retval = None
        self.__clickedOK = False
        parent.update()
        Dialog.__init__(self, parent, title) # enters main loop here

    def get_current_index(self):
        """ Return currently selected index (or -1) """

        # Need to convert to int; currently API returns a tuple of string
        cur_sel = self.__lb.curselection()
        if cur_sel and len(cur_sel) > 0:
            return int(cur_sel[0])
        else:
            return -1

    def getresult(self): return self.__retval

    def destroy(self):
        # first save the selected index before it is destroyed
        idx = self.get_current_index()
        # in PyRAF, assume they meant the first one if they clicked nothing,
        # since it is already active (underlined)
        if idx < 0: idx = 0
        # get the object at that index
        if self.__clickedOK and idx >= 0: # otherwise is None
            self.__retval = self.__choices[idx]
        if self.__retval and type(self.__retval) == str:
            self.__retval = self.__retval.strip()

        # now destroy
        self.__lb = None
        Dialog.destroy(self)

    def body(self, master):

        self.geometry("600x480")
        label = Label(master, text=self.__prompt, justify=LEFT)
#       label.grid(row=0, padx=8, sticky=W)
        label.pack(side=TOP, fill=X, padx=10, pady=8)

        frame = Frame(master)
#       frame.grid(row=1, padx=8, sticky=W+E)
        frame.pack(side=TOP, fill=X, padx=10, pady=8)

        vscrollbar = Scrollbar(frame, orient=VERTICAL)
        hscrollbar = Scrollbar(frame, orient=HORIZONTAL)
        self.__lb = Listbox(frame,
                            width=80, height=20,
                            selectmode=BROWSE,
                            xscrollcommand=hscrollbar.set,
                            yscrollcommand=vscrollbar.set)
#                           activestyle='none', # none = dont underline items
        hscrollbar.config(command=self.__lb.xview)
        hscrollbar.pack(side=BOTTOM, fill=X)
        vscrollbar.config(command=self.__lb.yview)
        vscrollbar.pack(side=RIGHT, fill=Y)
        self.__lb.pack(side=LEFT, fill=BOTH, expand=1)

        for itm in self.__choices:
            self.__lb.insert(END, str(itm))

        self.__lb.bind("<Double-Button-1>", self.ok) # dbl clk
#       self.__lb.selection_set(0,0)
        self.__lb.focus_set()

        return self.__lb

    def ok(self, val=None):
        self.__clickedOK = True # save that this wasn't a cancel
        Dialog.ok(self, val)

    def validate(self): return 1


class CustomNode(treelib.Node):
    """
    A node of the lineage tree.
    """
    def __init__(self, tag=None, identifier=None, expanded=True, data=None, treeref=None):
        super(CustomNode, self).__init__(tag, identifier, expanded, data)
        self.treeref = treeref

    def __len__(self):
        return len(self.treeref.nodes)

    def get_children(self):
        return self.treeref.children(self.identifier)


def appendChildren(treelib_node, zss_node):
    # Get the children of the node
    for child in treelib_node.get_children():
        zss_child_node = zss.Node(child.tag)
        appendChildren(child, zss_child_node)
        zss_node.addkid(zss_child_node)


def zssTreeFromTreeLibTree(tree):
    zss_root = zss.Node(tree.get_node(tree.root))
    appendChildren(tree.get_node(tree.root), zss_root)
    return zss_root    
    

def bsseLineageTracerTotalDistanceBetweenLineageTrees(aImarisId):
    """
    ImarisXT to calculate the total distance between lineage trees constructed
    from tracked Imaris ISpots objects.
    :param aImarisId: Imaris ID of pIceImarisConnector object
    :type aImarisId: {int|pIceImarisConnector}
    """
    tk_root = Tk()
    tk_root.withdraw()
    tk_root.update()
    
    # Instantiate IceImarisConnector object
    conn = pIce(aImarisId)

    # Is there something loaded?
    mDataSet = conn.mImarisApplication.GetDataSet()
    if mDataSet is None:
        return

    # Get the spots
    spots = conn.getAllSurpassChildren(True, 'Spots')
    if spots is None or len(spots) == 0:
        print('No Spots objects found in the scene!')
        return

    if len(spots) < 2:
        print('You need two Spots objects!')
        return

    # Get spot names
    spotNames = []
    for i in range(len(spots)):
        spotNames.append(spots[i].GetName())

    # Ask the user to specify the first Spots object
    s = ListSingleSelectDialog("Question",
                               "Please pick the first Spots object",
                               spotNames, None)
    
    if s.getresult() is None:
        return
    
    firstSpotName = s.getresult()
    spot1 = spots[spotNames.index(firstSpotName)]
    print("First spot object: '%s'." % firstSpotName)
    
    # Ask the user to specify the first Spots object
    s = ListSingleSelectDialog("Question",
                               "Please pick the second Spots object",
                               spotNames, None)
   
    if s.getresult() is None:
        return

    secondSpotName = s.getresult()
    spot2 = spots[spotNames.index(secondSpotName)]
    print("Second spot object: '%s'." % secondSpotName)

    # Build the lineage trees
    print("Building the lineage tree for '%s'..." % firstSpotName)
    lineage1 = buildLineage(spot1)
    print("Lineage for '%s' has %d independent trees." % (firstSpotName, len(lineage1)))
    
    print("Building the lineage tree for '%s'..." % secondSpotName)
    lineage2 = buildLineage(spot2)
    print("Lineage for '%s' has %d independent trees." % (secondSpotName, len(lineage2)))

    # Compare the trees
    print("Calculating total distance...")
    distance = calculateLineageDistance(lineage1, lineage2, spot1, spot2)
    
    if distance == 0:
        print("The lineage trees are identical (distance is %.0f)." % distance)
    else:
        print("Total distance is %.0f." % distance)

    # Ask the user to press a key to close the window (and return to Imaris)
    raw_input("Press any key to continue.")


def buildLineage(spot):
    """
    Build lineage tree for spot.
    :param spot: an Imaris ISpots object.
    :type spot: ISpots
    :return: lineage tree
    :rtype: list
    """
    
    # Get spot Ids
    spotIds = np.array(spot.GetIds())

    # Get spot information
    spotTimes = np.array(spot.GetIndicesT())

    # Get track ids and egdes
    trackIds = np.array(spot.GetTrackIds())
    uTrackIds = np.unique(trackIds)
    trackEdges = np.array(spot.GetTrackEdges())

    # We create a tree per track id and store the timepoint of origin in
    # the root of the tree
    trees = []

    for i in range(len(uTrackIds)):

        currTrackId = uTrackIds[i]
        currTrack,  = np.where(trackIds == currTrackId)
        currTrackEdges = trackEdges[currTrack, :]

        # Get the timepoint of origin
        spotIndx, = np.where(spotIds == currTrackEdges[0, 0])
        if spotIndx.shape[0] == 0:
            print("There is a mismatch between tracks and spot IDs: re-run the tracker!")
            continue
        originTimepoint = spotTimes[spotIndx][0]
    
        # We create a tree
        lineage = treelib.Tree()
        
        # As root node we set the timepoint
        root = CustomNode(tag=-1,
                         identifier=-1,
                         data=originTimepoint,
                         treeref=lineage)
        lineage.add_node(root)

        for j in range(currTrackEdges.shape[0]):
            
            currentEdge = currTrackEdges[j, :]
        
            if j == 0:
                parent = CustomNode(tag=currentEdge[0],
                         identifier=currentEdge[0],
                         treeref=lineage)
                lineage.add_node(parent, parent=root)
            else:
                parent = lineage.get_node(currentEdge[0])
                if parent is None:
                    parent = root

            # Add a node and update lastnode
            node = CustomNode(tag=currentEdge[1],
                              identifier=currentEdge[1],
                              treeref=lineage)
            lineage.add_node(node, parent=parent)

        trees.append(lineage)

    return trees


def calculateLineageDistance(lineage1, lineage2, spot1, spot2):
    """
    Calculate distance between lineage trees.
    :param lineage1: first lineage tree
    :type lineage1: list
    :param lineage2: second lineage tree
    :type lineage2: list
    :param spot1: first Imaris ISpots object
    :type spot1: ISpots
    :param spot2: second Imaris ISpots object
    :type spot2: ISpots
    :return: distance between lineage trees
    :rtype: float64
    """

    # Some constants
    numLineage1 = len(lineage1)
    numLineage2 = len(lineage2)
    
    if len(lineage1) == 0 or len(lineage2) == 0:
        return 0.0

    # Get spot ids
    spotId1 = np.array(spot1.GetIds())
    spotId2 = np.array(spot2.GetIds())

    # Get spot coordinates
    spotPosXYZ1 = np.array(spot1.GetPositionsXYZ())
    spotPosXYZ2 = np.array(spot2.GetPositionsXYZ())

    # Get the time points of origin and the tree length
    matrixLineage1 = np.zeros((numLineage1, 4))
    matrixLineage2 = np.zeros((numLineage2, 4))

    for i in range(numLineage1):
        tree = lineage1[i]
        node = tree.get_node(tree.root)
        timepoint = node.data
        s_id = node.get_children()[0].tag
        coords = spotPosXYZ1[spotId1 == s_id, :][0] 
        matrixLineage1[i, 0] = timepoint
        matrixLineage1[i, 1:] = coords

    for i in range(numLineage2):
        tree = lineage2[i]
        node = tree.get_node(tree.root)
        timepoint = node.data
        s_id = node.get_children()[0].tag
        coords = spotPosXYZ2[spotId2 == s_id, :][0] 
        matrixLineage2[i, 0] = timepoint
        matrixLineage2[i, 1:] = coords

    # Tree pairs to process
    nTreePairsToProcess = np.min([numLineage1, numLineage2])

    # Find the matching
    all_distances = np.zeros((np.max([numLineage1, numLineage2])))
    D = scipy.spatial.distance.cdist(matrixLineage1, matrixLineage2,
                                     'sqeuclidean')
    nAssigned = 0
    while nAssigned < nTreePairsToProcess:
        p = np.where(D == np.min(D))
        y = p[0][0]
        x = p[1][0]
        all_distances[nAssigned] = compareSingleTrees(lineage1[y], lineage2[x])
        D[y, x] = np.Inf
        nAssigned = nAssigned + 1
   
    # Total distance
    total_distance = np.sum(all_distances)
    
    return total_distance


def compareSingleTrees(tree1, tree2):
    """
    Compare two individual trees from a lineage tree (that may have several
    'trunks'.
    :param tree1: first tree.
    :type tree1: list
    :param tree2:  second tree
    :type tree2: list
    :return: ditance between the trees
    :rtype: float64
    """
    # Condense the trees
    if PARAMETERS['condense']:
        c_tree1 = condense_tree(tree1)
        c_tree2 = condense_tree(tree2)
        root1 = c_tree1.get_node(c_tree1.root)
        root2 = c_tree2.get_node(c_tree2.root)
        
        cmpFnc = compare_num_children

    else:
        root1 = tree1.get_node(tree1.root)
        root2 = tree2.get_node(tree2.root)
        
        cmpFnc = lambda a, b: 0
    
    return zss.compare.distance(root1,
                                root2,
                                CustomNode.get_children,
                                insert_cost=lambda node: 1,
                                remove_cost=lambda node: 1,
                                update_cost=cmpFnc)

    #    zss_tree1 = zssTreeFromTreeLibTree(tree1)
    #    zss_tree2 = zssTreeFromTreeLibTree(tree2)
    #    
    #    return zss.compare.distance(zss_tree1, zss_tree2,
    #                                zss.simple_tree.Node.get_children,
    #                                insert_cost=lambda node: 1,
    #                                remove_cost=lambda node: 1,
    #                                update_cost=lambda a, b: 0)


def is_branch(node):
    """
    Test is a node is a branch node.
    :param node: a node in the lineage tree.
    :type node: CustomNode
    :return: True if the node is a branch node.
    :rtype: Boolean
    """
    return len(node.treeref.is_branch(node.identifier)) > 1


def compare_num_children(n1, n2):
    """
    Compare the number of children.

    If PARAMETERS['use_diff_num_children'] is True, the function returns the
    difference in the number of children between the nodes.

    If PARAMETERS['use_diff_num_children'] is False, the function returns 1
    if the nodes have a different number of children and 0 otherwise.

    :param n1: node 1
    :type n1: CustomNode
    :param n2: node 2
    :type n2: CustomNode
    :return: distance between two nodes
    :rtype: int
    """

    if PARAMETERS['use_diff_num_children']:
        return abs(n1.data - n2.data)
    else:
        if n1.data == n2.data:
            return 0
        else:
            return 1


def condense_tree(tree):
    """
    Condense a tree by preserving only the branching events. All "linear
    stretches" are condensed; however, the information on how many nodes
    were present between two branch points is not discared. It is stored
    in the node .data property.
    :param tree: a tree
    :type tree: list
    :return: condensed tree
    :rtype: list
    """
    __DEBUG__ = False
    
    # Root of current tree
    c_root = tree.get_node(tree.root)
    
    # Prepare the new tree to build
    new_tree = treelib.Tree()
    n_root = CustomNode(tag=c_root.tag, 
                        identifier=c_root.identifier, 
                        treeref=new_tree)
    n_root.data = 0
    new_tree.add_node(n_root)  

    # Get the paths to the leaves
    paths = tree.paths_to_leaves()
    
    processed_non_branch_nodes = []
    
    # Process the paths
    for path in paths:
    
        if __DEBUG__:
            print("")
            print("* * * * * * * * * * * * * * * * * * * * * * * * * *")
            print("Processing path: " )
            print(path)

        # Process the nodes and add the branch points to the new
        # tree storing (and updating) the number of children to the next 
        # branch point or the leaf
        c_parent = c_root
        n_parent = n_root
        for i in range(len(path)):

            # Get the node id
            c_node_id = path[i]

            # Get the actual node
            c_node = tree.get_node(c_node_id)

            if __DEBUG__:
                print("")
                print("Processing node %s" % c_node_id)
            
            # Test if it is a branch
            if is_branch(c_node):

                if __DEBUG__:
                    print("Node %s is a branch point" % c_node_id)

                # Since it is a branch, it may have already been added
                if new_tree.contains(c_node.identifier):

                    # Get the node
                    n_node = new_tree.get_node(c_node.identifier)

                    if __DEBUG__:
                        print("Node %s is already in the tree (currently with %d children)" % (n_node.identifier, n_node.data))

                else:
                
                    # Create the node
                    n_node = CustomNode(tag=c_node.tag, 
                                       identifier=c_node.identifier, 
                                       treeref=new_tree)
                    n_node.data = 0
                    
                    # Add to the new tree
                    new_tree.add_node(n_node, parent=n_parent)    
    
                    if __DEBUG__:
                        print("Node %s is created and added to the tree" % n_node.identifier)

                    # Add a child to the parent count
                    if n_parent.data is None:
                        n_parent.data = 0
                    n_parent.data = n_parent.data + 1

                    if __DEBUG__:
                        print("One child is added to the parent node %s (that has now %d children)" % (n_parent.identifier, n_parent.data))
                    
                # Update the parents
                n_parent = n_node
                c_parent = c_node
            
            else:

                # If we haven't visited this node yet, we update the 
                # child count of its parent
                if c_node_id not in processed_non_branch_nodes:
                    
                    if n_parent.data is None:
                        n_parent.data = 0
                    n_parent.data = n_parent.data + 1
                    
                    if __DEBUG__:
                        print("Node %s is not a branch point. Added 1 child to its parent %s (that has now %d children)" % (c_node_id, n_parent.identifier, n_parent.data))
                    
                    # Mark this node as already processed
                    processed_non_branch_nodes.append(c_node_id)
                
                else:
                    
                    if __DEBUG__:
                        print("Node %s is not a branch point but was already processed" % (c_node_id))
                
    return new_tree


def test():
    
    from zss import Node
    import zss.compare

    #
    # Test
    #
    # The trees have the same structure. Only one label is different.
    A = (
            Node("f")
            .addkid(Node("a")
                .addkid(Node("d"))
                .addkid(Node("c")
                    .addkid(Node("b"))))
            .addkid(Node("e"))
        )
        
    B = (
            Node("f")
            .addkid(Node("m")            # <- Label change!
                .addkid(Node("d"))
                .addkid(Node("c")
                    .addkid(Node("b"))))
            .addkid(Node("e"))
        )
    
    # The simple_distance() function should report one update (distance = 1).
    zss_distance = zss.compare.simple_distance(A, B)
    print("Test: %d == 1" % (zss_distance))
    assert(zss_distance == 1)
    
    zss_distance = zss.compare.distance(A, B, zss.simple_tree.Node.get_children,
                                        insert_cost=lambda node: 1,
                                        remove_cost=lambda node: 1,
                                        update_cost=lambda a, b: 0)

    # The distance() function with no cost for one update (distance = 0).
    print("Test: %d == 0" % (zss_distance))
    assert(zss_distance == 0)

    #
    # Test
    #
    # The second tree is missing one node.
    A = (
            Node("f")
            .addkid(Node("a")
                .addkid(Node("d"))
                .addkid(Node("c")
                    .addkid(Node("b"))))
            .addkid(Node("e"))
        )

    B = (
            Node("f")
            .addkid(Node("a")
                .addkid(Node("c")
                    .addkid(Node("b"))))
            .addkid(Node("e"))
        )

    # The simple_distance() function should report one insert (distance = 1).
    zss_distance = zss.compare.simple_distance(A, B)
    print("Test: %d == 1" %(zss_distance))
    assert(zss_distance == 1)
    
    # The distance() function should report one insert (distance = 1).
    zss_distance = zss.compare.distance(A, B, zss.simple_tree.Node.get_children,
                                        insert_cost=lambda node: 1,
                                        remove_cost=lambda node: 1,
                                        update_cost=lambda a, b: 0)

    print("Test: %d == 1" % (zss_distance))
    assert(zss_distance == 1)

    #
    # Test
    #
    # The trees are identical as the ones in the previous test, but the 
    # distance() works on the libtree objects directly (NO conversion).
    A2 = treelib.Tree()
    a_root = CustomNode(tag="f", identifier="f", treeref=A2)
    A2.add_node(a_root)  
    a_an = CustomNode(tag="a", identifier="a", treeref=A2)
    A2.add_node(a_an, parent=a_root)    
    a_dn = CustomNode(tag="d", identifier="d", treeref=A2)
    A2.add_node(a_dn, parent=a_an)    
    a_cn = CustomNode(tag="c", identifier="c", treeref=A2)
    A2.add_node(a_cn, parent=a_an)    
    a_bn = CustomNode(tag="b", identifier="b", treeref=A2)
    A2.add_node(a_bn, parent=a_cn)    
    a_en = CustomNode(tag="e", identifier="e", treeref=A2)
    A2.add_node(a_en, parent=a_root)    

    B2     = treelib.Tree()
    b_root = CustomNode(tag="f", identifier="f", treeref=B2)
    B2.add_node(b_root)  
    b_an = CustomNode(tag="a", identifier="a", treeref=B2)
    B2.add_node(b_an, parent=b_root)    
    b_cn = CustomNode(tag="c", identifier="c", treeref=B2)
    B2.add_node(b_cn, parent=b_an)    
    b_bn = CustomNode(tag="b", identifier="b", treeref=B2)
    B2.add_node(b_bn, parent=b_cn)    
    b_en = CustomNode(tag="e", identifier="e", treeref=B2)
    B2.add_node(b_en, parent=b_root)    
    
    libtree_distance = zss.compare.distance(a_root, 
                                            b_root,
                                            CustomNode.get_children,
                                            insert_cost=lambda node: 1,
                                            remove_cost=lambda node: 1,
                                            update_cost=lambda a, b: 0)

    # The distance must the the same as in the previous test (distance = 1)    
    print("Test: %d == %d" %(libtree_distance, zss_distance))
    assert(libtree_distance == zss_distance)
    
    del A
    del B
    del a_root, a_an, a_dn, a_cn, a_bn, a_en
    del b_root, b_an, b_cn, b_bn, b_en
    
    
    #
    # Test
    #
    # Compare two identical trees, but one a zss tree and the other as 
    # converted treelib tree. The distance must be zero.
    A = (
            Node("f")
            .addkid(Node("a")
                .addkid(Node("d"))
                .addkid(Node("c")
                    .addkid(Node("b"))))
            .addkid(Node("e"))
        )

    A2 = treelib.Tree()
    a_root = CustomNode(tag="f", identifier="f", treeref=A2)
    A2.add_node(a_root)  
    a_an = CustomNode(tag="a", identifier="a", treeref=A2)
    A2.add_node(a_an, parent=a_root)    
    a_dn = CustomNode(tag="d", identifier="d", treeref=A2)
    A2.add_node(a_dn, parent=a_an)    
    a_cn = CustomNode(tag="c", identifier="c", treeref=A2)
    A2.add_node(a_cn, parent=a_an)    
    a_bn = CustomNode(tag="b", identifier="b", treeref=A2)
    A2.add_node(a_bn, parent=a_cn)    
    a_en = CustomNode(tag="e", identifier="e", treeref=A2)
    A2.add_node(a_en, parent=a_root)    
    
    converted_A2 = zssTreeFromTreeLibTree(A2)
    
    zss_distance = zss.compare.distance(A, converted_A2,
                                        zss.simple_tree.Node.get_children,
                                        insert_cost=lambda node: 1,
                                        remove_cost=lambda node: 1,
                                        update_cost=lambda a, b: 0)

    # The trees must be identical (distance = 0) 
    print("Test: %d == 0" %(zss_distance))
    assert(zss_distance == 0)

    del A
    del A2
    del a_root, a_an, a_dn, a_cn, a_bn, a_en

    #
    # Test
    #
    # Again a labvle difference, but this time in the treelib tree.
    A = treelib.Tree()
    a_root = CustomNode(tag="f", identifier="f", treeref=A)
    A.add_node(a_root)  
    a_an = CustomNode(tag="a", identifier="a", treeref=A)
    A.add_node(a_an, parent=a_root)    

    B = treelib.Tree()
    b_root = CustomNode(tag="f", identifier="f", treeref=B)
    B.add_node(b_root)  
    b_an = CustomNode(tag="b", identifier="b", treeref=B) # < Tag difference
    B.add_node(b_an, parent=b_root)    
    
    zss_distance = zss.compare.distance(zssTreeFromTreeLibTree(A),
                                        zssTreeFromTreeLibTree(B),
                                        zss.simple_tree.Node.get_children,
                                        insert_cost=lambda node: 1,
                                        remove_cost=lambda node: 1,
                                        update_cost=lambda a, b: 0)

    # The distance() function should report no cost for the update (distance = 0)
    print("Test: %d == 0" %(zss_distance))
    assert(zss_distance == 0)

    # Test
    A2 = treelib.Tree()
    a_root = CustomNode(tag="f", identifier="f", treeref=A2)
    A2.add_node(a_root)  
    a_an = CustomNode(tag="a", identifier="a", treeref=A2)
    A2.add_node(a_an, parent=a_root)    
    a_dn = CustomNode(tag="d", identifier="d", treeref=A2)
    A2.add_node(a_dn, parent=a_an)    
    a_cn = CustomNode(tag="c", identifier="c", treeref=A2)
    A2.add_node(a_cn, parent=a_an)    
    a_bn = CustomNode(tag="b", identifier="b", treeref=A2)
    A2.add_node(a_bn, parent=a_cn)    
    a_en = CustomNode(tag="e", identifier="e", treeref=A2)
    A2.add_node(a_en, parent=a_root)

    B2     = treelib.Tree()
    b_root = CustomNode(tag="f", identifier="f", treeref=B2)
    B2.add_node(b_root)  
    b_an = CustomNode(tag="a", identifier="a", treeref=B2)
    B2.add_node(b_an, parent=b_root)
    b_dn = CustomNode(tag="d", identifier="d", treeref=A2)
    B2.add_node(b_dn, parent=b_an)       
    b_en = CustomNode(tag="e", identifier="e", treeref=B2)
    B2.add_node(b_en, parent=b_root)
    b_cn = CustomNode(tag="c", identifier="c", treeref=B2)
    B2.add_node(b_cn, parent=b_en)    
    b_bn = CustomNode(tag="b", identifier="b", treeref=B2)
    B2.add_node(b_bn, parent=b_cn)
    g_bn = CustomNode(tag="g", identifier="g", treeref=B2)
    B2.add_node(g_bn, parent=b_bn)

    zss_distance = zss.compare.distance(zssTreeFromTreeLibTree(A2),
                                        zssTreeFromTreeLibTree(B2),
                                        zss.simple_tree.Node.get_children,
                                        insert_cost=lambda node: 1,
                                        remove_cost=lambda node: 1,
                                        update_cost=lambda a, b: 0)
    
    # Total distance is 5: 3 insertions and 2 removals
    print("Test: %d == 5" %(zss_distance))
    assert(zss_distance == 5)

    zss_distance = zss.compare.distance(zssTreeFromTreeLibTree(A2),
                                        zssTreeFromTreeLibTree(B2),
                                        zss.simple_tree.Node.get_children,
                                        insert_cost=lambda node: 1,
                                        remove_cost=lambda node: 0,
                                        update_cost=lambda a, b: 0)
    # Distance is 3: 3 insertions
    print("Test: %d == 3" %(zss_distance))
    assert(zss_distance == 3)

    zss_distance = zss.compare.distance(zssTreeFromTreeLibTree(A2),
                                        zssTreeFromTreeLibTree(B2),
                                        zss.simple_tree.Node.get_children,
                                        insert_cost=lambda node: 0,
                                        remove_cost=lambda node: 1,
                                        update_cost=lambda a, b: 0)
    # Distance is 2: 2 removals
    print("Test: %d == 2" %(zss_distance))
    assert(zss_distance == 2)

    # Test
    A2 = treelib.Tree()
    a_root = CustomNode(tag="f", identifier="f", treeref=A2)
    A2.add_node(a_root)  
    a_an = CustomNode(tag="a", identifier="a", treeref=A2)
    A2.add_node(a_an, parent=a_root)    
    a_dn = CustomNode(tag="d", identifier="d", treeref=A2)
    A2.add_node(a_dn, parent=a_an)    
    a_cn = CustomNode(tag="c", identifier="c", treeref=A2)
    A2.add_node(a_cn, parent=a_an)    
    a_bn = CustomNode(tag="b", identifier="b", treeref=A2)
    A2.add_node(a_bn, parent=a_cn)    
    a_en = CustomNode(tag="e", identifier="e", treeref=A2)
    A2.add_node(a_en, parent=a_root)    

    B2     = treelib.Tree()
    b_root = CustomNode(tag="f", identifier="f", treeref=B2)
    B2.add_node(b_root)  
    b_an = CustomNode(tag="a", identifier="a", treeref=B2)
    B2.add_node(b_an, parent=b_root)    
    b_dn = CustomNode(tag="d", identifier="d", treeref=B2)
    B2.add_node(b_dn, parent=b_an)  
    b_en = CustomNode(tag="e", identifier="e", treeref=B2)
    B2.add_node(b_en, parent=b_root)  
    b_cn = CustomNode(tag="c", identifier="c", treeref=B2)
    B2.add_node(b_cn, parent=b_en)    
    b_bn = CustomNode(tag="b", identifier="b", treeref=B2)
    B2.add_node(b_bn, parent=b_cn)    
    
    # Total distance is 4: 2 insertions and 2 removals
    libtree_distance = zss.compare.distance(a_root, 
                                            b_root,
                                            CustomNode.get_children,
                                            insert_cost=lambda node: 1,
                                            remove_cost=lambda node: 1,
                                            update_cost=lambda a, b: 0)
    
    # Total distance is 4: 2 insertions and 2 removals
    print("Test: %d == 4" %(libtree_distance))
    assert(libtree_distance == 4)
    
    libtree_distance = zss.compare.distance(a_root, 
                                            b_root,
                                            CustomNode.get_children,
                                            insert_cost=lambda node: 1,
                                            remove_cost=lambda node: 0,
                                            update_cost=lambda a, b: 0)
    
    # Distance is 2: 2 removals
    print("Test: %d == 2" %(libtree_distance))
    assert(zss_distance == 2)

    # Total distance is 4: 2 insertions and 2 removals
    libtree_distance = zss.compare.distance(a_root, 
                                            b_root,
                                            CustomNode.get_children,
                                            insert_cost=lambda node: 0,
                                            remove_cost=lambda node: 1,
                                            update_cost=lambda a, b: 0)
    
    # Distance is 2: 2 removal
    print("Test: %d == 2" %(libtree_distance))
    assert(libtree_distance == 2)
    
    
    # Test condensed tree
    A2 = treelib.Tree()
    a_root = CustomNode(tag="f", identifier="f", treeref=A2)
    A2.add_node(a_root)  
    a_an = CustomNode(tag="a", identifier="a", treeref=A2)
    A2.add_node(a_an, parent=a_root)    
    a_dn = CustomNode(tag="d", identifier="d", treeref=A2)
    A2.add_node(a_dn, parent=a_an)    
    a_cn = CustomNode(tag="c", identifier="c", treeref=A2)
    A2.add_node(a_cn, parent=a_an)    
    a_bn = CustomNode(tag="b", identifier="b", treeref=A2)
    A2.add_node(a_bn, parent=a_cn)
    a_gn = CustomNode(tag="g", identifier="g", treeref=A2)
    A2.add_node(a_gn, parent=a_bn)
    a_hn = CustomNode(tag="h", identifier="h", treeref=A2)
    A2.add_node(a_hn, parent=a_bn)
    a_en = CustomNode(tag="e", identifier="e", treeref=A2)
    A2.add_node(a_en, parent=a_root)    
    a_in = CustomNode(tag="i", identifier="i", treeref=A2)
    A2.add_node(a_in, parent=a_en)    
    a_jn = CustomNode(tag="j", identifier="j", treeref=A2)
    A2.add_node(a_jn, parent=a_in)  
    a_kn = CustomNode(tag="k", identifier="k", treeref=A2)
    A2.add_node(a_kn, parent=a_in) 
    a_ln = CustomNode(tag="l", identifier="l", treeref=A2)
    A2.add_node(a_ln, parent=a_kn)  
    a_mn = CustomNode(tag="m", identifier="m", treeref=A2)
    A2.add_node(a_mn, parent=a_kn) 
    a_on = CustomNode(tag="o", identifier="o", treeref=A2)
    A2.add_node(a_on, parent=a_jn) 
    a_pn = CustomNode(tag="p", identifier="p", treeref=A2)
    A2.add_node(a_pn, parent=a_jn) 
    a_qn = CustomNode(tag="q", identifier="q", treeref=A2)
    A2.add_node(a_qn, parent=a_jn) 
    a_rn = CustomNode(tag="r", identifier="r", treeref=A2)
    A2.add_node(a_rn, parent=a_jn) 
    a_sn = CustomNode(tag="s", identifier="s", treeref=A2)
    A2.add_node(a_sn, parent=a_hn) 

    c_A2 = condense_tree(A2)    

    A2.show()
    c_A2.show()
   
    na = c_A2.get_node('a')
    print('a has 3 children')
    assert(na.data == 3)
    nb = c_A2.get_node('b')
    print('b has 3 children')
    assert(nb.data == 3)
    nf = c_A2.get_node('f')
    print('f has 3 children')
    assert(nf.data == 3)
    ni = c_A2.get_node('i')
    print('i has 2 children')
    assert(ni.data == 2)
    nk = c_A2.get_node('k')
    print('k has 2 children')
    assert(nk.data == 2)
    nj = c_A2.get_node('j')
    print('j has 4 children')
    assert(nj.data == 4)

    # Test comparison of condensed trees
    A2 = treelib.Tree()
    a_root = CustomNode(tag="f", identifier="f", treeref=A2)
    A2.add_node(a_root)  
    a_an = CustomNode(tag="a", identifier="a", treeref=A2)
    A2.add_node(a_an, parent=a_root)    
    a_dn = CustomNode(tag="d", identifier="d", treeref=A2)
    A2.add_node(a_dn, parent=a_an)    
    a_cn = CustomNode(tag="c", identifier="c", treeref=A2)
    A2.add_node(a_cn, parent=a_an)    
    a_bn = CustomNode(tag="b", identifier="b", treeref=A2)
    A2.add_node(a_bn, parent=a_cn)
    a_gn = CustomNode(tag="g", identifier="g", treeref=A2)
    A2.add_node(a_gn, parent=a_bn)
    a_hn = CustomNode(tag="h", identifier="h", treeref=A2)
    A2.add_node(a_hn, parent=a_bn)
    a_en = CustomNode(tag="e", identifier="e", treeref=A2)
    A2.add_node(a_en, parent=a_root)    
    a_in = CustomNode(tag="i", identifier="i", treeref=A2)
    A2.add_node(a_in, parent=a_en)    
    a_jn = CustomNode(tag="j", identifier="j", treeref=A2)
    A2.add_node(a_jn, parent=a_in)  
    a_kn = CustomNode(tag="k", identifier="k", treeref=A2)
    A2.add_node(a_kn, parent=a_in) 
    a_ln = CustomNode(tag="l", identifier="l", treeref=A2)
    A2.add_node(a_ln, parent=a_kn)  
    a_mn = CustomNode(tag="m", identifier="m", treeref=A2)
    A2.add_node(a_mn, parent=a_kn) 
    a_on = CustomNode(tag="o", identifier="o", treeref=A2)
    A2.add_node(a_on, parent=a_jn) 
    a_pn = CustomNode(tag="p", identifier="p", treeref=A2)
    A2.add_node(a_pn, parent=a_jn) 
    a_qn = CustomNode(tag="q", identifier="q", treeref=A2)
    A2.add_node(a_qn, parent=a_jn) 
    a_rn = CustomNode(tag="r", identifier="r", treeref=A2)
    A2.add_node(a_rn, parent=a_jn) 
    a_sn = CustomNode(tag="s", identifier="s", treeref=A2)
    A2.add_node(a_sn, parent=a_hn) 

    # Test condensed tree
    B2 = treelib.Tree()
    b_root = CustomNode(tag="f", identifier="f", treeref=B2)
    B2.add_node(b_root)  
    b_an = CustomNode(tag="a", identifier="a", treeref=B2)
    B2.add_node(b_an, parent=b_root)    
    b_dn = CustomNode(tag="d", identifier="d", treeref=B2)
    B2.add_node(b_dn, parent=b_an)    
    b_cn = CustomNode(tag="c", identifier="c", treeref=B2)
    B2.add_node(b_cn, parent=b_an)    
    b_bn = CustomNode(tag="b", identifier="b", treeref=B2)
    B2.add_node(b_bn, parent=b_cn)
    b_gn = CustomNode(tag="g", identifier="g", treeref=B2)
    B2.add_node(b_gn, parent=b_bn)
    b_hn = CustomNode(tag="h", identifier="h", treeref=B2)
    B2.add_node(b_hn, parent=b_bn)
    b_en = CustomNode(tag="e", identifier="e", treeref=B2)
    B2.add_node(b_en, parent=b_root)    
    b_in = CustomNode(tag="i", identifier="i", treeref=B2)
    B2.add_node(b_in, parent=b_dn)    
    b_jn = CustomNode(tag="j", identifier="j", treeref=B2)
    B2.add_node(b_jn, parent=b_in)  
    b_kn = CustomNode(tag="k", identifier="k", treeref=B2)
    B2.add_node(b_kn, parent=b_in) 
    b_ln = CustomNode(tag="l", identifier="l", treeref=B2)
    B2.add_node(b_ln, parent=b_kn)  
    b_mn = CustomNode(tag="m", identifier="m", treeref=B2)
    B2.add_node(b_mn, parent=b_kn) 
    b_on = CustomNode(tag="o", identifier="o", treeref=B2)
    B2.add_node(b_on, parent=b_jn) 
    b_pn = CustomNode(tag="p", identifier="p", treeref=B2)
    B2.add_node(b_pn, parent=b_jn) 
    b_qn = CustomNode(tag="q", identifier="q", treeref=B2)
    B2.add_node(b_qn, parent=b_jn) 
    b_rn = CustomNode(tag="r", identifier="r", treeref=B2)
    B2.add_node(b_rn, parent=b_jn) 
    b_sn = CustomNode(tag="s", identifier="s", treeref=B2)
    B2.add_node(b_sn, parent=b_hn) 

    A2.show()
    B2.show()
    c_A2 = condense_tree(A2)
    c_B2 = condense_tree(B2)
    c_A2.show()
    c_B2.show()
    

    # Total distance is 14: 7 insertions and 7 removals
    libtree_distance = zss.compare.distance(A2.get_node(A2.root),
                                            B2.get_node(B2.root),
                                            CustomNode.get_children,
                                            insert_cost=lambda node: 1,
                                            remove_cost=lambda node: 1,
                                            update_cost=lambda a, b: 0)
    print("Test: %d == 14" %(libtree_distance))
    assert(libtree_distance == 14)
    
    # Total distance is 2: 1 insertion and 1 removal
    libtree_distance = zss.compare.distance(c_A2.get_node(c_A2.root),
                                            c_B2.get_node(c_B2.root),
                                            CustomNode.get_children,
                                            insert_cost=lambda node: 1,
                                            remove_cost=lambda node: 1,
                                            update_cost=lambda a, b: 0)
    print("Test: %d == 2" %(libtree_distance))
    assert(libtree_distance == 2)

    # Total distance is 3: 1 insertion, 1 removal and 1 update
    libtree_distance = zss.compare.distance(c_A2.get_node(c_A2.root),
                                            c_B2.get_node(c_B2.root),
                                            CustomNode.get_children,
                                            insert_cost=lambda node: 1,
                                            remove_cost=lambda node: 1,
                                            update_cost=compare_num_children)

    print("Test: %d == 3" %(libtree_distance))
    assert(libtree_distance == 3)


if __name__ == '__main__':
    """
    Entry point to run the XTension from the command line.
    """

    # Run a series of tests.
    #test()

    # Connect an existing and ready Imaris instance. Make sure that
    # Imaris is open, and all spots with their tracks and lineage
    # traces exist in the surpass scene. A dialog will pop up to
    # allow you to pick the ones to compare.
    conn = pIce(0)

    # Calculate the distance: make sure to set the
    # proper parameters at the top of this file!
    bsseLineageTracerTotalDistanceBetweenLineageTrees(conn)
