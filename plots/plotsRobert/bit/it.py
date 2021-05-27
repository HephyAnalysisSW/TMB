# Standard imports
import numpy as np
import copy
import cProfile
import operator 

from Analysis.Tools.helpers import chunk

class Node:
    def __init__( self, dataset, max_depth, min_size, depth=0):

        ## basic BDT configuration
        self.max_depth  = max_depth
        self.min_size   = min_size

        # data set
        self.dataset    = dataset
        self.size       = float(len(self.dataset))
        self.depth      = depth

        self.split(depth=depth)

    # compute the gini impurity from a set of booleans defining the 'left' box and (by negation) the 'right' box
    def gini( self, value, left_indices):
        ''' Calculate Gini index for split dataset
        '''
        gini = 0
        for group in [left_indices, ~left_indices]:
            group_size = float(np.count_nonzero(group))
            if group_size==0:
                continue
            score = 0.

            content = np.bincount(self.dataset[group][:,-1].astype('int64') )
            for class_val in self.classes:
                p = content[class_val]/group_size if class_val<len(content) else 0
                score += p*p
            gini += (1-score)*(group_size/self.size)
        return gini


    def get_split( self ):
        ''' determine where to split the dataset
        '''

        # number of classes in this node
        self.classes = list(set(self.dataset[:,-1].astype('int64')))

        # loop over the features ... assume the dataset consists of rows with [x1, x2, ...., xN, y]

        self.split_i_feature, self.split_value, self.split_score, self.split_left_group = float('inf'), float('inf'), float('inf'), None

        # loop over features, i.e., the columns in the dataset excluding the last one (which is the training target)
        for i_feature in range(len(self.dataset[0])-1):
            feature_values = self.dataset[:,i_feature]
            # get column & loop over all values
            for value in feature_values:
                #print (i_feature, value)
                left_indices = feature_values<value
                gini = self.gini( value, left_indices )
                if gini < self.split_score:
                    #print ("found better!", i_feature, value, "old score", self.split_score, "gini", gini)
                    self.split_i_feature = i_feature
                    self.split_value     = value
                    self.split_score     = gini
                    self.split_left_group= left_indices

        #print ("final:get_split", self.split_i_feature, self.split_value)
        #return {'split_i_feature':split_i_feature, 'split_value':split_value, 'split_score':split_score, 'split_left_group':split_left_group}

    def get_result(self, group):
        ''' Return the result. 
        In this example it is merely the maximum of the occurences of all classes, i.e, if class '0' have the
        relative majority, the result is '0'.
        '''
        outcomes = self.dataset[:,-1][group].astype('int64')
        counts = np.bincount(outcomes)
        #print ("get_result counts", counts, "argmax", np.argmax(counts) )
        return np.argmax(counts)

    # Create child splits for a node or make terminal
    def split(self, depth):

        # Find the best split
        self.get_split()

        # check for a 'no' split
        if not any(self.split_left_group) or all(self.split_left_group):
            #self.left = self.right = self.get_result(np.s_[:])
            #print ("Choice1", depth, self.get_result(np.s_[:]) )
            # If one of the groups after the split is empty, this wasn't a split. Return the result of this box, i.e., unsplit. 
            self.left = self.right = ResultNode(self.get_result(np.s_[:]))
            return
        # check for max depth
        if  max_depth <= depth+1: # Jason Brownlee starts counting depth at 1, we start counting at 0, hence the +1
            #print ("Choice2", depth, self.get_result(self.split_left_group), self.get_result(~self.split_left_group) )
            # The split was good, but we stop splitting further. Put the result of the split in the left/right boxes.
            self.left, self.right = ResultNode(self.get_result(self.split_left_group)), ResultNode(self.get_result(~self.split_left_group))
            return
        # process left child
        if np.count_nonzero(self.split_left_group) <= min_size:
            #print ("Choice3", depth, self.get_result(self.split_left_group) )
            # Too few events in the left box. We stop.
            self.left             = ResultNode(self.get_result(self.split_left_group))
        else:
            #print ("Choice4", depth )
            # Continue splitting left box. 
            self.left             = Node(self.dataset[self.split_left_group], max_depth=self.max_depth, min_size=self.min_size, depth=self.depth+1 )
        # process right child
        if np.count_nonzero(~self.split_left_group) <= min_size:
            #print ("Choice5", depth, self.get_result(~self.split_left_group) )
            # Too few events in the right box. We stop.
            self.right            = ResultNode(self.get_result(~self.split_left_group))
        else:
            #print ("Choice6", depth  )
            # Continue splitting right box. 
            self.right            = Node(self.dataset[~self.split_left_group], max_depth=self.max_depth, min_size=self.min_size, depth=self.depth+1 )

    # Prediction    
    def predict( self, row ):
        ''' obtain the result by recursively descending down the tree
        '''
        node = self.left if row[self.split_i_feature]<self.split_value else self.right
        if isinstance(node, ResultNode):
            return node.return_value
        else:
            return node.predict(row)

    # Print a decision tree
    def print_tree(self, depth=0):
        print('%s[X%d < %.3f]' % ((self.depth*' ', self.split_i_feature, self.split_value)))
        for node in [self.left, self.right]:
            node.print_tree(depth+1)

class ResultNode:
    ''' Simple helper class to store result value.
    '''
    def __init__( self, return_value ):
        self.return_value   = return_value
    def print_tree(self, depth=0):
        print('%s[%s]' % (((depth)*' ', self.return_value)))

from csv import reader
# Load a CSV file
def load_csv(filename = 'data_banknote_authentication.csv'):
    file = open(filename, "rt")
    lines = reader(file)
    dataset = list(lines)
    for row in dataset:
        for i_column, column in enumerate(row):
            row[i_column] = float(column.strip())

    return np.array(dataset)

max_depth = 5
min_size = 10

import uproot
import awkward
import numpy as np
import pandas as pd

input_file = "/eos/vbc/user/robert.schoefbeck/TMB/bit/MVA-training/ttG_WG_small/WGToLNu_fast/WGToLNu_fast.root"
upfile     = uproot.open( input_file )
tree       = upfile["Events"]
n_events   = len( upfile["Events"] )
n_events   = max(10000, n_events)
entrystart, entrystop = 0, n_events 

# Load features
features   = [ "mva_photon_pt", "mva_photon_eta", "mva_photonJetdR", "mva_photonLepdR", "mva_mT" ]
df         = tree.pandas.df(branches = features, entrystart=entrystart, entrystop=entrystop)
x          = df.values

# Load weights
from Analysis.Tools.WeightInfo import WeightInfo
w = WeightInfo("/eos/vbc/user/robert.schoefbeck/gridpacks/v6/WGToLNu_reweight_card.pkl")

# Load all weights and reshape the array according to ndof from weightInfo
weights    = tree.pandas.df(branches = ["p_C"], entrystart=entrystart, entrystop=entrystop).values.reshape((-1,w.nid))

assert len(x)==len(weights), "Need equal length for weights and features."

# shuffle!

