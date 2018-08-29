
import random
import os
import h5py
import zarr
import sys
import pandas as pd
import daiquiri
#import bsddb3
import time
import scipy
import pickle
import collections
import itertools
import operator
import tqdm
import shutil
import pprint
import numpy as np
import json

import matplotlib as mp
# Force matplotlib to not use any Xwindows backend.
mp.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import tsinfer
import msprime


old_simplify = msprime.TreeSequence.simplify
def new_simplify(ts, *args, **kwargs):
    tmp_tables = ts.dump_tables()
    tmp_tables.nodes.clear()
    tmp_tables.individuals.clear()
    for node in ts.nodes():
        if node.is_sample():
            try:
                metadata = ts.individual(node.individual).metadata
            except:
                metadata = None
        else:
            metadata = None
        tmp_tables.nodes.add_row(
            flags=node.flags, time=node.time, metadata=metadata)
    tmp_tables.nodes.set_columns(
       flags=tmp_tables.nodes.flags,
       time=tmp_tables.nodes.time,
       metadata=tmp_tables.nodes.metadata,
       metadata_offset=tmp_tables.nodes.metadata_offset)
    return old_simplify(tmp_tables.tree_sequence(), *args, **kwargs)

msprime.TreeSequence.simplify = new_simplify


def make_errors(v, p):
    """
    For each sample an error occurs with probability p. Errors are generated by
    sampling values from the stationary distribution, that is, if we have an
    allele frequency of f, a 1 is emitted with probability f and a
    0 with probability 1 - f. Thus, there is a possibility that an 'error'
    will in fact result in the same value.
    """
    w = np.copy(v)
    if p > 0:
        m = v.shape[0]
        frequency = np.sum(v) / m
        # Randomly choose samples with probability p
        samples = np.where(np.random.random(m) < p)[0]
        # Generate observations from the stationary distribution.
        errors = (np.random.random(samples.shape[0]) < frequency).astype(int)
        w[samples] = errors
    return w


def generate_samples(ts, error_p):
    """
    Returns samples with a bits flipped with a specified probability.

    Rejects any variants that result in a fixed column.
    """
    S = np.zeros((ts.sample_size, ts.num_mutations), dtype=np.int8)
    for variant in ts.variants():
        done = False
        # Reject any columns that have no 1s or no zeros
        while not done:
            S[:, variant.index] = make_errors(variant.genotypes, error_p)
            s = np.sum(S[:, variant.index])
            done = 0 < s < ts.sample_size
    return S.T

class RangeData:
    """
    A consecutive range of edges belonging to the same parent
    There should be n positions and n-1 parents
    """
    def __init__(self):
        self.pos_array = []
        self.parent_array = []

    @property
    def leftmost(self):
        return self.pos_array[0]

    @property
    def rightmost(self):
        return self.pos_array[-1]

    def left_parent(self, pos):
        if pos > 0:
            return self.parent_array[pos-1]
        else:
            return None

    def right_parent(self, pos):
        if pos < len(self.parent_array):
            return self.parent_array[pos]
        else:
            return None

    def contains(self, position):
        """
        Does this range contain this position. We deliberately allow both ends open
        as we will sometimes want to obtain the leftmost pos for a breakpoint which 
        occurs exactly at the rightmost position
        """
        return self.leftmost <= position <= self.rightmost

    def delete_intermediate_position(self, pos_idx, new_left_and_right_parent):
        assert pos_idx > 0 # must be intermediate
        assert pos_idx < len(self.parent_array), (pos_idx, len(self.parent_array))  # must be intermediate
        del self.pos_array[pos_idx]
        del self.parent_array[pos_idx]
        self.parent_array[pos_idx-1] = new_left_and_right_parent

SRB_data = collections.namedtuple("SRB_data", "pos left_parent right_parent")

class ChildData:
    """
    Store a list of consecutive ranges, and a dict mapping positions into those ranges.
    We need to:
    a) get the rightmost and leftmost position in this range
    b) delete a position (and replace the right and left parents with a new one)
    c) get the current left and right parent for a specific position
    d) insert an edge, 
    """
    def __init__(self, left=None, middle=None, right=None, parent1=None, parent2=None):
        """
        Can either be constructed with no args, to give a blank child, or with a pair 
        of consecutive edges
        """
        self.temp_edges = []
        self.ranges = [] # a list of RangeData objects
        if left!=None:
            self.ranges += [RangeData()]
            self.ranges[0].pos_array += [left, middle, right]
            self.ranges[0].parent_array += [parent1, parent2]
        
        
    def delete_position(self, position, new_left_and_right_parent, verbosity=0):
        """
        This should not allow us to delete positions at either end of a range,
        so we should not be able to change the length of ranges
        """
        assert len(self.ranges)
        for r in self.ranges:
            try:
                del_pos = r.pos_array.index(position)
                if verbosity:
                    printed_range = ["{:.3f}".format(pos) + ("" if par is None else " [{}]".format(par))
                        for pos, par in zip(r.pos_array, r.parent_array+[None])]
                    print("Deleting position #{} (@{}) (parents no longer {} and {}, but {}) from range\n {}".format(
                        del_pos, r.pos_array[del_pos], 
                        r.parent_array[del_pos-1], r.parent_array[del_pos], new_left_and_right_parent, 
                        " ".join(printed_range)))
                r.delete_intermediate_position(del_pos, new_left_and_right_parent)
                return 1
            except ValueError:
                continue
        assert False, "No such position in any range"
                
    def get_SRB(self, position):
        assert len(self.ranges)
        for r in self.ranges:
            try:
                pos = r.pos_array.index(position)
                return SRB_data(position, r.left_parent(pos), r.right_parent(pos))
            except ValueError:
                continue
        assert False, "No such position in any range"

    def rightmost_from(self, position):
        assert len(self.ranges)
        for r in self.ranges:
            if r.contains(position):
               return r.rightmost
        assert False, "Position {} not in any range".format(position)

    def leftmost_from(self, position):
        assert len(self.ranges)
        for r in self.ranges:
            if r.contains(position):
                return r.leftmost
        assert False, "Position {} not in any range: {}".format(position,
            [r.pos_array for r in self.ranges])

        
    def add_edge(self, left, right, parent_id):
        """
        Add an edge, potentially out of order. Call make_ranges() before using
        """
        assert len(self.ranges) == 0 #should only be used to initialize the edges
        self.temp_edges.append((left, right, parent_id))

    def make_ranges(self):
        """
        Once we have added a set of edges, this will determine how to cut them into 
        contiguous ranges
        """
        self.temp_edges.sort(key=operator.itemgetter(0))
        prev_right=None
        for left, right, parent_id in self.temp_edges:
            if left!=prev_right:
                #make a new range
                self.ranges.append(RangeData())
                self.ranges[-1].pos_array.append(left)
            self.ranges[-1].pos_array.append(right)
            self.ranges[-1].parent_array.append(parent_id)
            prev_right = right
        for r in self.ranges:
            assert len(r.pos_array)-1 == len(r.parent_array), \
                (len(r.pos_array), len(r.parent_array))
        self.temp_edges = [] #clear memory

def identify_SRBs(ts):
    breakpoints = collections.defaultdict(set)
    for (left, right), edges_out, edges_in in ts.edge_diffs():
        if len(edges_out) and len(edges_in):
            child_data = collections.defaultdict(list)
            for edge in edges_out:
                child_data[edge.child].append(edge.parent)
            for edge in edges_in:
                child_data[edge.child].append(edge.parent)
            for c, parents in child_data.items():
                assert 0 < len(parents) < 3
                if len(parents) == 2:
                    breakpoints[SRB_data(left, parents[0],parents[1])].add(c)
    # shared breakpoints have at least 2 children
    return {k:v for k,v in breakpoints.items() if len(v) > 1}


def SRB_replace_edges(new_id, time, SRB, SRB_children, child_to_parent, verbosity=0):
    """
    Insert 2 edges pointing to a new id into the child_to_parent list, and 
    delete the existing SRB, replacing it with a single edge pointing to the
    newly created edges.
    Return a tuple of the number of edges deleted and created
    """
    # find leftmost span of l parent and rightmost of r parent
    # this assumes contiguous edges on a parent    
    l_parent_lft = child_to_parent[SRB.left_parent].leftmost_from(SRB.pos)
    r_parent_rgt = child_to_parent[SRB.right_parent].rightmost_from(SRB.pos)
    #insert as replacement
    assert new_id != SRB.left_parent
    assert new_id != SRB.right_parent
    child_to_parent[new_id] = ChildData(
        l_parent_lft, SRB.pos, r_parent_rgt, SRB.left_parent, SRB.right_parent)
    edges_inserted = 2
    if verbosity:
        print("Inserted edge from child", new_id, "to parent", SRB.left_parent, 
            ". Location: [{}, {})".format(l_parent_lft, SRB.pos))
        print("Inserted edge from child", new_id, "to parent", SRB.right_parent, 
            ". Location: [{}, {})".format(SRB.pos, r_parent_rgt))
    
    # relabel existing children to point to this new parent either side of the
    # breakpoint
    edges_deleted = 0
    for child in SRB_children:
        if verbosity:
            print("Deleting at position {} by merging two edges in child {}".format(SRB.pos, child))
        edges_deleted += child_to_parent[child].delete_position(SRB.pos, new_id, verbosity)
    return edges_inserted, edges_deleted

def tsinfer_dev(
        n, L, seed, num_threads=1, recombination_rate=1e-8,
        error_rate=0, engine="C", log_level="WARNING",
        debug=True, progress=False, path_compression=True):

    np.random.seed(seed)
    random.seed(seed)
    L_megabases = int(L * 10**6)

    # daiquiri.setup(level=log_level)

    ts = msprime.simulate(
            n, Ne=10**4, length=L_megabases,
            recombination_rate=recombination_rate, mutation_rate=1e-8,
            random_seed=seed)
    if debug:
        print("num_sites = ", ts.num_sites)
    assert ts.num_sites > 0

    use_built_in_path_compression = True

    sample_data = tsinfer.SampleData.from_tree_sequence(ts)

    ancestor_data = tsinfer.generate_ancestors(
        sample_data, engine=engine, num_threads=num_threads)
    ancestors_ts = tsinfer.match_ancestors(
        sample_data, ancestor_data, engine=engine, 
        path_compression=use_built_in_path_compression, extended_checks=True)

    # Do not simplify the ancestors, since (1) we will not be able to map parents
    #  back to ancestors in the original ancestors_ts, and (2) we may lose where in the 
    #  timeslices to put the new ancestors. OTOH, we may miss some SRBs, if parents are 
    #  different between children
    # Additionally, we can't use path compression here, as it will create more ancestors
    full_ts_pc = tsinfer.match_samples(
        sample_data, ancestors_ts, engine=engine, simplify=True,
        path_compression=True, extended_checks=True)


    full_ts = tsinfer.match_samples(
        sample_data, ancestors_ts, engine=engine, simplify=False,
        path_compression=False, extended_checks=True)

    SRBs = identify_SRBs(full_ts)
    print("\nSRB count at start")
    # print how many SRBs we have
    ct = np.array([len(bp) for bp in SRBs.values()], dtype=np.int)
    edgecount = full_ts.num_edges, full_ts.simplify().num_edges
    print(
        ", ".join(["{}:{}".format(i,n) for i,n in enumerate(np.bincount(ct)) if i>1 and n>0])
        + "\n -- N edges in initial ts = {} ({} simplified, ".format(*edgecount)
        + "{} with path compression on samples)".format(full_ts_pc.num_edges))
    # sort so that we hit the most frequent SRBs first

    sorted_SRB_keys = sorted(list(SRBs.keys()), key=lambda x:len(SRBs[x]), reverse=True)
    
    tables = full_ts.dump_tables()
    samples = set()
    # make a new structure indexed by child rather than by parent
    # we need to keep track of samples too, as the parents of samples need to be
    # passed up
    child_to_parent = collections.defaultdict(ChildData)
    internal_node_times = {}
    for e in tables.edges:
        child_to_parent[e.child].add_edge(e.left, e.right, e.parent)
    for arr in child_to_parent.values():
        arr.make_ranges()
    for i, n in enumerate(tables.nodes):
        if n.flags & msprime.NODE_IS_SAMPLE:
            samples.add(i)
        else:
            internal_node_times[i] = n.time
    print("Samples", samples)

    edgelist = list(ancestors_ts.edges())
    edgelist.sort(key=operator.attrgetter('child', 'left'))
    print("\n".join([str(e) for e in edgelist]))
    #print("\n")

    #here we revert to the ancestors TS so that we can
    tables = ancestors_ts.dump_tables()
    #internal_node_times = {i:n.time for i,n in enumerate(tables.nodes)}

    #for t in ancestors_ts.trees():
    #    #print("== {} - {} ===".format(*t.interval))
    #    print(t.draw(format="unicode"))
    #allocate a new ancestor for every shared breakpoint site        
    delta=0.0000000001 #must be smaller than the delta used in path compression
    edges_differences = np.zeros(2, dtype=np.int) #count inserted, deleted
    for it, SRB in enumerate(sorted_SRB_keys):
        #print("== Collapsing SRB at position {} ==".format(SRB))
        new_SRBs = collections.defaultdict(list)
        for child in SRBs[SRB]:
            # Separate into SRBs with matching L & R parents 
            # Note that the true L or R parent may not be that described by SRB.left_parent or
            # SRB.right_parent, since the parents may have been replaced on previous
            # passes of the algorithm, so we need to recalculate the true L & R parents, and 
            # group them if necessary
            new_SRB = child_to_parent[child].get_SRB(SRB.pos)
            new_SRBs[new_SRB].append(child)
        if len(new_SRBs) > 1:
            print("What was a single SRB ({} with children {}) has been split into more than one: {}".format(
                SRB, SRBs[SRB], ", ".join(["{} with children {}".format(n, c) for n, c in new_SRBs.items()])))
            assert False
        for new_SRB, children in new_SRBs.items():
            if len(children) > 1:
                try:
                    youngest_parent_time = min(
                        internal_node_times[new_SRB.left_parent],
                        internal_node_times[new_SRB.right_parent])
                except KeyError:
                    print("Data for SRB", new_SRB, "with children", children, internal_node_times)
                    raise
                new_time = youngest_parent_time-delta
                #NB: to avoid collision with existing sample nodes, we should refer to these
                #newly created nodes by a *negative* number
                new_id = -tables.nodes.add_row(time=new_time, flags=tsinfer.SYNTHETIC_NODE_BIT)
                internal_node_times[new_id] = new_time
                #print("Inserted new node", new_id, "@", new_time)
                edges_differences += SRB_replace_edges(
                    new_id, youngest_parent_time, new_SRB, children, child_to_parent, 1)
    # remake the ancestors table, but don't use the samples: we will match them up later
    tables.edges.clear()
    for c, data in child_to_parent.items():
        if c in samples:
            #print("omitting edges from child {}".format(c))
            pass
        else:
            #print("Saving edges from child {}".format(c))
            for r in data.ranges:
                for i, parent_id in enumerate(r.parent_array):
                    try:
                        assert parent_id != c
                        assert internal_node_times[parent_id] > internal_node_times[c]
                    except:
                        print(
                            "Parent {}@{}, child {}@{}".format(parent_id,
                                internal_node_times[parent_id], c, internal_node_times[c]))
                        raise
                    tables.edges.add_row(
                        left=r.pos_array[i], right=r.pos_array[i+1], 
                        parent=abs(parent_id), child=abs(c))
    
    tables.sort()


    ancestors_ts = tables.tree_sequence()

    #print("\n".join([str(e) for e in ancestors_ts.nodes()]))
    edgelist = list(ancestors_ts.edges())
    edgelist.sort(key=operator.attrgetter('child', 'left'))
    #print("\n".join([str(e) for e in edgelist]))
    e_iter = iter(edgelist)
    #for c in sorted(child_to_parent.keys(), key=abs):
    #    if c not in samples:
    #        for r in child_to_parent[c].ranges:
    #            for p, l, r in zip(r.parent_array, r.pos_array[:-1], r.pos_array[1:]):
    #                print("[left={:.3f}, right={:.3f}, parent={}, child={}]".format(
    #                    l, r, p, c))
    #                e = next(e_iter)
    #                #assert e.left == l and e.right==r and e.parent==p and e.child==c

    #print(ancestors_ts.first().draw(format="unicode"))

    print("new ancestors tree sequence made")    
    full_ts = tsinfer.match_samples(
        sample_data, ancestors_ts, simplify=False,
        #engine="P",
        engine="C",
        path_compression=use_built_in_path_compression)


    

    SRBs = identify_SRBs(full_ts)
    print("SRB count at end")
    # print how many SRBs we have
    ct = np.array([len(bp) for bp in SRBs.values()], dtype=np.int)
    post_edgecount = full_ts.num_edges, full_ts.simplify().num_edges
    print(
        ", ".join(["{}:{}".format(i,n) for i,n in enumerate(np.bincount(ct)) if i>1])
        + " -- N edges in final ts = {} ({} simplified)".format(*post_edgecount))

    for c, pre, post in zip(["without", "with"], edgecount, post_edgecount):
        print("Reduction in edges {} simplification: {:.3f} %".format(
            c, (pre-post)/pre*100))

    for node in ts.nodes():
        if tsinfer.is_synthetic(node.flags):
            print("Synthetic node", node.id, node.time)
            parent_edges = [edge for edge in ts.edges() if edge.parent == node.id]
            child_edges = [edge for edge in ts.edges() if edge.child == node.id]
            child_edges.sort(key=lambda e: e.left)
            print("parent edges")
            for edge in parent_edges:
                print("\t", edge)
            print("child edges")
            for edge in child_edges:
                print("\t", edge)

def subset_sites(ts, position):
    """
    Return a copy of the specified tree sequence with sites reduced to those
    with positions in the specified list.
    """
    tables = ts.dump_tables()
    lookup = frozenset(position)
    tables.sites.clear()
    tables.mutations.clear()
    for site in ts.sites():
        if site.position in lookup:
            site_id = tables.sites.add_row(
                site.position, ancestral_state=site.ancestral_state,
                metadata=site.metadata)
            for mutation in site.mutations:
                tables.mutations.add_row(
                    site_id, node=mutation.node, parent=mutation.parent,
                    derived_state=mutation.derived_state,
                    metadata=mutation.metadata)
    return tables.tree_sequence()

def minimise(ts):
    tables = ts.dump_tables()

    out_map = {}
    in_map = {}
    first_site = 0
    for (_, edges_out, edges_in), tree in zip(ts.edge_diffs(), ts.trees()):
        for edge in edges_out:
            out_map[edge.child] = edge
        for edge in edges_in:
            in_map[edge.child] = edge
        if tree.num_sites > 0:
            sites = list(tree.sites())
            if first_site:
                x = 0
                first_site = False
            else:
                x = sites[0].position
            print("X = ", x)
            for edge in out_map.values():
                print("FLUSH", edge)
            for edge in in_map.values():
                print("INSER", edge)

            # # Flush the edge buffer.
            # for left, parent, child in edge_buffer:
            #     tables.edges.add_row(left, x, parent, child)
            # # Add edges for each node in the tree.
            # edge_buffer.clear()
            # for root in tree.roots:
            #     for u in tree.nodes(root):
            #         if u != root:
            #             edge_buffer.append((x, tree.parent(u), u))

    # position = np.hstack([[0], tables.sites.position, [ts.sequence_length]])
    # position = tables.sites.position
    # edges = []
    # print(position)
    # tables.edges.clear()
    # for edge in ts.edges():
    #     left = np.searchsorted(position, edge.left)
    #     right = np.searchsorted(position, edge.right)

    #     print(edge, left, right)
    #     # if right - left > 1:
    #         # print("KEEP:", edge, left, right)
    #         # tables.edges.add_row(
    #         #     position[left], position[right], edge.parent, edge.child)
    #         # print("added", tables.edges[-1])
    #     # else:
    #         # print("SKIP:", edge, left, right)

    # ts = tables.tree_sequence()
    # for tree in ts.trees():
    #     print("TREE:", tree.interval)
    #     print(tree.draw(format="unicode"))





def minimise_dev():
    ts = msprime.simulate(5, mutation_rate=1, recombination_rate=2, random_seed=3)
    # ts = msprime.load(sys.argv[1])

    position = ts.tables.sites.position[::2]
    subset_ts = subset_sites(ts, position)
    print("Got subset")

    ts_new = tsinfer.minimise(subset_ts)
    for tree in ts_new.trees():
        print("TREE:", tree.interval)
        print(tree.draw(format="unicode"))
    # print(ts_new.tables)
    print("done")
    other = minimise(subset_ts)


def run_build():

    sample_data = tsinfer.load(sys.argv[1])
    ad = tsinfer.generate_ancestors(sample_data)
    print(ad)


if __name__ == "__main__":

    # run_build()

    # np.set_printoptions(linewidth=20000)
    # np.set_printoptions(threshold=20000000)

    # tutorial_samples()

    # build_profile_inputs(10, 10)
    # build_profile_inputs(100, 10)
    # build_profile_inputs(1000, 100)
    # build_profile_inputs(10**4, 100)
    # build_profile_inputs(10**5, 100)

    for j in range(30,31):
        print(j)
        tsinfer_dev(4, 0.03, seed=j, num_threads=0, engine="P", recombination_rate=1e-8)
    # copy_1kg()
    #tsinfer_dev(4, 0.05, seed=123, num_threads=0, engine="C", recombination_rate=1e-8)
    #tsinfer_dev(10, 5, seed=456, num_threads=0, engine="C", recombination_rate=1e-8)
    #tsinfer_dev(10, 5, seed=689, num_threads=0, engine="C", recombination_rate=1e-8)
    #tsinfer_dev(10, 5, seed=101112, num_threads=0, engine="C", recombination_rate=1e-8)

    # minimise_dev()

#     for seed in range(1, 10000):
#         print(seed)
#         # tsinfer_dev(40, 2.5, seed=seed, num_threads=1, genotype_quality=1e-3, engine="C")
