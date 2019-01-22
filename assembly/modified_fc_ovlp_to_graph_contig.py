#!/usr/bin/env python
#coding: utf-8

#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
# Copyright (c) 2019, Yoshihiko Suzuki
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$

#from pbcore.io import FastaReader
import networkx as nx
import os
import shlex
import sys
import subprocess
import re
from collections import defaultdict

DEBUG_LOG_LEVEL = 0

class SGNode(object):
    """
    class representing a node in the string graph
    """
    def __init__(self, node_name):
        self.name = node_name
        self.out_edges = []
        self.in_edges = []
    def add_out_edge(self, out_edge):
        self.out_edges.append(out_edge)
    def add_in_edge(self, in_edge):
        self.in_edges.append(in_edge)

class SGEdge(object):
    """
    class representing an edge in the string graph
    """
    def __init__(self, in_node, out_node):
        self.in_node = in_node
        self.out_node = out_node
        self.attr = {}
    def set_attribute(self, attr, value):
        self.attr[attr] = value

def reverse_end( node_id ):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end

class StringGraph(object):
    """
    class representing the string graph
    """
    def __init__(self):
        self.nodes = {}
        self.edges = {}
        self.n_mark = {}
        self.e_reduce = {}
        self.repeat_overlap = {}
        
    def add_node(self, node_name):
        """ 
        add a node into the graph by given a node name
        """
        if node_name not in self.nodes:
            self.nodes[node_name] = SGNode(node_name)
    
    def add_edge(self, in_node_name, out_node_name, **attributes):
        """ 
        add an edge into the graph by given a pair of nodes
        """
        if (in_node_name, out_node_name) not in self.edges:
        
            self.add_node(in_node_name)
            self.add_node(out_node_name)
            in_node = self.nodes[in_node_name]
            out_node = self.nodes[out_node_name]    
            
            edge = SGEdge(in_node, out_node)
            self.edges[ (in_node_name, out_node_name) ] = edge
            in_node.add_out_edge(edge)
            out_node.add_in_edge(edge)
        edge =  self.edges[ (in_node_name, out_node_name) ]
        for k, v in attributes.items():
            edge.attr[k] = v

    def init_reduce_dict(self):
        for e in self.edges:
            self.e_reduce[e] = False

    def mark_chimer_edges(self):

        for n_name in self.nodes:
            n = self.nodes[n_name]
            
            out_nodes = set( [ e.out_node for e in n.out_edges ] )
            in_nodes = [e.in_node for e in n.in_edges ] 
            is_chimer = True
            for in_node in in_nodes:
                for v in [e.out_node for e in in_node.out_edges]:
                    if v in out_nodes:
                        is_chimer = False
                        break

            if is_chimer == True:
                for e in n.out_edges:
                    v, w =  e.in_node.name, e.out_node.name
                    self.e_reduce[ (v, w) ] = True
                for e in n.in_edges:
                    v, w =  e.in_node.name, e.out_node.name
                    self.e_reduce[ (v, w) ] = True


            # need to remove the node from the graph rather than just mark the edges are "reduced"?


    def mark_spur_edge(self):

        removed_edges = set()
        for  v in self.nodes:
            if len(self.nodes[v].out_edges) > 1:
                for out_edge in self.nodes[v].out_edges:
                    w = out_edge.out_node.name
                    
                    if len(self.nodes[w].out_edges) == 0 and self.e_reduce[ (v, w) ] != True:
                        self.e_reduce[(v, w)] = True
                        removed_edges.add( (v, w) )
                        v2, w2 = reverse_end(w), reverse_end(v)
                        self.e_reduce[(v2, w2)] = True
                        removed_edges.add( (v2, w2) )

            if len(self.nodes[v].in_edges) > 1:
                for in_edge in self.nodes[v].in_edges:
                    w = in_edge.in_node.name
                    if len(self.nodes[w].in_edges) == 0 and self.e_reduce[ (w, v) ] != True:
                        self.e_reduce[(w, v)] = True
                        removed_edges.add( (w, v) )
                        v2, w2 = reverse_end(w), reverse_end(v)
                        self.e_reduce[(w2, v2)] = True
                        removed_edges.add( (w2, v2) )
        return removed_edges

    def mark_tr_edges(self):
        """
        transitive reduction
        """
        n_mark = self.n_mark
        e_reduce = self.e_reduce
        FUZZ = 500
        for n in self.nodes:
            n_mark[n] = "vacant"
    
        for n_name, node in self.nodes.items():

            out_edges = node.out_edges
            if len(out_edges) == 0:
                continue
            
            out_edges.sort(key=lambda x: x.attr["length"])
            
            for e in out_edges:
                w = e.out_node
                n_mark[ w.name ] = "inplay"
            
            max_len = out_edges[-1].attr["length"]
                
            max_len += FUZZ
            
            for e in out_edges:
                e_len = e.attr["length"]
                w = e.out_node
                if n_mark[w.name] == "inplay":
                    w.out_edges.sort( key=lambda x: x.attr["length"] )
                    for e2 in w.out_edges:
                        if e2.attr["length"] + e_len < max_len:
                            x = e2.out_node
                            if n_mark[x.name] == "inplay":
                                n_mark[x.name] = "eliminated"
            
            for e in out_edges:
                e_len = e.attr["length"]
                w = e.out_node
                w.out_edges.sort( key=lambda x: x.attr["length"] )
                if len(w.out_edges) > 0:
                    x = w.out_edges[0].out_node
                    if n_mark[x.name] == "inplay":
                        n_mark[x.name] = "eliminated"
                for e2 in w.out_edges:
                    if e2.attr["length"] < FUZZ:
                        x = e2.out_node
                        if n_mark[x.name] == "inplay":
                            n_mark[x.name] = "eliminated"
                            
            for out_edge in out_edges:
                v = out_edge.in_node
                w = out_edge.out_node
                if n_mark[w.name] == "eliminated":
                    e_reduce[ (v.name, w.name) ] = True
                    v_name, w_name = reverse_end(w.name), reverse_end(v.name)
                    e_reduce[(v_name, w_name)] = True
                n_mark[w.name] = "vacant"
                

    def mark_best_overlap(self):
        """
        find the best overlapped edges
        """

        best_edges = set()
        removed_edges = set()

        for v in self.nodes:

            out_edges = self.nodes[v].out_edges
            if len(out_edges) > 0:
                out_edges.sort(key=lambda e: e.attr["score"])
                e = out_edges[-1]
                best_edges.add( (e.in_node.name, e.out_node.name) )

            in_edges = self.nodes[v].in_edges
            if len(in_edges) > 0:
                in_edges.sort(key=lambda e: e.attr["score"])
                e = in_edges[-1]
                best_edges.add( (e.in_node.name, e.out_node.name) )

        if DEBUG_LOG_LEVEL > 1:
            print "X", len(best_edges)

        for e_n, e in self.edges.items():
            v = e_n[0]
            w = e_n[1]
            if self.e_reduce[ (v, w) ] != True:
                if (v, w) not in best_edges:
                    self.e_reduce[(v, w)] = True
                    removed_edges.add( (v, w) )
                    v2, w2 = reverse_end(w), reverse_end(v)
                    self.e_reduce[(v2, w2)] = True
                    removed_edges.add( (v2, w2) )
                
        return removed_edges

    def resolve_repeat_edges(self):

        #nxsg = nx.DiGraph()
        #for v, w in self.edges:
        #    if self.e_reduce[(v, w)] != True:
        #        nxsg.add_edge(v, w)
        #nxsg_r = nxsg.reverse()

        edges_to_reduce = []
        nodes_to_test = set()
        for v_n, v in self.nodes.items():
            
            out_nodes = []
            for e in v.out_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    out_nodes.append( e.out_node.name )

            in_nodes = []
            for e in v.in_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    in_nodes.append( e.in_node.name )

            if len(out_nodes) == 1 and len(in_nodes)  == 1:
                nodes_to_test.add(v_n)
        
        for v_n in list( nodes_to_test ):
            
            v = self.nodes[v_n]

            out_nodes = []
            for e in v.out_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    out_nodes.append( e.out_node.name )

            in_nodes = []
            for e in v.in_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    in_nodes.append( e.in_node.name )

            in_node_name = in_nodes[0] 

            for out_edge in self.nodes[in_node_name].out_edges:
                vv = out_edge.in_node.name
                ww = out_edge.out_node.name

                ww_out = self.nodes[ww].out_edges
                v_out = self.nodes[v_n].out_edges
                ww_out_nodes = set( [ n.out_node.name for n in ww_out] )
                v_out_nodes = set(  [ n.out_node.name for n in v_out] )
                o_overlap = len( ww_out_nodes & v_out_nodes )

                ww_in_count = 0
                for e in self.nodes[ww].in_edges:
                    if self.e_reduce[ ( e.in_node.name, e.out_node.name ) ] == False:
                        ww_in_count += 1
                #print ww, ww_in_count 
                #G1 = nx.ego_graph( nxsg,  ww, 3, undirected=False)
                #G2 = nx.ego_graph( nxsg,  v_n, 3, undirected=False)
                #o_overlap = len( set(G1.nodes()) & set(G2.nodes()) )

                if ww != v_n and\
                   self.e_reduce[ (vv, ww) ] == False and\
                   ww_in_count > 1 and\
                   ww not in nodes_to_test and\
                   o_overlap == 0:
                    edges_to_reduce.append( (vv, ww) )

            out_node_name = out_nodes[0]

            for in_edge in self.nodes[out_node_name].in_edges:
                vv = in_edge.in_node.name
                ww = in_edge.out_node.name

                vv_in = self.nodes[vv].in_edges
                v_in = self.nodes[v_n].in_edges
                vv_in_nodes = set( [ n.in_node.name for n in vv_in] )
                v_in_nodes = set(  [ n.in_node.name for n in v_in] )
                i_overlap = len( vv_in_nodes & v_in_nodes )

                vv_out_count = 0
                for e in self.nodes[vv].out_edges:
                    if self.e_reduce[ ( e.in_node.name, e.out_node.name )] == False:
                        vv_out_count += 1
                #print vv, vv_out_count 
                #G1 = nx.ego_graph( nxsg_r,  vv, 3, undirected=False)
                #G2 = nx.ego_graph( nxsg_r,  v_n, 3, undirected=False)
                #i_overlap = len( set(G1.nodes()) & set(G2.nodes()) )

                if vv != v_n and\
                   self.e_reduce[ (vv, ww) ] == False and\
                   vv_out_count > 1 and\
                   vv not in nodes_to_test and\
                   i_overlap == 0:
                    edges_to_reduce.append( (vv, ww) )

        removed_edges = set()
        for e in edges_to_reduce:
            self.e_reduce[e] = True
            removed_edges.add(e)

        return removed_edges


    def sub_in(self, prev_nodes):   # 与えられたnodesから1つin_nodeを辿って到達可能なnodesを求める
        rtn_nodes = set()
        for n in prev_nodes:
            #print "prev node:", n
            for e in self.nodes[n].in_edges:
                #print "in edge:", e
                #print e.in_node.name
                if self.e_reduce[ (e.in_node.name, e.out_node.name) ] == False:   # reduction済みのものは考慮に入れるか
                    rtn_nodes.add(e.in_node.name)
        return rtn_nodes

    def traverse_in_nodes(self, name):
        #print "node name:", name
        rtn = set()
        current_nodes = set()
        current_nodes.add(name)

        for i in range(3):   # traverse 3 or less nodes
            current_nodes = self.sub_in(current_nodes)
            rtn |= current_nodes
            #print current_nodes

        #print "in_rtn:", rtn
        return rtn

    def sub_out(self, prev_nodes):
        rtn_nodes = set()
        for n in prev_nodes:
            #print "prev node:", n
            for e in self.nodes[n].out_edges:
                #print "in edge:", e
                #print e.in_node.name
                if self.e_reduce[ (e.in_node.name, e.out_node.name) ] == False:
                    rtn_nodes.add(e.out_node.name)
        return rtn_nodes

    def traverse_out_nodes(self, name):
        #print "node name:", name
        rtn = set()
        current_nodes = set()
        current_nodes.add(name)

        for i in range(3):   # traverse 3 or less nodes
            current_nodes = self.sub_out(current_nodes)
            rtn |= current_nodes
            #print current_nodes

        #print "out_rtn:", rtn
        return rtn
    
    def mark_putative_fp_edges(self):   # new algorithm for detecting putative FP overlaps
        edges_to_reduce = []

        for e_n, e in self.edges.items():
            v = e_n[0]
            w = e_n[1]
            if self.e_reduce[ (v, w) ] == True:
                continue

            v_in_count = 0
            for e in self.nodes[v].in_edges:
                if self.e_reduce[ (e.in_node.name, e.out_node.name) ] == False:
                    v_in_count += 1
            v_out_count = 0
            for e in self.nodes[v].out_edges:
                if self.e_reduce[ (e.in_node.name, e.out_node.name) ] == False:
                    v_out_count += 1
            w_in_count = 0
            for e in self.nodes[w].in_edges:
                if self.e_reduce[ (e.in_node.name, e.out_node.name) ] == False:
                    w_in_count += 1
            w_out_count = 0
            for e in self.nodes[w].out_edges:
                if self.e_reduce[ (e.in_node.name, e.out_node.name) ] == False:
                    w_out_count += 1
            
            if v_in_count == 0 or v_out_count == 0 or w_in_count == 0 or w_out_count == 0:
                continue

            #print "\n--- candidate:", v, w, "---"

            v_in_nodes = self.traverse_in_nodes(v)
            v_out_nodes = self.traverse_out_nodes(v)
            w_in_nodes = self.traverse_in_nodes(w)
            w_out_nodes = self.traverse_out_nodes(w)
            #print v_in_nodes
            #print v_out_nodes
            #print w_in_nodes
            #print w_out_nodes
            
            if len((v_in_nodes | v_out_nodes) & (w_in_nodes | w_out_nodes)) == 0:
                print "*****putative FP edge:", v, w, "*****"
                edges_to_reduce.append( (v, w) )

        fp_edges = set()
        for e in edges_to_reduce:
            self.e_reduce[e] = True
            fp_edges.add(e)

        print fp_edges

        return fp_edges
    

    def get_out_edges_for_node(self, name, mask=True):
        rtn = []
        for e in self.nodes[name].out_edges:
            v = e.in_node
            w = e.out_node
            if self.e_reduce[ (v.name, w.name) ] == False:
                rtn.append(e)
        return rtn
        
        
    def get_in_edges_for_node(self, name, mask=True):
        rtn = []
        for e in self.nodes[name].in_edges:
            v = e.in_node
            w = e.out_node
            if self.e_reduce[ (v.name, w.name) ] == False:
                rtn.append(e)
        return rtn

    def get_best_out_edge_for_node(self, name, mask=True):
        rtn = []
        for e in self.nodes[name].out_edges:
            v = e.in_node
            w = e.out_node
            if self.e_reduce[ (v.name, w.name) ] == False:
                rtn.append(e)
        rtn.sort(key=lambda e: e.attr["score"])

        return rtn[-1]

    def get_best_in_edge_for_node(self, name, mask=True):
        rtn = []
        for e in self.nodes[name].in_edges:
            v = e.in_node
            w = e.out_node
            if self.e_reduce[ (v.name, w.name) ] == False:
                rtn.append(e)
        rtn.sort(key=lambda e: e.attr["score"])
        return rtn[-1]
        

RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))
def generate_seq_from_path(sg, seqs, path):
    subseqs = []
    r_id, end = path[0].split(":")
    
    count = 0
    for i in range( len( path ) -1 ):
        w_n, v_n = path[i:i+2]
        edge = sg.edges[ (w_n, v_n ) ]
        read_id, coor = edge.attr["label"].split(":")
        b,e = coor.split("-")
        b = int(b)
        e = int(e)
        if b < e:
            subseqs.append( seqs[read_id][b:e] )
        else:
            subseqs.append( "".join( [RCMAP[c] for c in seqs[read_id][b:e:-1]] ) )

    return "".join(subseqs)

def reverse_end( node_id ):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end

def reverse_edge( e ):
    e1, e2 = e
    return reverse_end(e2), reverse_end(e1)

def reverse_path( p ):
    p = p[::-1]
    return [reverse_end(n) for n in p]

    
def find_bundle(ug, u_edge_data, start_node, depth_cutoff, width_cutoff, length_cutoff):

    tips = set()
    bundle_edges = set()
    bundle_nodes = set()

    local_graph = nx.ego_graph(ug, start_node, depth_cutoff, undirected=False)
    length_to_node = {start_node:0}
    score_to_node = {start_node:0}

    v = start_node
    end_node = start_node

    if DEBUG_LOG_LEVEL > 1: 
        print
        print 
        print "start", start_node

    bundle_nodes.add(v)
    for vv, ww, kk in local_graph.out_edges(v, keys = True):
        max_score = 0
        max_length = 0

        if (vv, ww, kk) not in bundle_edges and\
                reverse_end(ww) not in bundle_nodes:

            bundle_edges.add( (vv, ww, kk) )
            tips.add(ww)

    for v in list(tips):
        bundle_nodes.add(v)

    depth = 1
    width = 1.0
    converage = False


    while 1:
        if DEBUG_LOG_LEVEL > 1:
            print "# of tips", len(tips)

        if len(tips) > 4:
            converage = False
            break

        if len(tips) == 1:
            end_node = tips.pop()

            if DEBUG_LOG_LEVEL > 1:
                print "end", end_node

            if end_node not in length_to_node:
                v = end_node
                max_score_edge = None
                max_score = 0
                for uu, vv, kk in local_graph.in_edges(v, keys=True):
                    if uu not in length_to_node:
                        continue

                    score = u_edge_data[ (uu, vv, kk) ][1]

                    if score > max_score:

                        max_score = score
                        max_score_edge = (uu, vv, kk)

                length_to_node[v] = length_to_node[max_score_edge[0]] +  u_edge_data[ max_score_edge ][0]
                score_to_node[v] = score_to_node[max_score_edge[0]] +  u_edge_data[ max_score_edge ][1]
                

            converage = True
            break
        

        depth += 1
        width = 1.0 * len(bundle_edges) / depth

        if depth > 10 and width > width_cutoff:
            converage = False
            break

        if depth > depth_cutoff:
            converage = False
            break
        
        tips_list = list(tips)

        tip_updated = False
        loop_detect = False
        length_limit_reached = False

        for v in tips_list:
            if DEBUG_LOG_LEVEL > 1:
                print "process", v

            if len(local_graph.out_edges(v, keys=True)) == 0: # dead end route
                #print "no out edge", v
                continue

            max_score_edge = None
            max_score = 0

            extend_tip = True

            for uu, vv, kk in local_graph.in_edges(v, keys=True):
                if DEBUG_LOG_LEVEL > 1: 
                    print "in_edges", uu, vv, kk
                    print uu, "in length_to_node",  uu in length_to_node

                if uu not in length_to_node:
                    extend_tip = False
                    break

                score = u_edge_data[ (uu, vv, kk) ][1]

                if score > max_score:

                    max_score = score
                    max_score_edge = (uu, vv, kk)
            
            if extend_tip:
            
                length_to_node[v] = length_to_node[max_score_edge[0]] +  u_edge_data[ max_score_edge ][0]
                score_to_node[v] = score_to_node[max_score_edge[0]] +  u_edge_data[ max_score_edge ][1]

                if length_to_node[v] > length_cutoff:
                    length_limit_reached = True
                    converage = False
                    break

                v_updated = False
                for vv, ww, kk in local_graph.out_edges(v, keys=True):

                    if DEBUG_LOG_LEVEL > 1:
                        print "test", vv, ww, kk

                    if ww in length_to_node:
                        loop_detect = True
                        if DEBUG_LOG_LEVEL > 1:
                            print "loop_detect", ww
                        break

                    if (vv, ww, kk) not in bundle_edges and\
                            reverse_end(ww) not in bundle_nodes:

                        if DEBUG_LOG_LEVEL > 1:
                            print "add", ww

                        tips.add(ww)
                        bundle_edges.add( (vv, ww, kk) )
                        tip_updated = True
                        v_updated = True

                if v_updated:

                    if DEBUG_LOG_LEVEL > 1:
                        print "remove", v

                    tips.remove(v)

                    if len(tips) == 1:
                        break

            if loop_detect:
                converage = False
                break

        if length_limit_reached:
            converage = False
            break

        if loop_detect:
            converage = False
            break

        if not tip_updated:
            converage = False
            break

        for v in list(tips):
            bundle_nodes.add(v)

        

    data = start_node, end_node, bundle_edges, length_to_node[end_node], score_to_node[end_node], depth
    
    """
    start_node_r, end_node_r = reverse_end(end_node), reverse_end(start_node)
    
    bundle_edge_r = set()
    for v, w, k in list(bundle_edges):
        vv, ww, kk = reverse_end(w), reverse_end(v), reverse_end(k)
        bundle_edge_r.add( (vv, ww, kk) )

    data_r = start_node_r, end_node_r, bundle_edge_r, length_to_node[end_node], score_to_node[end_node], depth
    """
    data_r = None

    if DEBUG_LOG_LEVEL > 1:
        print converage, data, data_r
    return converage, data, data_r


def SGToNXG(sg):   # from falcon_asm.py, attr is reduced
    G=nx.DiGraph()

    #max_score = max([ sg.edges[ e ].attr["score"] for e in sg.edges if sg.e_reduce[e] != True ])   # 全てのエッジのうち最大のスコアを求める
    #out_f = open("edges_list","w")   # ついでにエッジの情報も出力
    for v, w in sg.edges:   # in_node_nameとout_node_name
        if sg.e_reduce[(v, w)] != True:   # 除去されていないエッジについてのみ出力
        ##if 1:
            out_degree = len(sg.nodes[v].out_edges)
            G.add_node( v, size = out_degree )   # sizeはノードを表示するときの大きさ
            G.add_node( w, size = out_degree )   # out_nodeの方のsizeはout_nodeから出るエッジが存在しないときのためにとりあえず設定しているだけ？
            #label = sg.edges[ (v, w) ].attr["label"]
            #score = sg.edges[ (v, w) ].attr["score"]
            #print >>out_f, v, w, label, score 
            G.add_edge( v, w )   # weightとn_weightは何のために付加している？
            #print in_node_name, out_node_name
    #out_f.close()
    return G


def generate_string_graph(args):

    overlap_file = args.overlap_file

    contained_reads = set()
    chimer_ids = set()
 
    filter_reads = False
    
    seqs = set()

    G=nx.Graph()
    edges =set()
    overlap_data = []
    contained_reads = set()
    overlap_count = {}


    # loop through the overlapping data to load the data in the a python array
    # contained reads are identified 

    with open(overlap_file) as f:
        for l in f:
            l = l.strip().split()

            #work around for some ill formed data recored
            #if len(l) != 13:
            #    continue
            
            f_id, g_id, score, identity = l[:4]

            if f_id == g_id:  # don't need self-self overlapping
                continue
            
            if filter_reads:

                if g_id not in seqs: 
                    continue

                if f_id not in seqs:
                    continue

            score = int(score)
            identity = float(identity)
            contained = l[12]
            if contained == "contained":
                contained_reads.add(f_id)
                continue
            if contained == "contains":
                contained_reads.add(g_id)
                continue
            if contained == "none":
                continue

            if identity < args.min_idt: # only take record with >96% identity as overlapped reads
                continue
            #if score > -2000:
            #    continue
            f_strain, f_start, f_end, f_len = (int(c) for c in l[4:8])
            g_strain, g_start, g_end, g_len = (int(c) for c in l[8:12])

            # only used reads longer than the 4kb for assembly
            if f_len < args.min_len: continue
            if g_len < args.min_len: continue
            
            # double check for proper overlap
            """
            if f_start > 24 and f_len - f_end > 24:  # allow 24 base tolerance on both sides of the overlapping
                continue
            
            if g_start > 24 and g_len - g_end > 24:
                continue
            
            if g_strain == 0:
                if f_start < 24 and g_len - g_end > 24:
                    continue
                if g_start < 24 and f_len - f_end > 24:
                    continue
            else:
                if f_start < 24 and g_start > 24:
                    continue
                if g_start < 24 and f_start > 24:
                    continue
            """

            overlap_data.append( (f_id, g_id, score, identity,
                                  f_strain, f_start, f_end, f_len,
                                  g_strain, g_start, g_end, g_len) )

            overlap_count[f_id] = overlap_count.get(f_id,0)+1
            overlap_count[g_id] = overlap_count.get(g_id,0)+1
            
    #print "###", len(overlap_data), len(contained_reads)
    overlap_set = set()
    sg = StringGraph()
    for od in overlap_data:
        f_id, g_id, score, identity = od[:4]
        if f_id in contained_reads:
            continue
        if g_id in contained_reads:
            continue
        f_s, f_b, f_e, f_l = od[4:8]
        g_s, g_b, g_e, g_l = od[8:12]
        overlap_pair = [f_id, g_id]
        overlap_pair.sort()
        overlap_pair = tuple( overlap_pair )
        if overlap_pair in overlap_set:  # don't allow duplicated records
            continue
        else:
            overlap_set.add(overlap_pair)

        
        if g_s == 1: # revered alignment, swapping the begin and end coordinates
            g_b, g_e = g_e, g_b
        
        # build the string graph edges for each overlap
        if f_b > 1:
            if g_b < g_e:
                """
                     f.B         f.E
                  f  ----------->
                  g         ------------->
                            g.B           g.E
                """
                if f_b == 0 or g_e - g_l == 0:
                    continue
                sg.add_edge( "%s:B" % g_id, "%s:B" % f_id, label = (f_id, f_b, 0), 
                                                           length = abs(f_b-0),
                                                           score = -score, 
                                                           identity = identity )
                sg.add_edge( "%s:E" % f_id, "%s:E" % g_id, label = (g_id, g_e, g_l), 
                                                           length = abs(g_e-g_l),
                                                           score = -score,
                                                           identity = identity)
            else:
                """
                     f.B         f.E
                  f  ----------->
                  g         <-------------
                            g.E           g.B           
                """
                if f_b == 0 or g_e == 0:
                    continue
                sg.add_edge( "%s:E" % g_id, "%s:B" % f_id, label = (f_id, f_b, 0), 
                                                           length = abs(f_b -0),
                                                           score = -score,
                                                           identity = identity)
                sg.add_edge( "%s:E" % f_id, "%s:B" % g_id, label = (g_id, g_e, 0), 
                                                           length = abs(g_e- 0),
                                                           score = -score,
                                                           identity = identity)
        else:
            if g_b < g_e:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         ------------->
                            g.B           g.E
                """
                if g_b == 0 or f_e - f_l == 0:
                    continue
                sg.add_edge( "%s:B" % f_id, "%s:B" % g_id, label = (g_id, g_b, 0), 
                                                           length = abs(g_b - 0),
                                                           score = -score,
                                                           identity = identity)
                sg.add_edge( "%s:E" % g_id, "%s:E" % f_id, label = (f_id, f_e, f_l), 
                                                           length = abs(f_e-f_l),
                                                           score = -score,
                                                           identity = identity)
            else:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         <-------------
                            g.E           g.B           
                """
                if g_b - g_l == 0 or f_e - f_l ==0:
                    continue
                sg.add_edge( "%s:B" % f_id, "%s:E" % g_id, label = (g_id, g_b, g_l), 
                                                           length = abs(g_b - g_l),
                                                           score = -score,
                                                           identity = identity)
                sg.add_edge( "%s:B" % g_id, "%s:E" % f_id, label = (f_id, f_e, f_l), 
                                                           length = abs(f_e - f_l),
                                                           score = -score,
                                                           identity = identity)


    sg.init_reduce_dict()
    #G = SGToNXG(sg)
    #nx.write_gexf(G, "overlap_graph.gexf")

    #if not args.disable_chimer_prediction:
    #    sg.mark_chimer_edges()
    #sg.mark_spur_edge()
    

    sg.mark_tr_edges() # mark those edges that transitive redundant
    #G = SGToNXG(sg)
    #nx.write_gexf(G, "transitive_reduced.gexf")

    if DEBUG_LOG_LEVEL > 1:
        print sum( [1 for c in sg.e_reduce.values() if c == True] )
        print sum( [1 for c in sg.e_reduce.values() if c == False] )


    removed_edges = set()
    if args.lfc == True:
        removed_edges = sg.resolve_repeat_edges()  
    else:
        removed_edges = sg.mark_best_overlap() # mark those edges that are best overlap edges

    #G = SGToNXG(sg)
    #nx.write_gexf(G, "heuristic_reduced.gexf")

    ## new FP overlap prediction
    #fp_edges = set()
    #if args.predict_fp_edges == True:
    #    fp_edges = sg.mark_putative_fp_edges()
    #    G = SGToNXG(sg)
    #    nx.write_gexf(G, "FP_reduced.gexf")
    
    spur_edges = sg.mark_spur_edge()
    #G = SGToNXG(sg)
    #nx.write_gexf(G, "spur_reduced.gexf")

    if DEBUG_LOG_LEVEL > 1:
        print sum( [1 for c in sg.e_reduce.values() if c == False] )

    #max_score = max([ sg.edges[ e ].attr["score"] for e in sg.edges if sg.e_reduce[e] != True ])
    out_f = open("sg_edges_list", "w")
    nxsg = nx.DiGraph()
    edge_data = {}
    for v, w in sg.edges:
        e = sg.edges[ (v, w) ]
        rid, sp, tp = e.attr["label"]
        score = e.attr["score"]
        identity = e.attr["identity"]
        length = abs(sp-tp)

        if  sg.e_reduce[(v, w)] != True:
            type_ = "G"
            label = "%s:%d-%d" % (rid, sp, tp)
            nxsg.add_edge(v, w, label = label, length = length, score = score)
            edge_data[ (v, w) ] = (rid, sp, tp, length, score, identity, type_)
        elif (v, w) in removed_edges:
            type_ = "R"
        elif (v, w) in spur_edges:
            type_ = "S"
        #elif (v, w) in fp_edges:   # new condition
        #    type_ = "FP"
        elif sg.e_reduce[(v, w)] == True:
            type_ = "TR"

        print >>out_f, v, w, rid, sp, tp, score, identity, type_
        
    out_f.close()
    nxsg_r = nxsg.reverse()    

    return nxsg, nxsg_r, edge_data



def construct_compound_paths(ug, u_edge_data):

    source_nodes = set()
    sink_nodes = set()
    simple_nodes = set()
    branch_nodes = set()

    all_nodes = ug.nodes()
    for n in all_nodes:
        in_degree = len( ug.in_edges(n) )
        out_degree = len( ug.out_edges(n) )
        if in_degree == 0:
            source_nodes.add(n)
        if out_degree == 0:
            sink_nodes.add(n)
        if in_degree == 1 and out_degree == 1:
            simple_nodes.add(n)
        if in_degree > 1 or out_degree > 1:
            branch_nodes.add(n)

    #print "#", len(all_nodes),len(source_nodes), len(sink_nodes), len(simple_nodes), len(branch_nodes)
    compound_paths_0 = []
    for p in list(branch_nodes):
        if ug.out_degree(p) > 1:
            coverage, data, data_r =  find_bundle(ug, u_edge_data, p, 32, 8, 500000)
            if coverage == True:
                start_node, end_node, bundle_edges, length, score, depth = data
                compound_paths_0.append(  (start_node, "NA", end_node, 1.0*len(bundle_edges)/depth, length, score, bundle_edges ) )

    compound_paths_0.sort( key=lambda x: -x[5] )


    edge_to_cpath = {}
    compound_paths_1 = {}
    for s, v, t, width, length, score, bundle_edges in compound_paths_0:
        overlapped = False
        for vv, ww, kk in list(bundle_edges):
            if (vv, ww, kk) in edge_to_cpath:
                overlapped = True
                break
            rvv = reverse_end(vv)
            rww = reverse_end(ww)
            rkk = reverse_end(kk)
            if (vv, ww, kk) in edge_to_cpath:
                overlapped = True
                break
            

        if not overlapped:
            if DEBUG_LOG_LEVEL > 1:
                print "constructing", s,v, t

            bundle_edges_r = []
            rs = reverse_end(t)
            rt = reverse_end(s)

            for vv, ww, kk in list(bundle_edges):
                edge_to_cpath.setdefault( (vv, ww, kk), set() )
                edge_to_cpath[ (vv, ww, kk) ].add( ( s, t, v) )
                rvv = reverse_end(ww)
                rww = reverse_end(vv)
                rkk = reverse_end(kk)
                edge_to_cpath.setdefault( (rvv, rww, rkk), set() )
                edge_to_cpath[ (rvv, rww, rkk) ].add( (rs, rt, v) ) #assert v == "NA"
                bundle_edges_r.append(  (rvv, rww, rkk) )
            
            compound_paths_1[ ( s, v, t) ] = width, length, score, bundle_edges
            compound_paths_1[ ( rs, v, rt) ] = width, length, score, bundle_edges_r

             
    compound_paths_2 = {}
    edge_to_cpath = {}
    for s, v, t in compound_paths_1:
        rs = reverse_end(t)
        rt = reverse_end(s)
        if (rs, "NA", rt) not in compound_paths_1:
            if DEBUG_LOG_LEVEL > 1:
                print "non_compliment bundle", s, v, t, len(compound_paths_1[( s, v, t)][-1])
            continue
        width, length, score, bundle_edges = compound_paths_1[ (s, v, t) ]
        compound_paths_2[ (s, v, t) ] = width, length, score, bundle_edges
        for vv, ww, kk in list(bundle_edges):
            edge_to_cpath.setdefault( (vv, ww, kk), set() )
            edge_to_cpath[ (vv, ww, kk) ].add( ( s, t, v) )


    compound_paths_3 = {}
    for k, val in compound_paths_2.items():
        
        start_node, NA, end_node = k
        rs = reverse_end(end_node)
        rt = reverse_end(start_node)
        assert (rs, "NA", rt) in compound_paths_2
        
        contained = False
        for vv, ww, kk in ug.out_edges(start_node, keys=True):
            if len(edge_to_cpath.get( (vv, ww, kk), [] )) > 1: 
                contained = True

        if not contained:
            compound_paths_3[k] = val
            if DEBUG_LOG_LEVEL > 1:
                print "compound", k 

    compound_paths = {}
    for s, v, t in compound_paths_3:
        rs = reverse_end(t)
        rt = reverse_end(s)
        if (rs, "NA", rt) not in compound_paths_3:
            continue
        compound_paths[ (s, v, t) ] = compound_paths_3[ (s, v, t) ]

    return compound_paths


def jaccard(s, t):

    if len(s) == 0 or len(t) == 0:
        return 0.0
    
    intersection_len = len(s.intersection(t))
    union_len = len(s) + len(t) - intersection_len
    return float(intersection_len) / union_len


if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser(description='a example string graph assembler that is desinged for handling diploid genomes')
    parser.add_argument('overlap_file', help='a file that contains the overlap information.')
    parser.add_argument('unitig_path_file', help='a file that contains untiig path file. that is, ctg_paths file of assembly with simple_nodes')

    parser.add_argument('--min_len', type=int, default=2200, 
                        help='minimum length of the reads to be considered for assembling')
    parser.add_argument('--min_idt', type=float, default=96,
                        help='minimum alignment identity of the reads to be considered for assembling')
    parser.add_argument('--disable_chimer_prediction', action="store_true", default=True,
                        help='you may want to disable this as some reads can be falsely identified as chimers in low coverage case')
    parser.add_argument('--lfc', action="store_true", default=False,
                        help='use local flow constraint method rather than best overlap method to resolve knots in string graph')
    parser.add_argument('--predict_fp_edges', action="store_true", default=False,
                        help='use FP overlaps prediction')

    args = parser.parse_args()
    
    unitig_path_file = args.unitig_path_file

    # transitivity reduction, remove spurs, remove putative edges caused by repeats
    sg, sg_r, edge_data = generate_string_graph(args)


    simple_paths = {}
    dual_path = {}


    sg2 = nx.DiGraph()

    for v, w in edge_data:

        assert (reverse_end(w), reverse_end(v)) in edge_data
        
        #if (v, w) in masked_edges:
        #    continue

        rid, sp, tp, length, score, identity, type_ = edge_data[ (v, w) ]
        if type_ != "G":
            continue

        label = "%s:%d-%d" % (rid, sp, tp)
        sg2.add_edge( v, w, label = label, length = length, score = score)

        
    # utg construction phase 1, identify all simple paths
    s_nodes = set()
    t_nodes = set()
    simple_nodes = set()

    all_nodes = sg2.nodes()
    for n in all_nodes:
        in_degree = len( sg2.in_edges(n) )
        out_degree = len( sg2.out_edges(n) )
        if in_degree == 1 and out_degree == 1:
            simple_nodes.add(n)
        else:
            if out_degree != 0:
                s_nodes.add(n)
            if in_degree != 0:
                t_nodes.add(n)


    free_edges = set(sg2.edges())

    if DEBUG_LOG_LEVEL > 1: 
        for s in list(simple_nodes):
            print "simple_node", s
        for s in list(s_nodes):
            print "s_node", s
        for s in list(t_nodes):
            print "t_node", s

        for v,w in free_edges:
            if (reverse_end(w), reverse_end(v) ) not in free_edges:
                print "bug", v,w
                print oreverse_end(w), reverse_end(v)

    while len(free_edges) != 0:
        if len(s_nodes) != 0:
            n = s_nodes.pop()
            if DEBUG_LOG_LEVEL > 1:
                print "initial utg 1", n
        else:
            e = free_edges.pop()
            free_edges.add(e)
            n = e[0]
            if DEBUG_LOG_LEVEL > 1:
                print "initial utg 2", n

        path = []
        path_length =0
        path_score = 0 
        for v, w in sg2.out_edges(n):
            #print "test",v,w,path
            if (v, w) not in free_edges:
                continue
            rv = reverse_end(v)
            rw = reverse_end(w)
            #if (rw, rv) not in free_edges:
            #    continue

            path_length = 0
            path_score = 0
            v0 = v
            w0 = w
            path = [v, w]
            path_edges = set()
            path_edges.add( (v, w) )
            path_length += edge_data[ (v, w) ][3]
            path_score += edge_data[ (v, w) ][4]
            free_edges.remove( (v, w) )

            r_path_length = 0
            r_path_score = 0
            rv0 = rv
            rw0 = rw
            r_path = [rv, rw] # need to reverse again
            r_path_edges = set()
            r_path_edges.add( (rw, rv) )
            r_path_length += edge_data[ (rw, rv) ][3]
            r_path_score += edge_data[ (rw, rv) ][4]
            free_edges.remove( (rw, rv) )

            while w in simple_nodes:
                w, w_ = sg2.out_edges(w)[0]
                if (w, w_) not in free_edges:
                    break
                rw_, rw = reverse_end(w_), reverse_end(w)

                if ( rw_, rw ) in path_edges:
                    break

                path.append(w_)
                path_edges.add( (w, w_) )
                path_length += edge_data[ (w, w_) ][3]
                path_score += edge_data[ (w, w_) ][4]
                free_edges.remove( (w, w_) )
                
                r_path.append(rw_)
                r_path_edges.add( (rw_, rw) )
                r_path_length += edge_data[ (rw_, rw) ][3]
                r_path_score += edge_data[ (rw_, rw) ][4]
                free_edges.remove( (rw_, rw) )
                

                w = w_

            simple_paths[ (v0, w0, path[-1]) ] = path_length, path_score, path
            r_path.reverse()
            assert r_path[0] == reverse_end(path[-1])
            simple_paths[ (r_path[0], rw0, rv0) ] = r_path_length, r_path_score, r_path

            if DEBUG_LOG_LEVEL > 1:
                print  path_length, path_score, path

            dual_path[ (r_path[0], rw0, rv0) ] = (v0, w0, path[-1])
            dual_path[ (v0, w0, path[-1]) ] = (r_path[0], rw0, rv0)



    ug = nx.MultiDiGraph()
    u_edge_data = {}
    circular_path = set()

    for s, v, t in simple_paths:
        length, score, path = simple_paths[ (s, v, t) ]
        u_edge_data[ (s, t, v) ] = (length, score, path, "simple")
        if s != t:
            ug.add_edge(s, t, key = v, type_ = "simple", via = v, length = length, score = score)
        else:
            circular_path.add( (s, t, v) )

    #nx.write_gexf(ug, "unitig_graph_1.gexf")

    if DEBUG_LOG_LEVEL > 1:
        with open("utg_data0","w") as f:
            for s, t, v in u_edge_data:
                rs = reverse_end(t)
                rt = reverse_end(s)
                rv = reverse_end(v)
                assert (rs, rt, rv) in u_edge_data
                length, score, path_or_edges, type_ = u_edge_data[ (s, t, v) ]
                
                if type_ == "compound":
                    path_or_edges = "|".join( [ ss+"~"+vv+"~"+tt for ss, tt, vv in path_or_edges ] )
                else:
                    path_or_edges = "~".join( path_or_edges )
                print >>f, s, v, t, type_, length, score, path_or_edges

    utg_spurs = set()
    all_nodes = ug.nodes()
    excessive_link_nodes = set()

    for n in all_nodes:
        in_degree = len( set( e[0] for e in ug.in_edges(n))  ) # ignore mutli-edges
        out_degree = len( set( e[1] for e in ug.out_edges(n)) )
        if in_degree == 1 and out_degree == 1:
            simple_nodes.add(n)
        else:
            if out_degree != 0:
                s_nodes.add(n)
            if in_degree != 0:
                t_nodes.add(n)

            if out_degree >= 2 and in_degree == 0:
                excessive_link_nodes.add(n)
            if in_degree >= 2 and out_degree == 0:
                excessive_link_nodes.add(n)
            if out_degree > 3 and in_degree < 2:
                excessive_link_nodes.add(n)
            if in_degree > 3 and out_degree < 2:
                excessive_link_nodes.add(n)
            if out_degree >= 3 and in_degree >= 3:
                excessive_link_nodes.add(n)
            if in_degree >= 3 and out_degree >= 3:
                excessive_link_nodes.add(n)

    spur_edges = set()

    ug2 = ug.copy()
    edges_to_remove = set()

    for n in list(excessive_link_nodes):
        in_degree = len( set( e[0] for e in ug.in_edges(n))  ) # ignore mutli-edges
        out_degree = len( set( e[1] for e in ug.out_edges(n)) )

        for s, t, v in ug.in_edges(n, keys=True):
            length, score, edges, type_ = u_edge_data[ (s, t, v) ]
            if length > 30000:
                continue
            u_edge_data[ (s, t, v) ] = length, score, edges, "spur:1"
            if DEBUG_LOG_LEVEL > 1:
                print "spur:1", s,t,v, in_degree, out_degree
            spur_edges.add( (s, t, v) )
            edges_to_remove.add( (s, t, v) )
            rs = reverse_end(t)
            rt = reverse_end(s)
            rv = reverse_end(v)
            edges_to_remove.add( (rs, rt, rv) )
            length, score, edges, type_ = u_edge_data[ (rs, rt, rv) ]
            u_edge_data[ (rs, rt, rv) ] = length, score, edges, "spur:1"

        for s, t, v in ug.out_edges(n, keys=True):
            length, score, edges, type_ = u_edge_data[ (s, t, v) ]
            if length > 30000:
                continue
            u_edge_data[ (s, t, v) ] = length, score, edges, "spur:1"
            if DEBUG_LOG_LEVEL > 1:
                print "spur:1", s,t,v, in_degree, out_degree
            spur_edges.add( (s, t, v) )
            edges_to_remove.add( (s, t, v) )
            rs = reverse_end(t)
            rt = reverse_end(s)
            rv = reverse_end(v)
            edges_to_remove.add( (rs, rt, rv) )
            length, score, edges, type_ = u_edge_data[ (rs, rt, rv) ]
            u_edge_data[ (rs, rt, rv) ] = length, score, edges, "spur:1"

    for n in s_nodes:
        if ug.in_degree(n) != 0:
            continue
        for s, t, v in ug.out_edges(n, keys=True):
            length, score, edges, type_ = u_edge_data[ (s, t, v) ]
            if length > 30000:
                continue
            in_degree = len( set( e[0] for e in ug.in_edges(t))  ) # ignore mutli-edges
            out_degree = len( set( e[1] for e in ug.out_edges(t)) )
            if in_degree > 1 and out_degree > 0:
                spur_edges.add( (s, t, v) )
                edges_to_remove.add( (s, t, v) )
                u_edge_data[ (s, t, v) ] = length, score, edges, "spur:2"
                rs = reverse_end(t)
                rt = reverse_end(s)
                rv = reverse_end(v)
                edges_to_remove.add( (rs, rt, rv) )
                length, score, edges, type_ = u_edge_data[ (rs, rt, rv) ]
                u_edge_data[ (rs, rt, rv) ] = length, score, edges, "spur:2"

    for n in t_nodes:
        if ug.out_degree(n) != 0:
            continue
        for s, t, v in ug.in_edges(n, keys=True):
            length, score, edges, type_ = u_edge_data[ (s, t, v) ]
            if length > 30000:
                continue
            in_degree = len( set( e[0] for e in ug.in_edges(s))  ) # ignore mutli-edges
            out_degree = len( set( e[1] for e in ug.out_edges(s)) )
            if in_degree > 0 and out_degree > 1:
                spur_edges.add( (s, t, v) )
                edges_to_remove.add( (s, t, v) )
                u_edge_data[ (s, t, v) ] = length, score, edges, "spur:2"
                rs = reverse_end(t)
                rt = reverse_end(s)
                rv = reverse_end(v)
                edges_to_remove.add( (rs, rt, rv) )
                length, score, edges, type_ = u_edge_data[ (rs, rt, rv) ]
                u_edge_data[ (rs, rt, rv) ] = length, score, edges, "spur:2"

    for s, t, v in list(edges_to_remove):
        ug2.remove_edge( s, t, key= v)

    #nx.write_gexf(ug2, "unitig_graph_2.gexf")

    #phase 2, finding all "consistent" compound paths
    compound_paths = construct_compound_paths(ug2, u_edge_data)
    compound_path_file = open("c_path","w")

    ug2_edges = set(ug2.edges(keys = True))
    edges_to_remove  = set()
    for s, v, t in compound_paths:
        width, length, score, bundle_edges =  compound_paths[ (s, v, t) ] 
        print >> compound_path_file, s,v,t, width, length, score, "|".join( [e[0]+"~"+e[2]+"~"+e[1] for e in bundle_edges] )
        for ss, tt, vv in bundle_edges:
            if (ss, tt, vv) in ug2_edges:
                edges_to_remove.add( (ss, tt, vv) )

    
    for s, t, v in edges_to_remove:
        ug2.remove_edge( s, t ,v )
        length, score, edges, type_ = u_edge_data[ (s, t, v) ]
        if type_ != "spur":
            u_edge_data[ (s, t, v) ] = length, score, edges, "contained"


    for s, v, t in compound_paths:
        width, length, score, bundle_edges =  compound_paths[ (s, v, t) ] 
        u_edge_data[ (s, t, v) ] = (length, score, bundle_edges, "compound")
        ug2.add_edge( s, t, key = v, via = v, type_="compound", length = length, score = score)

        assert v == "NA"
        rs = reverse_end(t)
        rt = reverse_end(s)
        assert (rs, v, rt) in compound_paths
        dual_path[ (s, v, t) ] = (rs, v, rt)
        dual_path[ (rs, v, rt) ] = (s, v, t)

    compound_path_file.close()

    ug = ug2
    #nx.write_gexf(ug, "unitig_graph_3.gexf")   # spur reduction, compound reductionが終わった後のunitig graph

    # identify spurs in the utg graph
    # Currently, we use ad-hoc logic filtering out shorter utg, but we ca
    # add proper alignment comparison later to remove redundant utgs 


    with open("utg_data","w") as f:
        for s, t, v in u_edge_data:
            length, score, path_or_edges, type_ = u_edge_data[ (s, t, v) ]
            
            if v == "NA":
                path_or_edges = "|".join( [ ss+"~"+vv+"~"+tt for ss, tt, vv in path_or_edges ] )
            else:
                path_or_edges = "~".join( path_or_edges )
            print >>f, s, v, t, type_, length, score, path_or_edges


    # contig construction from utgs

    s_nodes = set()
    t_nodes = set()
    simple_nodes = set()
    simple_out = set()
    simple_in = set()

    all_nodes = ug.nodes()
    for n in all_nodes:
        in_degree = len( ug.in_edges(n) )
        out_degree = len( ug.out_edges(n) )
        if in_degree == 1 and out_degree == 1:
            simple_nodes.add(n)
        else:
            if out_degree != 0:
                s_nodes.add(n)
            if in_degree != 0:
                t_nodes.add(n)
        if out_degree == 1:
            simple_out.add(n)
        if in_degree == 1:
            simple_in.add(n)

    all_nodes = set(all_nodes)
    c_path = []
    
    free_edges = set()
    for s, t, v in ug.edges(keys=True):
        free_edges.add( (s, t, v) )

    while len(free_edges) != 0:

        if len(s_nodes) != 0:
            n = s_nodes.pop()
            #print "from s_nodes"
        else:
            e = free_edges.pop()
            n = e[0]
            print "from free_edges"
            print ug.in_edges(n)
            traced_nodes = set()
            while len(ug.in_edges(n)) == 1:   # c_pathがすべてsimple_pathになるように細工
                print "CAUTION: starting to trace in simple path!", n
                print ug.in_edges(n)[0]
                if ug.in_edges(n)[0][0] in traced_nodes:
                    break
                else:
                    traced_nodes.add(ug.in_edges(n)[0][0])
                n = ug.in_edges(n)[0][0]
            

            #print "no more s_nodes. extract randomly from free_edges.", n

        
     
        for s, t, v in ug.out_edges(n, keys=True):

            path_start = n
            path_end = None
            path_key = None
            path = []
            path_length = 0
            path_score = 0
            path_nodes = set()
            path_nodes.add(s)
            if DEBUG_LOG_LEVEL > 1:
                print "check 1", s, t, v
            path_key = t
            t0 = s
            while t in simple_nodes:   # simple_out or simple_nodes?
                if t in path_nodes:
                    break
                rt = reverse_end(t)
                if rt in path_nodes:
                    break

                path.append( (t0, t, v) )
                path_nodes.add(t)
                length, score, path_or_edges, type_ = u_edge_data[ (t0, t, v) ]
                path_length += length
                path_score += score
                assert len( ug.out_edges( t, keys=True ) ) == 1 # t is "simple_out" node
                t0, t, v = ug.out_edges( t, keys=True )[0] 

            path.append( (t0, t, v) )
            length, score, path_or_edges, type_ = u_edge_data[ (t0, t, v) ]
            path_length += length
            path_score += score
            path_nodes.add(t)
            path_end = t

            c_path.append( (path_start, path_key, path_end, path_length, path_score, path, len(path)) ) 
            #if DEBUG_LOG_LEVEL > 1:
            #print "c_path", path_start, path_key, path_end, path_length, path_score, len(path)
            for e in path:
                if e in free_edges:
                    free_edges.remove( e )
                #else:
                #    print e, "is not in free_edges"

        #print c_path[len(c_path) - 1]
 
    if DEBUG_LOG_LEVEL > 1:
        print "left over edges:", len(free_edges)


    ug4 = nx.DiGraph()   # unitig_graph_3をsimple path reductionしたdirected graph
    for path_start, path_key, path_end, p_len, p_score, path, n_edges in c_path:
        ug4.add_edge(path_start, path_end)
    nx.write_gexf(ug4, "unitig_graph_4.gexf")




    ## bin情報をmappingしてunitig graphをuntangleする

    #################################################################################################################################################################################どちらか1つをコメントイン

    ## MetaBATのbin情報を取得(fofn)
    
    metabat_bin = {}   # {u_name : bin_id}
    with open("metabat_bins.fofn", 'r') as f:
        for bin_file in f:
            bin_id = re.search(r'\.bin\.(.*?)\.', bin_file).groups()[0]

            u_names = set()
            with open(bin_file.strip(), 'r') as g:
                for line in g:
                    if re.search(r'^>', line):
                        u_names.add(re.search(r'^>(.*?)\|', line).groups()[0])
                
            for u_name in u_names:
                if "F" not in u_name:
                    #print "circular contig"
                    continue
                elif "-" in u_name:
                    #print "associated contig"
                    continue
                metabat_bin[u_name] = bin_id
                metabat_bin[u_name.replace("F", "R")] = bin_id   # rc pathもbin情報を貼り付け、しかしsimple_nodesなので意味はない？
                # Rのパスが出力contigに含まれる可能性も一応考える
    
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    """
    ## MetaBATのbin情報を取得(csv)
    metabat_bin = {}
    with open("metabat_bin_assignment.csv", 'r') as f:
        f.readline()
        for line in f:
            u_name_all, bin_id = line.strip().split(',')
            u_name = re.search(r'(.*?)\|', u_name_all).groups()[0]
            if "F" not in u_name or "-" in u_name:   # circular contig or associated contig
                continue
            metabat_bin[u_name] = bin_id
    """
    ##################################################################################################################################################################################
    """
    ## モチーフ情報を取得、使う予定無し、使うとしてもbinningの枠組みに取り込む
    motifs = []
    with open("summary_motifs.mc", 'r') as f:
        for line in f:
            if re.search(r'^[0-9]', line):   # contig名
                contig_name = line.strip()
                motifs.append( (contig_name, set()) )
            else:   # motif
                if re.search(r'(.*?)\t.*', line):
                    motif = re.search(r'(.*?)\t.*', line).groups()[0]
                elif re.search(r'(.*?) .*', line):
                    motif = re.search(r'(.*?) .*', line).groups()[0]
                else:
                    motif = re.search(r'(.*?)\n', line).groups()[0]
                if motif != "":
                    motifs[-1][1].add(motif)
    motif_info = {}   # {u_name : motif_set}
    for data in motifs:
        u_name, motif_set = data
        motif_info[u_name] = motif_set
    """
    ## unitigの情報を取得
    unitig_label_metabat = {}   # {(path_start, path_key, path_end) : bin_id}
    #unitig_label_motifs = {}   # {(path_start, path_key, path_end) : motif_set}
    bin_id_tmp = -1
    with open(unitig_path_file, 'r') as f:   # simple_nodesによって生成されたctg_paths
        for line in f:
            u_name, type_, path_start, path_key, path_end = re.search(r'(.*?) (.*?) (.*?)~.*?~(.*?) (.*?) ', line).groups()

            if type_ != "ctg_linear":   # ctg_pathsにassociated contigの情報はない
                print "[CAUTION] non-linear unitig is in the ctg_paths file!!"
                continue

            if u_name in metabat_bin:
                unitig_label_metabat[ (path_start, path_key, path_end) ] = metabat_bin[u_name]   # unitigの情報にbinのidを紐付ける
            else:
                unitig_label_metabat[ (path_start, path_key, path_end) ] = bin_id_tmp   # binningされなかったunitigはそれぞれ異なるidを割り振る
                bin_id_tmp -= 1
            """
            if u_name in motif_info:
                unitig_label_motifs[ (path_start, path_key, path_end) ] = motif_info[u_name]   # モチーフの情報も紐付けているが使わない予定
            else:
                unitig_label_motifs[ (path_start, path_key, path_end) ] = set()
            """
    ## 情報付きunitig graphを構成
    binned_edge_data = {}
    binned_graph = nx.MultiDiGraph()   # pathをつなげることによってmulti-edgeが生じる可能性がある
    for path_start, path_key, path_end, path_length, path_score, path, n_edges in c_path:   # c_pathはunitig graphに属するedgeの集合
        #print "in"
        binned_edge_data[ (path_start, path_end, path_key) ] = (path_length, path_score, path)   # edgeの体裁を整えたものに変形

        if (path_start, path_key, path_end) not in unitig_label_metabat:   # simple_nodesによるctg_pathsに含まれていない、すなわちunitigとして出力されなかったedgeが存在することがある
            if path_start == path_end:   # simple loopの場合：unitigとして出力されたものとedgeとでloopの起点の扱いが違うことがあるので、これを解決する。simple loop = circular unitig???
                print "simple loop"
                for s, t, v in path:   # unitig_info中に格納されている同じsimple loopの起点を見つける
                    if (s, t, s) in unitig_label_metabat:
                        print "simple loop resolved"
                        #binned_graph.add_edge( s, s, key = t, bin = unitig_label_metabat[ (s, t, s) ], motif = unitig_label_motifs[ (s, t, s) ])
                        binned_graph.add_edge( s, s, key = t, bin = unitig_label_metabat[ (s, t, s) ])
                        break
                print "[CAUTION] simple loop not resolved!!"   # R contig?
            else:
                print "[CAUTION] not simple loop but cannot detect in unitig_label_metabat!!"
        else:
            #binned_graph.add_edge( path_start, path_end, key = path_key, bin = unitig_label_metabat[ (path_start, path_key, path_end) ], motif = unitig_label_motifs[ (path_start, path_key, path_end) ])   # ラベルを追加してedge multiplexityを許したグラフを再構成
            binned_graph.add_edge( path_start, path_end, key = path_key, bin = unitig_label_metabat[ (path_start, path_key, path_end) ])   # ラベルを追加してedge multiplexityを許したグラフを再構成

    #circular_contigs = []   # unitigをbin情報で繋げたことによって生成された環状コンティグのnx.MultiDiGraph()を格納するリスト
    circular_contig_edges = set()

    ## binとモチーフの情報からunitig graphをentangle、現在の実装は完全に隣り合ったedge同士だけを連結する
    for n in binned_graph.nodes():
        if len(binned_graph.in_edges(n)) >= 1 and len(binned_graph.out_edges(n)) >= 1:   # in_nodeまたはout_nodeどちらかが無い場合は見ない

            bin_list_in = defaultdict(list)   # {bin ID : [現在対象としているin_edgeのうちそのbin IDに属するもの]}
            for e in binned_graph.in_edges(n, keys=True, data=True):
                bin_list_in[e[3]['bin']].append(e)   # e = (path_start, path_end, path_key, {'bin': bin_id, 'motif' : motif_set}), bin_list_in = {bin_id : e_list}

            bin_list_out = defaultdict(list)   # 同じくout_edge
            for e in binned_graph.out_edges(n, keys=True, data=True):
                bin_list_out[e[3]['bin']].append(e)

            print "@node", n
            print "*bin_list_in*"
            for bin, e_list in bin_list_in.items():
                print "bin", bin, ":", e_list
            print "*bin_list_out*"
            for bin, e_list in bin_list_out.items():
                print "bin", bin, ":", e_list

            for bin_id, e_list in bin_list_in.items():   # 各bin IDに対して分離判定を行っており、binのmultiplexityにも対応している
                if len(e_list) == 1 and len(bin_list_out[bin_id]) == 1:   # あるbin IDに属するin_edgeが1本かつout_edgeが1本である時、

                    #print "****", e_list[0], bin_list_out[bin_id][0]

                    ### モチーフによるconnectingの閾値を設けたい場合はここをコメントイン
                    #if jaccard(e_list[0][3]['motif'], bin_list_out[bin_id][0][3]['motif']) < 0.6:   # Jaccard similarityの閾値は適宜変更(パラメタにするべき)
                    #    continue

                    print "[connecting] at node", n, ", bin_id =", bin_id
                    in_e = e_list[0]   # つなげるin_edge
                    in_edge = (in_e[0], in_e[1], in_e[2])   # s, t, v、edgeとして扱うにはこちらの方が都合がよいのでこの形式を取っている
                    out_e = bin_list_out[bin_id][0]   # つなげるout_edge
                    out_edge = (out_e[0], out_e[1], out_e[2])
                    print "[connecting] in_edge =", in_edge, ", out_edge =", out_edge   # これはunitig edge
                    
                    if in_edge == out_edge:   # self-loop (1つのnodeに1つのedgeがin_edgeかつout_edgeとして存在している)である場合は何もしない。本当に何もしなくてよいのか？？？？
                        print "[connecting] this is a self loop. skip."
                        ############################################## path中のself-loopだったら分離せず、in_nodeまたはout_nodeしかなかったら分離？ #########################################################
                        # binによるlabelを通過しているわけだから、無条件に分離してもよいのでは？
                        #binned_graph.remove_edge(*in_edge)
                        #del binned_edge_data[in_edge]
                        #circular_path.add( (s, t, v) )   # 環状コンティグに追加、s, t, vの指定が正しいかは不明、u_edge_data (初期のunitigを縮約する前)に入ってない可能性があるからおそらくダメ？


                        circular_contig_edges.add(in_edge)


                        ## in_edge, out_edgeともに存在する場合はcontinueするなどの条件を付けた方が良い？

                        #in_length, in_score, in_path = binned_edge_data[in_edge]   # binned_edge_dataはc_pathを更新するために必要
                        #out_length, out_score, out_path = binned_edge_data[out_edge]
                        
                        #print in_path, out_path   # これはoverlap edgeのリスト
                        
                        #del binned_edge_data[in_edge]
                        #del binned_edge_data[out_edge]
                        
                        #binned_edge_data[ (in_edge[0] + str(bin_id), in_edge[0] + str(bin_id), in_edge[2]) ] = (in_length, in_score, in_path)   # in_pathは解消前と同じになっているが大丈夫か？？？

                        #binned_graph.remove_edge(*in_edge)
                        #binned_graph.add_node(in_edge[0] + str(bin_id), size = 1)   # ノードを複製してself-loopを生成
                        #binned_graph.add_edge( in_edge[0] + str(bin_id), in_edge[0] + str(bin_id), key = in_edge[2], bin = bin_id )

                        ### わざわざここで環状の判定をしなくても、下でcontigを構成する時にstart_node = end_nodeであれば環状としてしまえばよいのでは？？
                        ### 本当に期待通りに環状となるか調べる

                    else:
                        binned_graph.remove_edge(*in_edge)   # 分離するedgeをグラフから除去
                        binned_graph.remove_edge(*out_edge)
                        
                        binned_graph.add_edge( in_edge[0], out_edge[1], key = in_edge[2], bin = bin_id )   # (in_edgeのs, out_edgeのt, in_edgeのv)というedge、すなわち間のnode nを消して繋げたedgeをグラフに追加
                        
                        in_length, in_score, in_path = binned_edge_data[in_edge]   # binned_edge_dataはc_pathを更新するために必要
                        out_length, out_score, out_path = binned_edge_data[out_edge]
                        
                        #print in_path, out_path   # これはoverlap edgeのリスト
                        
                        del binned_edge_data[in_edge]
                        del binned_edge_data[out_edge]
                        
                        binned_edge_data[ (in_edge[0], out_edge[1], in_edge[2]) ] = (in_length + out_length, in_score + out_score, in_path + out_path)

    print circular_contig_edges

    free_edges = set()
    for s, t, v in ug.edges(keys=True):   # binned_graph上でsimple_nodesでcontigを構成するのならcontig overlapを考慮する必要はないはずなので、free_edges等は不要のはず
        free_edges.add( (s, t, v) )


    c_path = []   # c_pathの更新
    for path_node, path_info in binned_edge_data.items():
        c_path.append( (path_node[0], path_node[2], path_node[1]) + path_info + (len(path_info[2]),) )
#        print c_path[len(c_path) - 1]


    ctg_id = 0

    ctg_paths = open("ctg_paths","w")

    c_path.sort( key=lambda x: -x[3] )

    #print c_path[0]

    cg = nx.DiGraph()
    
    for path_start, path_key, path_end, p_len, p_score, path, n_edges in c_path:
        length = 0
        score = 0
        length_r = 0
        score_r = 0

        non_overlapped_path = []
        non_overlapped_path_r = []
        for s, t, v in path:
            #print s, t, v
            if v != "NA": 
                rs, rt, rv = reverse_end(t), reverse_end(s), reverse_end(v)
            else:
                rs, rt, rv = reverse_end(t), reverse_end(s), "NA"
            if (s, t, v) in free_edges and (rs, rt, rv) in free_edges:
                non_overlapped_path.append( (s,t,v) )
                non_overlapped_path_r.append( (rs, rt, rv)  )
                length += u_edge_data[ (s, t, v) ][0]
                score += u_edge_data[ (s, t, v) ][1]
                length_r += u_edge_data[ (rs, rt, rv) ][0]
                score_r += u_edge_data[ (rs, rt, rv) ][1]
            else:
                break

        #print non_overlapped_path

        if len(non_overlapped_path) == 0:
            #print "rc path?"
            continue
        s0, t0, v0 = non_overlapped_path[0]   # non_overlapped_pathはコンティグを構成するedgeのリスト
        end_node = non_overlapped_path[-1][1]

        cg.add_edge(s0, end_node)   # reverse complement pathは追加しなくてよいか？

        if s0 == end_node:
            print "[INFO] circular edge:", s0, s0, t0

        if s0 == end_node and ((len(binned_graph.in_edges(s0)) == 1 and len(binned_graph.out_edges(s0)) == 1) or (s0, s0, t0) in circular_contig_edges):   # 環状コンティグの判定、完全に環状もしくは解消されたself-loop
            print "[INFO] circular contig detected:", "%06dF" % ctg_id, "#in_edges =", len(binned_graph.in_edges(s0)), "#out_edges =", len(binned_graph.out_edges(s0))
            print >> ctg_paths, "%06dF" % ctg_id, "ctg_circular", s0+"~"+v0+"~"+t0, end_node, length, score, "|".join([ c[0]+"~"+c[2]+"~"+c[1] for c in non_overlapped_path ] )
        else:
            print >> ctg_paths, "%06dF" % ctg_id, "ctg_linear", s0+"~"+v0+"~"+t0, end_node, length, score, "|".join([ c[0]+"~"+c[2]+"~"+c[1] for c in non_overlapped_path ] )
        non_overlapped_path_r.reverse()
        s0, t0, v0 = non_overlapped_path_r[0]
        end_node = non_overlapped_path_r[-1][1]

        cg.add_edge(s0, end_node)   # どうせforwardかreverse complementのcontigしか出力されないのだからこっちは出さなくてよい？

        if s0 == end_node and ((len(binned_graph.in_edges(s0)) == 1 and len(binned_graph.out_edges(s0)) == 1) or (s0, s0, t0) in circular_contig_edges):   # 環状コンティグ、上のものとペアで見つかるはず
            print "[INFO] circular contig (R) detected"
            print >> ctg_paths, "%06dR" % ctg_id, "ctg_circular", s0+"~"+v0+"~"+t0, end_node, length_r, score_r, "|".join([ c[0]+"~"+c[2]+"~"+c[1] for c in non_overlapped_path_r ] )
        else:
            print >> ctg_paths, "%06dR" % ctg_id, "ctg_linear", s0+"~"+v0+"~"+t0, end_node, length_r, score_r, "|".join([ c[0]+"~"+c[2]+"~"+c[1] for c in non_overlapped_path_r ] )
        ctg_id += 1
        for e in non_overlapped_path:
            if e in free_edges:
                free_edges.remove(e)
            else:
                print "[CAUTION] contig path does not exist:", ctg_id
        for e in non_overlapped_path_r:
            if e in free_edges:
                free_edges.remove(e)
            else:
                print "[INFO] contig path (rc) does not exist (maybe F and R are on the same contig?):", ctg_id   # FとRが同じpath上に存在するとき？

    nx.write_gexf(cg, "contig_graph.gexf")

    for s, t, v in list(circular_path):
        length, score, path, type_ = u_edge_data[ (s, t, v) ]
        print >> ctg_paths, "%6d" % ctg_id, "utg_circular", s+"~"+v+"~"+t, t, length, score, s+"~"+v+"~"+t   # 名称をctg_circularからutg_circularに変更
        ctg_id += 1

    ctg_paths.close()

