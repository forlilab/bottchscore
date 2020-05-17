#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright: 2019, Stefano Forli (forli@scripps.edu)

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.


import sys
from collections import OrderedDict
try:
    from openbabel import openbabel as ob
except:
    print("Error: OpenBabel >v.3 is required!")
    sys.exit(1)
from math import log

# ref: https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00723

# TODO 
# - missing atropisomers (https://en.wikipedia.org/wiki/Atropisomer)
# - according to the original implementation, t-butyl and isopropyl have same complexity?

class BottchScore:
    def __init__(self, verbose=False):
        """ 
            TERMS:
            
            d_i :   bonds with different chemical groups
            e_i :   unique list of non-H chemical elements involved in bonds
            s_i :   chirality bit
            v_i :   valence electrons (calculated as octet(8) - max number of bonds)
            b_i :   sum of all bond orders
        
        """
        self.verbose=verbose

    def score(self, mol):
        """ """
        if mol.NumAtoms()==0:
            print("*** WARNING: NO ATOMS FOR MOLECULE [%s] ***" % mol.GetTitle())
            return
        self.mol = mol
        # create the storage for the non-hydrogen atoms
        self._build_automorphism()
        self._indices = OrderedDict()
        self._equivalents = {}
        self._calculate_terms()
        self._calculate_score()
        if self.verbose:
            self.print_table()
        return self._intrinsic_complexity

    def _calculate_terms(self):
        """ calculate the different terms contribution"""
        for idx in range(1, mol.NumAtoms()+1):
            # obabel numbering is 1-based
            atom = mol.GetAtom(idx)
            # hydrogen
            if self._is_hydrogen(idx):
                continue
            self._indices[idx]={}
            self._calc_di(idx, atom)
            self._calc_Vi(idx, atom)
            self._calc_bi_ei_si(idx, atom)
        

    def print_table(self):
        seq = [ 'di', 'ei', 'si', 'Vi', 'bi', 'complexity']
        atom_list = list(self._equivalents.keys())
        print("\t", "\t".join([str(x) for x in atom_list]))
        print("===========================================================")
        for prop in seq:
            if prop=='complexity':
                string = 'cmplx'
            else:
                string = prop
            print("%s\t"% string, end=' ')
            for i in atom_list:
                if prop=='complexity':
                    print("%2.1f\t" % self._indices[i][prop], end=' ')
                else:
                    print("%d\t" % self._indices[i][prop], end=' ')
            print("")
        print('---------------------------------------')
        print("Intrinsic complexity : %2.1f" % self._intrinsic_complexity)
        print("Complexity/atom      : %2.1f" % (self._intrinsic_complexity/len(self._indices)))
        print("===========================================================")
        

    def _calculate_score(self):
        """ calculate total complexity and adjust it to the symmetry factor
            applying eq.3 of ref.1
        """
        total_complexity = 0
        # eq.3 of ref.1, left side
        for idx in list(self._indices.keys()):
            complexity = self._calculate_complexity(idx)
            self._indices[idx]['complexity'] = complexity
            total_complexity += complexity
        # eq.3 of ref.1, left side
        for idx, eq_groups in list(self._equivalents.items()):
            for e in eq_groups:
                total_complexity -= 0.5 *self._indices[idx]['complexity'] / (len(eq_groups))
        self._intrinsic_complexity = total_complexity 

    def _calculate_complexity(self,idx):
        """ perform calculation of single complexity on a given index"""
        data = self._indices[idx]
        if self.verbose:
            print(("complexty [%3d]:  %d * %d * %d log2( %d * %d) = %2.1f"  % (idx, 
                data['di'], data['ei'],data['si'], data['Vi'], data['bi'], 
                data['di']*data['ei']*data['si']*log(data['Vi']*data['bi'], 2)
                )))

        return data['di']*data['ei']*data['si']*log(data['Vi']*data['bi'], 2)


    def _calc_Vi(self, idx, atom):
        """ """
        num = atom.GetAtomicNum()
        val = 8-ob.GetMaxBonds(num)
        self._indices[idx]['Vi'] = val


    def _calc_bi_ei_si(self, idx, atom):
        """ """
        bi = 0
        ei = [atom.GetAtomicNum()]
        #this_element = atom.GetAtomicNum()
        si = 1
        if atom.IsChiral():
            si+=1
        for neigh in ob.OBAtomAtomIter(atom):
            if self._is_hydrogen(neigh.GetIdx()):
                continue
            bond = self.mol.GetBond(atom, neigh)
            # get bond order
            bi += bond.GetBondOrder()
            # get neighbor element
            ei.append(neigh.GetAtomicNum())
        self._indices[idx]['bi'] = bi
        self._indices[idx]['ei'] = len(set(ei))
        self._indices[idx]['si'] = si

    def _calc_di(self, idx, atom):
        """ """
        # keep track of equivalent groups, and mark the first one
        # to be used as main
        self._equivalents[idx] = []
        if not idx in self._equivalents:
            if idx in self.automorphs:
                if len(set(self._equivalents.keys()) & self.automorphs[idx])==0:
                    self._equivalents[idx]=self.automorphs[idx]
            else:
                self._equivalents[idx]=set()
            
            #self._equivalents[idx] = []
        if idx in self.automorphs:
            for u in self.automorphs[idx]:
                self._equivalents[idx].append(u)
        # count how many non-equivalent neighbors there are for atom
        groups = []
        for neigh in ob.OBAtomAtomIter(atom):
            if self._is_hydrogen(neigh.GetIdx()):
                continue
            neigh_idx = neigh.GetIdx()
            if (neigh_idx in self.automorphs):
                if len(set(groups) & self.automorphs[neigh_idx])>0:
                    continue
            groups.append(neigh_idx)
        self._indices[idx]['di'] = len(groups)


    def _build_automorphism(self):
        """ automorphisms are 0-based
            prune the automorphism map to remove
            identities and compact multiple
            mappings
        """
        automorphs = ob.vvpairUIntUInt()
        ob.FindAutomorphisms(self.mol, automorphs)
        self.automorphs = {}
        for am in automorphs:
            for i,j in am:
                if i==j:
                    continue
                k=i+1
                l=j+1
                if self._is_hydrogen(k):
                    continue
                if self._is_hydrogen(l):
                    continue
                if not k in self.automorphs:
                    self.automorphs[k] = []
                self.automorphs[k].append(l)
        for k,v in list(self.automorphs.items()):
            self.automorphs[k] = set(v)
        if self.verbose:
            from pprint import pprint
            print("AUTOMORPHS")
            pprint(self.automorphs)

    def _is_hydrogen(self, idx):
        return self.mol.GetAtom(idx).GetAtomicNum()==1


if __name__=='__main__':
    import sys
    import os
    import argparse
    usage = "Usage: %s -i molfile.ext [-p] [-v]" % sys.argv[0]
    epilog = ""
    parser = argparse.ArgumentParser(description='Tool to calculate Böttcher score on small molecules.', usage=usage, 
                epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', action="store", help='input structure; support all input formats supported by OB, including multi-structure formats', metavar='filename.ext', required=True)
    parser.add_argument('-p', action="store_true", help='generate PNG image of the structure', default=False)
    parser.add_argument('-v', action="store_true", help='verbose mode; print the full table of the terms used to estimate the score, as described in the paper', default=False)
    if len(sys.argv)<2:
        parser.print_help()
        sys.exit(1)
    ARGS = parser.parse_args()

    # input file
    infile = ARGS.i
    # check if verbose
    verbose = ARGS.v
    # save image
    save_png = ARGS.p
    # Parse format
    name, ext = os.path.splitext(infile)
    ext = ext[1:].lower()
    counter=1
    # initialize mol parser
    mol_parser = ob.OBConversion()
    mol_parser.SetInAndOutFormats(ext,'smi')
    # initialize class
    bottch = BottchScore(verbose)
    # load the first molecule found in the file
    mol = ob.OBMol()
    more = mol_parser.ReadFile(mol, infile)
    # score the molecule
    score=bottch.score(mol)
    if not verbose:
        print("%d: %2.1f"% (counter,score))
    if save_png:
        name = mol.GetTitle()
        if name.strip()=="":
            name = "mol_%d" % counter
        mol_parser.SetOutFormat('png')
        mol_parser.WriteFile(mol, '%s_image.png' % name)
    # check if there are more molecules in the file and score them
    while more:
        mol = ob.OBMol()
        more = mol_parser.Read(mol)
        if not more:
            break
        counter+=1
        score = bottch.score(mol)
        if not verbose:
            print("%d: %2.1f"% (counter,score))
        if save_png:
            name = mol.GetTitle()
            if name.strip()=="":
                name = "mol_%d" % counter
            mol_parser.SetOutFormat('png')
            mol_parser.WriteFile(mol, '%s_image.png' % name)



