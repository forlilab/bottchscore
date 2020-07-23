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
# - according to the original implementation, t-butyl and isopropyl have sam complexiy?


class BottchScore:
    def __init__(self, verbose=False, debug=False):
        """ 
            TERMS:
            
            d_i :   bonds with different chemical groups
            e_i :   unique list of non-H chemical elements involved in bonds
            s_i :   chirality bit
            v_i :   valence electrons (calculated as octet(8) - max number of bonds)
            b_i :   sum of all bond orders
        
        """
        self.converter = ob.OBConversion()
        self.converter.SetOutFormat('smi')
        self.verbose=verbose
        self.debug=debug
        self._mesomery_patterns = {
                # SMARTS_pattern : [equivalent atoms idx list, contribution ]
                '[$([#8;X1])]=*-[$([#8;X1])]' : [ [[0,2]], 1.5], # carboxylate, nitrate
                #'[$([#7;X2](=*)(-*=*-*))]': [ [[0,1],[0,2],[2,3],[3,4]], 1.5 ],
                '[$([#7;X2])](=*)(-*=[$(*-*)])': [ [[0,1],[0,2],[2,3]], 1.5 ],


                #'[$([#8;X1])]~*~[$([#8;X1])]' : [ [0,2], 1.5], # carboxylate, nitrate
                #'[$([#8;X1])]~*~[$([#8;X2,H1])]' : [ [0,2], 1.5], # carboxylate, nitrate
                #'[$([#7;X2](=*)[*=*])]'
                # '[$([#7;X2](=*)(-*=*(-*)))]' :  [ [[]],  1.5 ], # guanidine, azire

                    
                #'[CH2X4]' # histidine
                '[#6X3]1:' # imidazole
                '[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:'
                '[#6X3H]:'
                '[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:'
                '[#6X3H]1' :  [[[1,3]], 2.5],
                }

    def score(self, mol, disable_mesomeric=False):
        """ """
        if self.debug:
            print("==========================================")
        if mol.NumAtoms()<2:
            #print("*** WARNING: NO ATOMS FOR MOLECULE [%s] ***" % mol.GetTitle())
            return 0
            #return
        self._initialize_mol(mol, disable_mesomeric)
        self._calculate_terms()
        self._calculate_score()
        if self.verbose:
            self.print_table()
        return self._intrinsic_complexity

    def _initialize_mol(self, mol, disable_mesomeric):
        """ perform initial operations on the molecule caching chemical
        information for aromatic and some tautomeric perceptions"""
        self.mol = mol
        #if self.verbose:
        full = self.converter.WriteString(self.mol)
        self._smiles = full.split()[0]
        if self.debug: print("DEBUG> SMILES: %s" % self._smiles)
        # create the storage for the non-hydrogen atoms
        self._build_automorphism()
        self._indices = OrderedDict()
        self._mesomery_equivalence = {}
        if not disable_mesomeric:
            self._calc_mesomery()

    def _calc_mesomery(self):
        """calculate atoms for which b_i value needs to be corrected"""
        matcher = ob.OBSmartsPattern()
        for patt, idx_info in self._mesomery_patterns.items():
            #print("IDXINFO", idx_info)
            idx_pairs, contribution = idx_info
            matcher.Init(patt)
            found = matcher.Match(self.mol)
            if not found:
                #print ("NADA")
                continue
            if self.debug: print("DEBUG> Matched pattern |%s|"% patt)
            found = [list(x) for x in matcher.GetUMapList()]
            for f in found:
                if self.debug: print("DEBUG> found pattern: ",f)
                for pair in idx_pairs:
                    for idx in pair:
                        if self.debug: print("DEBUG> Assigning mesomery:",f, "->", f[idx])
                        self._mesomery_equivalence[f[idx]] = contribution
                    p0 = f[pair[0]]
                    p1 = f[pair[1]]
                    if self.debug: print("DEBUG> Updating automorphs:", p0,p1 )
                    if not p0 in self.automorphs:
                        self.automorphs[p0] = set()
                    self.automorphs[p0].add(p1)
                    if not p1 in self.automorphs:
                        self.automorphs[p1] = set()
                    self.automorphs[p1].add(p0)
        if self.debug:
            print("DEBUG> MESO-EQUIV", self._mesomery_equivalence)
            print("DEBUG> AUTO-EQUIV", self.automorphs)

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
        print("Intrinsic complexity : %2.2f" % self._intrinsic_complexity)
        print("Complexity/atom      : %2.2f" % (self._intrinsic_complexity/len(self._indices)))
        print("===========================================================")
        

    def _calculate_score(self):
        """ calculate total complexity and adjust it to the symmetry factor
            applying eq.3 of ref.1
        """
        total_complexity = 0
        # eq.3 of ref.1, left side
        if self.verbose:
            print(">> %s" % self._smiles)
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
            try:
                print(("complexity [%3d]:  %d * %d * %d log2( %d * %d) = %2.1f"  % (idx, 
                    data['di'], data['ei'],data['si'], data['Vi'], data['bi'], 
                    data['di']*data['ei']*data['si']*log(data['Vi']*data['bi'], 2)
                    )))
            except:
                print("[ *** Error calculating complexity: atom_idx[%d] *** ]" % idx)
        try:
            return data['di']*data['ei']*data['si']*log(data['Vi']*data['bi'], 2)
        except:
            print("[ *** Error calculating complexity: atom_idx[%d] *** ]" % idx)
            return 0

    def _calc_Vi(self, idx, atom):
        """ """
        valence = atom.GetTotalValence()
        charge = atom.GetFormalCharge()
        if self.debug: print("DEBUG> Vi[%3d]: v = %d | charge = %d"% (idx, valence, charge))
        vi = 8 - valence + charge
        self._indices[idx]['Vi'] = vi

    def _calc_bi_ei_si(self, idx, atom):
        """ """
        bi = 0
        ei = [self._get_atomic_num(atom)]
        atom_idx = atom.GetIdx()
        si = 1
        if atom.IsChiral():
            si+=1
        for neigh in ob.OBAtomAtomIter(atom):
            neigh_idx = neigh.GetIdx()
            if self._is_hydrogen(neigh_idx):
                continue
            if atom_idx in self._mesomery_equivalence:
                # use mesomeric-corrected bond order
                #print("FOUNX", atom_idx, atom_idx in self._mesomery_equivalence, self._mesomery_equivalence)
                contribution = self._mesomery_equivalence[atom_idx]
            else:
                # get bond order
                bond = self.mol.GetBond(atom, neigh)
                contribution = bond.GetBondOrder()
            #print("BONDORDER B_i", contribution)
            bi += contribution
            # get neighbor element
            ei.append(self._get_atomic_num(neigh))
        self._indices[idx]['bi'] = bi
        self._indices[idx]['ei'] = len(set(ei))
        self._indices[idx]['si'] = si

    def _get_atomic_num(self, atom):
        """ function to return the atomic number taking into account 
            isotope
        """
        isotope = atom.GetIsotope()
        if not isotope==0:
            return isotope
        return atom.GetAtomicNum()
        

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
        #itab = ob.OBIsitopeTable()
        self._equivalents = {}
        automorphs = ob.vvpairUIntUInt()
        mol_copy = ob.OBMol(self.mol)
        for i in ob.OBMolAtomIter(mol_copy):
            isotope = i.GetIsotope()
            if not isotope==0:
                n = i.GetAtomicNum()
                i.SetAtomicNum(n+isotope)
        ob.FindAutomorphisms(mol_copy, automorphs)
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
        if self.debug:
            from pprint import pprint
            print("AUTOMORPHS", end='')
            pprint(self.automorphs)
        #elif self.debug:
        #    print("DEBUG> Automorphs:", self.automorphs)

    def _is_hydrogen(self, idx):
        return self.mol.GetAtom(idx).GetAtomicNum()==1


if __name__=='__main__':
    import sys
    import os
    import argparse
    debug = False
    usage = "%s -i molfile.ext [-p] [-v] " % sys.argv[0]
    epilog = ""
    parser = argparse.ArgumentParser(description='Tool to calculate Böttcher score on small molecules.', usage=None, 
                epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', action="store", help='input structure; support all input formats supported by OB, including multi-structure formats', metavar='filename.ext', required=True)
    parser.add_argument('-m', action="store_true", help='disable mesomeric effect estimate', default=False)
    parser.add_argument('-p', action="store_true", help='generate PNG image of the structure', default=False)
    parser.add_argument('-v', action="store_true", help='verbose mode; print the full table of the terms used to estimate the score, as described in the paper', default=False)
    #parser.add_argument('-s', action="store", help='SMILES input mode')
    if len(sys.argv)<2:
        parser.print_help()
        sys.exit(1)
    if "-d" in sys.argv:
        debug=True
        sys.argv.remove('-d')
    ARGS = parser.parse_args()
    # input file
    infile = ARGS.i
    # check if verbose
    verbose = ARGS.v
    # save image
    save_png = ARGS.p
    # mesomeric effect
    disable_mesomeric = ARGS.m
    # Parse format
    name, ext = os.path.splitext(infile)
    ext = ext[1:].lower()
    counter=1
    # initialize mol parser
    mol_parser = ob.OBConversion()
    mol_parser.SetInAndOutFormats(ext,'smi')
    # initialize class
    bottch = BottchScore(verbose, debug)
    # load the first molecule found in the file
    mol = ob.OBMol()
    more = mol_parser.ReadFile(mol, infile)
    while more:
        ob.PerceiveStereo(mol)
        # score the molecule
        score=bottch.score(mol, disable_mesomeric)
        if not verbose:
            print("%d: %4.2f\t %s"% (counter,score,mol.GetTitle()))
        if save_png:
            name = mol.GetTitle()
            if name.strip()=="":
                name = "mol_%d" % counter
            mol_parser.SetOutFormat('png')
            mol_parser.WriteFile(mol, '%s_image.png' % name)
        more = mol_parser.Read(mol)
        counter+=1
