#!/user/bin/env python
'''phiPsiExtractor.py

Creates a dataset of Phi/Psi Angles. The input to this class must be a single
protein chain.

Examples
--------
get dataset of phi/psi angles:

>>> pdb.flatMapToPair(new StructureToPolymerChains())
...    .filter(new ContainsLProteinChain())
>>> phiPsi = phiPsiExtractor.getDataset(pdb)
>>> PhiPsi.show(10)

'''
from pyspark import SparkContext
from mmtfPyspark.ml import pythonRDDToDataset
from mmtfPyspark.utils import DsspSecondaryStructure, mmtfStructure
from pyspark.sql import Row
import numpy as np
import itertools as it
from pyspark.sql import DataFrame
from functools import reduce

from Bio.PDB.Polypeptide import three_to_index, aa3

def get_dihedral(p0, p1, p2, p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def get_dataset(structure):
    '''Returns a dataset with protein sequence and secondary structure assignments.

    Parameters
    ----------
    structure : mmtfStructure
       single protein chain

    Returns
    -------
    dataset
       dataset with sequence and secondary structure assignments
    '''
    print("RUNNING")
    rows = structure.flatMap(
        lambda x: _get_phi_psi(x))  # Map or flatMap
    print("MAPPED")
    # convert to dataset
    colNames = ["pdbId", "chain", "resi", "resn", "phi", "psi"]#+["is"+aa for aa in aa3]
    #sc = SparkContext.getOrCreate()

    #newdf = it.chain.from_iterable(rows)
    #allrows = sc.union(rows)
    #return reduce(DataFrame.unionAll, rows.collect())
    return pythonRDDToDataset.get_dataset(rows, colNames)


def get_python_rdd(structure):
    '''Returns a pythonRDD of 3-state secondary structure

    Parameters
    ----------
    structure : mmtfStructure
    '''

    return structure.map(lambda x: _get_phi_psi(x))


def _get_phi_psi(t):
    '''Get Phi and Psi angles per residue from a chain
    '''
#     with open("/home/ec2-user/types.log", "a") as f:
#         f.write("type: {0}\n".format(type(t[1])))
    
    print("RUNNING", t[0])
    key = t[0]
    structure = mmtfStructure.MmtfStructure(t[1].encode_data())
    #structure = t[1]
#     with open("/home/ec2-user/types.log", "a") as f:
#         f.write("new type: {0}\n".format(type(structure)))
    
    
    if structure.num_chains != 1:
        raise Exception(
            "This method can only be applied to single polyer chain.")

#     with open("/home/ec2-user/types.log", "a") as f:
#         f.write("num_chains: {0}\n".format(structure.num_chains))
        
    _pdbId, _chainId = key.split(".")
    _pdbId = _pdbId.lower()
    
#     with open("/home/ec2-user/types.log", "a") as f:
#         f.write("pdbId: {0}\n".format(_pdbId))

    chainIndex = 0
    groupIndex = 0
    atomIndex = 0
    torsion = []
    
#     with open("/home/ec2-user/types.log", "a") as f:
#         f.write("num_models: {0}\n".format(structure.num_models))
    
#    foo = 0
    for i in range(0, structure.num_models):
#        if foo == 0:
        #print("model: " + str(i+1))
#            foo = 1

        with open("/home/ec2-user/types.log", "a") as f:
            f.write("model: {0}\n".format(str(i+1)))
        for j in range(0, structure.chains_per_model[i]):
            chainName = structure.chain_name_list[chainIndex]
            chainId = structure.chain_id_list[chainIndex]
            groups = structure.groups_per_chain[chainIndex]

#                 with open("/home/ec2-user/types.log", "a") as f:
#                     f.write("chainName: {0}\n".format(chainName))
#                     f.write("chainId: {0}\n".format(chainId))

            prev_coords = None
            #coords = None
            for k in range(0, structure.groups_per_chain[chainIndex]):
                groupId = structure.group_id_list[groupIndex]
                #insertionCode = structure.ins_code_list[groupIndex]
                secStruct = structure.sec_struct_list[groupIndex]
                seqIndex = structure.sequence_index_list[groupIndex]
                groupType = structure.group_type_list[groupIndex]
                groupName = structure.group_list[groupType]["groupName"]

#                     with open("/home/ec2-user/types.log", "a") as f:
#                         f.write("groupId: {0}\n".format(groupId))
#                         f.write("groupType: {0}\n".format(groupType))
#                         f.write("groupName: {0}\n".format(groupName))
#                         #f.write("insertionCode: {0}\n".format(insertionCode))

                coords = {name.strip(): np.array((
                    structure.x_coord_list[atomIndex+i],
                    structure.y_coord_list[atomIndex+i],
                    structure.z_coord_list[atomIndex+i])) \
                        for i, name in enumerate(structure.group_list[groupType]["atomNameList"])}

#                     with open("/home/ec2-user/types.log", "a") as f:
#                         f.write("coords was None, now is: {0}\n".format(coords))

                atomIndex += len(structure.group_list[groupType]["atomNameList"])

#                     with open("/home/ec2-user/types.log", "a") as f:
#                         f.write("atomIndex: {0}\n".format(atomIndex))

                try:
                    next_group_type = structure.group_type_list[groupIndex+1]
                    next_coords = {name.strip(): np.array((
                        structure.x_coord_list[atomIndex+i],
                        structure.y_coord_list[atomIndex+i],
                        structure.z_coord_list[atomIndex+i])) \
                            for i, name in enumerate(structure.group_list[next_group_type]["atomNameList"])} \
                    if groupIndex < structure.num_groups-1 else None
                except IndexError:
                    next_coords = None

                try:
                    if prev_coords:

#                             with open("/home/ec2-user/types.log", "a") as f:
#                                 f.write("prev_coords was true, prev_coords: {0}\n".format(prev_coords))

                        phi = get_dihedral(prev_coords["C"], coords["N"], coords["CA"], coords["C"])  
                    else:
                        print(prev_coords)
                        phi = np.NaN
                    #phi = get_dihedral(prev_coords["C"], coords["N"], coords["CA"], coords["C"]) if prev_coords else np.NaN
                except KeyError:
                    phi = np.NaN
                try:
                    psi = get_dihedral(coords["N"], coords["CA"], coords["C"], next_coords["N"]) if next_coords else np.NaN
                except KeyError:
                    psi = np.NaN

                prev_coords = coords

                feats = Row(_pdbId, _chainId, int(seqIndex), str(groupName), float(phi), float(psi))#+get_residue(groupName)
                torsion.append(feats)
                groupIndex += 1

#                     with open("/home/ec2-user/types.log", "a") as f:
#                         #f.write("finally, torsion is: {0}\n".format(str(torsion)))
#                         f.write("finally, coords is: {0}\n".format(coords))

            chainIndex += 1

    return torsion

one_hot = np.eye(20)
def get_residue(res):
    if res == "MSE": res="MET"
    try:
        index = three_to_index(res)
        mask = [int(i==index) for i in range(20)]
        return mask
    except (KeyError, IndexError):
        return [0]*20
