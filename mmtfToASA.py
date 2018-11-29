#!/user/bin/env python
'''mmtfToASA.py

Creates a dataset of accesible surface areas. The input to this class must be a single
protein chain.

Examples
--------
get dataset of phi/psi angles:

>>> pdb.flatMapToPair(new StructureToPolymerChains())
...    .filter(new ContainsLProteinChain())
>>> asa = mmtfToASA.getDataset(pdb)
>>> asa.show(10)

'''

from mmtfPyspark.ml import pythonRDDToDataset
from pyspark.sql import Row

import freesasa

def get_dataset(structure, parameters=None, classifier=None, options=None):
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

    rows = structure.map(
        lambda x: _get_free_sasa(x, parameters, classifier, options))  # Map or flatMap

    # convert to dataset
    colNames = ["structureChainId", "totalArea"]

    return pythonRDDToDataset.get_dataset(rows, colNames)


def get_python_rdd(structure):
    '''Returns a pythonRDD of 3-state secondary structure

    Parameters
    ----------
    structure : mmtfStructure
    '''

    return structure.map(lambda x: _get_phi_psi(x))


def _get_free_sasa(t, parameters=None, classifier=None, options=None):
    '''Get factions of alpha, beta and coil within a chain
    '''

    key = t[0]
    structure = t[1]
    if structure.num_chains != 1:
        raise Exception(
            "This method can only be applied to single polyer chain.")

    dsspQ8, dsspQ3 = '', ''

    groupIndex = 0
    atomIndex = 0

    freesasaStructure = freesasa.Structure()
    if (classifier is None):
        classifier = freesasa.Classifier()
    optbitfield = freesasa.Structure._get_structure_options(options or freesasa.Structure.defaultOptions)

    for i in range(0, structure.num_models):
        print("model: " + str(i+1))

        for j in range(0, structure.chains_per_model[i]):
            chainName = structure.chain_name_list[chainIndex]
            chainId = structure.chain_id_list[chainIndex]
            groups = structure.groups_per_chain[chainIndex]

            entityType = structure.entity_list[chainToEntityIndex[chainIndex]]["type"]

            #if not entityType == "polymer": continue

            prev_coords = None
            coords = None
            for k in range(0, structure.groups_per_chain[chainIndex]):
                groupId = structure.group_id_list[groupIndex]
                insertionCode = structure.ins_code_list[groupIndex]
                secStruct = structure.sec_struct_list[groupIndex]
                seqIndex = structure.sequence_index_list[groupIndex]

                groupType = structure.group_type_list[groupIndex]
                groupName = structure.group_list[groupType]["groupName"]

                for i, name in enumerate(structure.group_list[groupType]["atomNameList"]):
                    if (classifier.classify(groupName, name) is 'Unknown'):
                        if (optbitfield & freesasa.FREESASA_SKIP_UNKNOWN):
                            continue
                        if (optbitfield & freesasa.FREESASA_HALT_AT_UNKNOWN):
                            raise Exception("Halting at unknown atom")

                    freesasaStructure.addAtom(name, groupName, seqIndex, chainName,
                                      structure.x_coord_list[atomIndex+i],
                                      structure.y_coord_list[atomIndex+i],
                                      structure.z_coord_list[atomIndex+i])

                atomIndex += len(structure.group_list[groupType]["atomNameList"])
                groupIndex += 1

    freesasaStructure.setRadiiWithClassifier(classifier)
    freesasaResult = freesasa.calc(freesasaStructure, parameters or freesasa.Parameters())
    sasa_classes = classifyResults(freesasaResult, freesasaStructure, classifier)

    return Row([key, sasa_classes.totalArea])
