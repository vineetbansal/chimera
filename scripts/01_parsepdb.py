from prody import parsePDB, flagDefinition
from prody.atomic.hierview import HierView
from prody.atomic.residue import Residue


if __name__ == '__main__':

    pdbfile = '/media/vineetb/t5-vineetb/biolip/downloaded_data/receptor/2lueA.pdb'
    atoms = parsePDB(pdbfile)

    print(atoms.numAtoms())  # 1981
    print(flagDefinition('backbone'))
    sel = atoms.select('stdaa')  # standard amino acids
    print(sel.numAtoms())  # 1981

    sel = atoms.select('backbone')  # N CA C O
    print(sel.numAtoms())  # 119 * 4 = 476

    sel = atoms.select('not backbone')
    print(sel.numAtoms())  # 1981 - 476 = 1505

    sel = atoms.select('not hydrogen')
    print(sel.numAtoms())  # 981

    sel = atoms.select('stdaa and not backbone and not hydrogen')
    print(sel.numAtoms())  # 505
    print('----------------')

    # The 505 we have selected correspond to 113 different residues
    for i, res in enumerate(HierView(sel).iterResidues()):
        print(i, res.getResname())
        coords = res._getCoords()

    # ----------------------------------------- #

    ligand_id = 'III'
    pdbfile = '/media/vineetb/t5-vineetb/biolip/downloaded_data/ligand/2lue_III_B_1.pdb'

    print('debug')

