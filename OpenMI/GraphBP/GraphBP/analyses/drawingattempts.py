from rdkit import Chem
# from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
# IPythonConsole.ipython_useSVG=True  #< set this to False if you want PNGs instead of SVGs
# IPythonConsole.molSize = 500,500



mol = Chem.MolFromSmiles("O=C(NC1=CNN/C1=C1/N=c2ccc(CN3CCOCC3)cc2=N1)NC1CC1")
# mol = Chem.AddHs(mol)
mol

# sppl = Chem.SDMolSupplier(mol)
# outname = mol.replace(".sdf", ".txt")
# out_file = open(outname, "w")
# out_file.write(f"{smi}\t{name}\n")
# out_file.close()

writer = Chem.SDWriter('newmol.sdf')  #(os.path.join(f'{write_dir}', 'Hmol.sdf'))
writer.write(mol)
writer.close()


img = Draw.MolsToGridImage(mol)

# with open('molgrid.png', 'wb') as png:
#     png.write(img.data)