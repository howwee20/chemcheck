# engine.py
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def smiles_to_png(smiles: str, w=300, h=300) -> bytes:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    AllChem.Compute2DCoords(mol)
    drawer = Draw.MolDraw2DCairo(w, h)
    Draw.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

if __name__ == "__main__":
    png = smiles_to_png("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
    with open("aspirin.png", "wb") as f:
        f.write(png)
    print("Wrote aspirin.png")
