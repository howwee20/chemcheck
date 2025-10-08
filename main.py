# main.py
import io, logging
from functools import lru_cache

import requests
import pubchempy as pcp
from fastapi import FastAPI, Form, Query
from fastapi.responses import StreamingResponse, JSONResponse, HTMLResponse, PlainTextResponse
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

APP_VERSION = "0.5.2"

# ----- app & logging -----
logging.basicConfig(level=logging.INFO)
log = logging.getLogger("chemcheck")
app = FastAPI(title="ChemCheck Engine", version=APP_VERSION)

# ---------- name → SMILES resolver ----------
@lru_cache(maxsize=1024)
def resolve_to_smiles(name_or_smiles: str) -> str:
    """
    Resolve a common name or SMILES into a valid canonical SMILES.
    Order: direct SMILES → OPSIN → PubChem.
    """
    s = (name_or_smiles or "").strip()
    if not s:
        raise ValueError("Empty input")

    # 1) Already a SMILES?
    mol = Chem.MolFromSmiles(s)
    if mol:
        can = Chem.MolToSmiles(mol)
        log.info("Resolver: input is already SMILES -> %s", can)
        return can

    # 2) OPSIN (name → SMILES)
    try:
        url = f"https://opsin.ch.cam.ac.uk/opsin/{requests.utils.quote(s)}.json"
        r = requests.get(url, timeout=5)
        if r.ok:
            data = r.json()
            smi = data.get("smiles")
            if smi and Chem.MolFromSmiles(smi):
                can = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
                log.info("Resolver: OPSIN hit for %r -> %s", s, can)
                return can
        else:
            log.info("Resolver: OPSIN no hit (%s) for %r", r.status_code, s)
    except Exception as e:
        log.info("Resolver: OPSIN failed for %r: %s", s, e)

    # 3) PubChem (fallback)
    try:
        comps = pcp.get_compounds(s, "name")
        if comps:
            smi = getattr(comps[0], "connectivity_smiles", None) or getattr(comps[0], "canonical_smiles", None)
            if smi and Chem.MolFromSmiles(smi):
                can = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
                log.info("Resolver: PubChem hit for %r -> %s", s, can)
                return can
        log.info("Resolver: PubChem empty for %r", s)
    except Exception as e:
        log.info("Resolver: PubChem failed for %r: %s", s, e)

    raise ValueError(f"Could not resolve '{s}' to a structure")

# ---------- renderer ----------
def render_smiles(
    smiles: str,
    w: int = 500,
    h: int = 400,
    fmt: str = "png",
    show_indices: bool = False,
    show_lone_pairs: bool = False,
    show_formal_charges: bool = True,
) -> bytes:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES passed to renderer")

    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(w, h) if fmt == "svg" else Draw.MolDraw2DCairo(w, h)
    dopts = drawer.drawOptions()

    # Safe options across RDKit versions
    if hasattr(dopts, "addAtomIndices"):
        dopts.addAtomIndices = bool(show_indices)
    if hasattr(dopts, "addStereoAnnotation"):
        dopts.addStereoAnnotation = True
    if hasattr(dopts, "explicitMethyl"):
        dopts.explicitMethyl = True
    if hasattr(dopts, "addAtomCharge"):
        dopts.addAtomCharge = bool(show_formal_charges)
    # ✅ correct name for the padding knob:
    if hasattr(dopts, "additionalAtomLabelPadding"):
        dopts.additionalAtomLabelPadding = 0.10

    # crude lone-pair hint: annotate O/N/S with ":" label (visual cue only)
    if show_lone_pairs:
        for a in mol.GetAtoms():
            if a.GetSymbol() in ("O", "N", "S"):
                a.SetProp("atomLabel", f"{a.GetSymbol()}:")
        # only set atomLabels if this RDKit has the attribute
        if hasattr(dopts, "atomLabels"):
            dopts.atomLabels = {
                i: a.GetProp("atomLabel")
                for i, a in enumerate(mol.GetAtoms())
                if a.HasProp("atomLabel")
            }

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    data = drawer.GetDrawingText()
    return data.encode("utf-8") if fmt == "svg" else data

# ---------- routes ----------
@app.get("/", response_class=HTMLResponse)
def root():
    return """
<!doctype html>
<html>
  <body style="font-family: system-ui; margin:40px auto; max-width:720px; line-height:1.5;">
    <h1>ChemCheck</h1>
    <p>Instantly generate clean 2D chemical structures from names or SMILES.</p>
    <p><a href="/ui" style="display:inline-block;padding:10px 16px;border:1px solid #222;border-radius:8px;text-decoration:none;">Launch Validator →</a></p>
    <p><small>Free during beta • For learning/validation only • Built by an MSU student</small></p>
  </body>
</html>
"""

@app.get("/__status", response_class=PlainTextResponse)
def status():
    return f"ChemCheck {APP_VERSION}"

@app.get("/test-resolve", response_class=PlainTextResponse)
def test_resolve(q: str = Query(..., description="name or SMILES")):
    return resolve_to_smiles(q)

@app.post("/draw")
def draw(
    name_or_smiles: str = Form(...),
    fmt: str = Form("png"),
    show_indices: bool = Form(False),
    show_lone_pairs: bool = Form(False),
    show_formal_charges: bool = Form(True),
    w: int = Form(500),
    h: int = Form(400),
):
    try:
        # Resolve name → SMILES
        try:
            smiles = resolve_to_smiles(name_or_smiles)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, status_code=400)

        # Render
        try:
            blob = render_smiles(
                smiles, w=w, h=h, fmt=fmt,
                show_indices=show_indices,
                show_lone_pairs=show_lone_pairs,
                show_formal_charges=show_formal_charges,
            )
        except Exception as e:
            log.info("Renderer failed for %r (%s): %s", name_or_smiles, smiles, e)
            return JSONResponse({"error": "Failed to render structure."}, status_code=400)

        media = "image/svg+xml" if fmt == "svg" else "image/png"
        return StreamingResponse(io.BytesIO(blob), media_type=media)
    except Exception as e:
        log.info("draw endpoint unexpected failure for %r: %s", name_or_smiles, e)
        return JSONResponse(
            {"error": "Invalid input. Try a common name (e.g., caffeine) or a valid SMILES."},
            status_code=400,
        )

@app.get("/ui", response_class=HTMLResponse)
def ui():
    return """
<!doctype html>
<html>
  <body style="font-family: system-ui; max-width:760px; margin:40px auto;">
    <h1>ChemCheck</h1>
    <form id="f" style="display:flex; gap:8px; align-items:center; flex-wrap:wrap;">
      <input name="name_or_smiles" placeholder="Name or SMILES" style="flex:1; min-width:360px; padding:8px;" />
      <select name="fmt"><option>png</option><option>svg</option></select>
      <label><input type="checkbox" name="show_indices"> indices</label>
      <label><input type="checkbox" name="show_lone_pairs"> lone pairs</label>
      <label><input type="checkbox" name="show_formal_charges" checked> charges</label>
      <input name="w" type="number" value="500" style="width:90px;">w
      <input name="h" type="number" value="400" style="width:90px;">h
      <button type="submit">Generate</button>
    </form>
    <div id="out" style="margin-top:16px;"></div>
    <div id="dl" style="margin-top:8px; display:none;">
      <a id="save" download="chemcheck.png">Download</a>
    </div>
    <script>
      const f=document.getElementById('f'),out=document.getElementById('out'),dl=document.getElementById('dl'),save=document.getElementById('save');
      f.onsubmit=async e=>{
        e.preventDefault();
        out.innerHTML='Rendering...'; dl.style.display='none';
        const fd=new FormData(f);
        const resp=await fetch('/draw',{method:'POST',body:fd});
        if(!resp.ok){ out.innerHTML='<div style="color:#b00;">'+ await resp.text() +'</div>'; return; }
        const blob=await resp.blob();
        const url=URL.createObjectURL(blob);
        const fmt=fd.get('fmt')||'png';
        out.innerHTML = fmt==='svg'
          ? await blob.text().then(t=>t)
          : `<img src="${url}" style="max-width:100%;">`;
        save.href=url; save.download=`chemcheck.${fmt}`; dl.style.display='block';
      };
    </script>
    <small>For learning/validation only. Always do your own work.</small>
  </body>
</html>
"""
