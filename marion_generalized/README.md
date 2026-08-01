# A one-parameter generalization of Marion Walter's theorem

Marion Walter's theorem: the central hexagon cut from a triangle by the six cevians to the two
**trisection** points of each side has area 1/10 of the triangle. Replacing the trisection points
by the two points dividing each side in ratios `t:1` and `1:t` gives a one-parameter family
(Marion = `t = 2`). This project works out the whole 19-cell subdivision:

- closed rational area formulas in `t` for all five cell types and a proof they tile the triangle;
- an incremental *genealogy* of the subdivision, with its splitting grammar;
- a rigidity proof (no three cevians concurrent for `t > 1`) via concurrence determinants;
- the arithmetic of the involution `t -> 1/t`: signed areas, orientation, and the `u = t + 1/t`
  double cover with its invariant / `sqrt(u^2 - 4)` splitting.

The write-up is in `latex/` (elegant `amsart` + Palatino template).

## Layout

```
.
├── latex/          LaTeX source and compiled PDF
│   ├── marion.tex
│   └── marion.pdf
├── code/           Python: figure generation, core library, verification scripts
│   ├── common.py                     shared library (arrangement, classify, palette, paths)
│   ├── extract_representatives.py    -> code/representatives.json
│   ├── fig_family.py                 -> figures/family_panels.png, marion_arrangement.png
│   ├── fig_tree_and_areas.py         -> figures/genealogy_tree.png, cell_areas.png
│   ├── fig_orientation_and_signed.py -> figures/orientation.png, signed_cancellation.png
│   ├── fig_involution.py             -> figures/involution_cevians.png
│   ├── fig_involution_by_type.py     -> figures/involution_by_type.png
│   ├── palette.py                    how the blue->red gradient palette was sampled
│   ├── compare_palettes.py           renders alternative palettes
│   ├── verify_areas.py               explicit area formulas + signed-sum reconciliation
│   ├── verify_concurrence.py         the eight concurrence determinants
│   └── verify_split_rules.py         per-stage split rules of the genealogy
├── figures/        generated PNGs (referenced by the LaTeX via \graphicspath)
├── build.sh
├── requirements.txt
└── README.md
```

## Requirements

- Python 3.10+ with `numpy`, `matplotlib`, `sympy`, `shapely`
  (`pip install -r requirements.txt`).
- A LaTeX distribution with `pdflatex` and the packages `amsart, mathpazo, microtype, amsmath,
  amssymb, mathtools, graphicx, float, booktabs, enumitem, tikz, hyperref` (TeX Live / MiKTeX).

## Build

```bash
pip install -r requirements.txt
./build.sh
```

`build.sh` runs the figure scripts in dependency order (each writes into `figures/`), then calls
`pdflatex` three times in `latex/` (for the table of contents and cross-references). The scripts
resolve their own paths, so they can be run from anywhere.

## Palette

Cell colours run on a light blue->red gradient (outer to inner): corner triangle `#5B8FD6`,
inner triangle `#48B0A5`, quadrilateral `#74C078`, pentagon `#E58C6A`, hexagon `#D9534F`.
Edit `COL` in `code/common.py` and rerun to change them.

## Notes

- All area *ratios* are affine invariants, so the triangle's shape is irrelevant; the code works
  in a fixed reference triangle and displays in an equilateral frame.
- The closed forms are *signed* areas; under `t -> 1/t` they are not individually invariant
  (corner triangles grow, quadrilaterals/pentagons flip sign) but the signed total stays 1.
