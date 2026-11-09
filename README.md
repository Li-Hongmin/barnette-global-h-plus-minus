# Exponentially Many Barnette Graphs with the Global $H^{+-}$ Property

[![verify](https://github.com/Li-Hongmin/barnette-global-h-plus-minus/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/barnette-global-h-plus-minus/actions/workflows/verify.yml) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054750-blue)](https://doi.org/10.5281/zenodo.24054750)

**Hongmin Li** (Institute of Science Tokyo; The University of Tokyo) · ORCID [0000-0003-0228-0600](https://orcid.org/0000-0003-0228-0600)

Main result of the manuscript: an explicit construction gives exponentially many Barnette graphs with the global H^{+-} property (for every ordered pair of distinct edges some Hamilton cycle contains the first and avoids the second).

- Manuscript on rStation: <https://rstation.ai/problems/math/papers/1/v5> (version 5)
- Archived on Zenodo: concept DOI [10.5281/zenodo.24054750](https://doi.org/10.5281/zenodo.24054750) (all versions); this repository mirrors version 5, DOI [10.5281/zenodo.24054768](https://doi.org/10.5281/zenodo.24054768)

## Contents

- `paper/paper.md`, `paper/barnette-global-h-plus-minus.pdf` — the manuscript (version 5).
- `verification/` — the supplement of that version, unpacked unchanged; its own README describes each file.
- `.github/workflows/verify.yml` — reruns every check below on each push.

## What the programs check

| Check | Command (run inside `verification/`) |
|---|---|
| Table 1, minimal covers, cube witnesses, cap isomorphism | `python check_local.py` |
| 40-state transfer closure and separations | `python verify_reflection_40_state_h_plus_minus.py` |
| Final layer of Theorem 8.2 | `python check_final.py` |
| Seed conditions of Theorem 8.2 | `python check_seed.py` |
| Family check to depth 6 (standard) / 5 (reflected) | `python check_family.py` |

The programs check the finite and exact computations used in the manuscript; the mathematical arguments are those written in the paper.

## Reproducing

```sh
cd verification
pip install networkx==3.4.2
python check_local.py
python verify_reflection_40_state_h_plus_minus.py
python check_final.py
python check_seed.py
python check_family.py
```

Requirements: Python 3.9+. Run the commands from a copy if you want to keep the recorded outputs untouched.

## Cite

```bibtex
@misc{li2026barnetteglobalhplusminus,
  author    = {Li, Hongmin},
  title     = {Exponentially Many Barnette Graphs with the Global $H^{+-}$ Property},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24054750},
  url       = {https://doi.org/10.5281/zenodo.24054750}
}
```

## License

Manuscript and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`).
