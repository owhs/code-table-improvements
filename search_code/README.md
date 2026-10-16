# search_code

Scripts used for the searches, copied unchanged. They expect the original folder layout (`tools/`, `engine/problems/`, `h4/`) and the
parsed Brouwer tables (`data/brouwer/*.json`, not included).

- `tools/qary_search.py`: one search run for one cell. `qary_hunt.py`, `qary_hunt_focus.py`, `qary_hunt_burst.py`: drivers over all cells.
- `engine_problems/qary.py`, `h4/clique.py`, `h4/cw_ops.py`: orbit enumeration, conflict graph, clique search.
- `h4/sweep5.py`: constant-weight sweep over block-cyclic groups.
- `tools/qary_verify.py`, `tools/qary_check3.py`: independent checkers.
- `tools/propagate.py`, `tools/root_map.py`, `tools/cw_roots.py`: derive neighbouring table cells from a new code; find which cells a value comes from.
- `tools/cw_improve.py`, `tools/cw_improve_all.py`: 1-swap and 2-swap local search from a known record code.
