# Code table improvements

New lower bounds for entries in A. E. Brouwer's tables of binary constant-weight, ternary, mixed binary/ternary and 5-ary codes,
a proposed table of spectral-signature codes, and the search scripts. Oscar Sanchez (GitHub: owhs), October 2026, with Claude (Anthropic).

Each result folder has the codes as plain text and a checker that uses only the Python standard library.

| folder | result | table value |
|---|---|---|
| `constant_weight_A31_10_12_1015/` | A(31,10,12) >= 1015 | 961 |
| `constant_weight_A31_10_12/` | A(31,10,12) >= 961 (earlier code, superseded by 1015; on Brouwer's page since 2026-10-14) | 930 before |
| `constant_weight_A32_10_10/` | A(32,10,10) >= 558 | 500 |
| `ternary_mixed/` | 26 codes (24 found by search, 2 derived from them), for example A_3(14,8) >= 123, N(5,8,5) >= 513, N(10,4,4) >= 1104 | 108, 486, 1050 |
| `quinary/A5_10_4_18750/` | A_5(10,4) >= 18750 | 17500 |
| `quinary/A5_9_4_4500/` | A_5(9,4) >= 4500 | 4375 |
| `quinary/reproductions_LO2019/` | independent codes for A_5(8,5) >= 257 and A_5(7,5) >= 57; these values are due to Laaksonen and Ostergard (2019, Table 3) and are not yet on the 5-ary page | 165, 53 |
| `derived_from_LO2019/` | five mixed binary/ternary bounds that follow from Laaksonen and Ostergard's ternary codes, for example N(1,13,4) >= 18711; not yet on the mixed page | 16767 |
| `spectral_signature/` | proposed table: exact sizes of binary constant-weight codes whose words have distinct power spectra, n <= 20 | none |
| `search_code/` | the search and checking scripts | |

The table values were read from Brouwer's pages on 2026-10-16. The mixed codes were compared with Brouwer's page only; the 5-ary values
were also compared with Table 3 of Laaksonen and Ostergard (2019) and with Laaksonen's results page.

Related: [DNA_GC_codes_2026](https://github.com/owhs/DNA_GC_codes_2026), and the Game of Sloanes pull request
https://github.com/gnikylime/GameofSloanes/pull/18.

## Method

A group of coordinate symmetries is fixed, the orbits of words that are valid codes are enumerated, and orbits are chosen by a clique search
(prescribed-automorphism search). Details and references: `METHOD_prescribed_automorphism_search.md`.
