# Third-party code and data

Our [MIT License](LICENSE) covers our own code and the short section adapted
from fido (see below). The bundled copy of MiMIX keeps its own MIT license. The
datasets are not covered by our license; each keeps its original terms. Being
publicly available does not by itself mean a file may be redistributed.

Sources were last checked on October 29, 2026. Some questions are still open;
they are listed in the "Still to check" column below.

## Code

### fido

The function `ssrv.mln` in `simulation/src/01a_main_functions.R` contains a
section adapted from `fido::pibble`. The section ends at the comment
`End of copy from "pibble"`. fido ([jsilve24/fido](https://github.com/jsilve24/fido))
is released under GPL-3.

The author of that code, Justin D. Silverman, gave permission on October 29,
2026 to release this section under our MIT License. This applies only to the
section included here. The fido package itself is still GPL-3. We do not know
which fido version the section was copied from.

### MiMIX

`simulation/src/mimix/MicrobiomeMixedModels.jl/`, `scripts/`, and
`nutnet-analysis/configs/` come from
[nsgrantham/mimix](https://github.com/nsgrantham/mimix), which is MIT
licensed, Copyright (c) 2017 Neal Grantham. A copy of that license is in
[`simulation/src/mimix/LICENSE`](simulation/src/mimix/LICENSE).

We do not know which upstream commit this copy came from, or exactly what was
changed locally. `julia_env/` lists the Julia packages MiMIX needs; those
packages keep their own licenses.

## Data

If you reuse a dataset, cite the original study and follow its license.

| Files | Source and license | Still to check |
|---|---|---|
| Prasse: `2015_prasse_appliedEnvMicrobiology_wetlands/Data/TFW_comdat.RDS`, `TFW_metadata.csv` | Files with the same names are in [QSeq_Model/Data](https://github.com/Djeppschmidt/QSeq_Model/tree/main/Data), from [Epp Schmidt et al. (2023)](https://doi.org/10.1007/s00248-023-02273-z). That repository has no license for the data. (The separate QSeq package is MIT, but that does not cover these files.) | Confirm the files are the same. Get permission to redistribute, or replace the files with download instructions. |
| Stämmler: `2016_stammler_microbiome_micehuman/40168_2016_175_MOESM{1,7,9,11}_ESM.*` | Supplementary files from [Stämmler et al. (2016)](https://doi.org/10.1186/s40168-016-0175-0). The [article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4915049/) is [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) and its data are [CC0](https://creativecommons.org/publicdomain/zero/1.0/) unless stated otherwise. | Check each file for its own notice. Note any changes made to the files. |
| Zemb: `2020_zemb_microOpen_spike/zemb_counts.csv`, `zemb_metadata.csv`, `zemb_qPCR.csv` | [Zemb et al. (2020)](https://doi.org/10.1002/mbo3.977) links to [dataset 10.15454/YWMMH4](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi%3A10.15454%2FYWMMH4), whose count table (`Dada2Table.tab`) is under the Etalab Open License 2.0. `zemb_parse.R` shows how the data were prepared. | Confirm the counts match that dataset. Find where the metadata and qPCR files came from; the dataset license does not cover them. |
| Jin: `2022_Jin_natureComm_technicalReplicates/Supplementary Data {3,4,5,7,8,9}.xlsx` | Supplementary files from [Jin et al. (2022)](https://doi.org/10.1038/s41467-022-28426-1), published under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) except where other credit is given. | Check each file for its own notice. The publisher lists Supplementary Data 1–8 and Source Data, so confirm what "Supplementary Data 9" is (see below). |
| Vandeputte: `Vandeputte/data/41586_2017_BFnature24460_MOESM10_ESM.csv`, `OTU_nochim.csv`, `taxa_assignments.csv` | [Vandeputte et al. (2017)](https://doi.org/10.1038/nature24460). MOESM10 is the source data for Figure 4. We found no open license on the article page. | Find the terms for the source data and where the OTU and taxonomy files came from. If they cannot be redistributed, replace them with download instructions or get permission. |
| SELEX: counts from `data("selex")` in ALDEx2; `selex/selex-truth.txt` | The [ALDEx2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html) credits the data to [McMurrough et al. (2014)](https://doi.org/10.1073/pnas.1322352111). The counts come from the installed ALDEx2 package and are not bundled here. | Find out who made `selex-truth.txt`, from what measurements, and under what terms. |
| `2016_stammler_microbiome_micehuman/stammler_mrckf_out_filter5.RDS` | An older saved MRCKF fit, used only for a plot that is not in the paper. | Confirm who made it and from which data and code. |

### Jin download links

All Jin files can be downloaded from
`https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-28426-1/MediaObjects/`
plus the file name below.

| Local file | Publisher file name |
|---|---|
| Supplementary Data 3.xlsx | `41467_2022_28426_MOESM5_ESM.xlsx` |
| Supplementary Data 4.xlsx | `41467_2022_28426_MOESM6_ESM.xlsx` |
| Supplementary Data 5.xlsx | `41467_2022_28426_MOESM7_ESM.xlsx` |
| Supplementary Data 7.xlsx | `41467_2022_28426_MOESM9_ESM.xlsx` |
| Supplementary Data 8.xlsx | `41467_2022_28426_MOESM10_ESM.xlsx` |
| Supplementary Data 9.xlsx | Probably `41467_2022_28426_MOESM12_ESM.xlsx` (Source Data); not confirmed |

The R and Julia packages you install to run the code are not part of this
repository and keep their own licenses.
