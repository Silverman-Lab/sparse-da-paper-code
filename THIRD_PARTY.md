# Third-party material and source attribution

The root [MIT License](LICENSE) covers the authors' original code and the
fido-derived excerpt described below. Bundled MiMIX code retains its upstream
MIT license and copyright notice. Datasets and other third-party material are
not relicensed by the root license. Public availability alone does not establish
permission to redistribute material.

Sources and license evidence were checked on October 29, 2026. Matching source
filenames do not establish identical file contents. Outstanding checks below
remain unresolved; this document does not grant permissions on behalf of
third-party rights holders.

## Code

### fido-derived excerpt

`simulation/src/01a_main_functions.R`, in `ssrv.mln`, contains a block adapted
from `fido::pibble`, ending at the comment `End of copy from "pibble"`.
Source: [jsilve24/fido](https://github.com/jsilve24/fido).
The [upstream package](https://jsilve24.github.io/fido/) uses GPL-3.

Justin D. Silverman, the author of the copied code, authorized distribution of
this included excerpt under the repository's MIT License on October 29, 2026.
This permission applies to the fido-derived code included here; it does not
change the license of the separately installed fido package or other fido code.
The local adaptation incorporates the block into the Sparse SSRV sampler.
The exact upstream version was not recorded in the supplied bundle and remains
to be identified.

### MiMIX

The bundled implementation in `simulation/src/mimix/MicrobiomeMixedModels.jl/`
and associated `scripts/` and `nutnet-analysis/configs/` originate from
[nsgrantham/mimix](https://github.com/nsgrantham/mimix).
The [upstream MIT license](https://github.com/nsgrantham/mimix/blob/master/LICENSE)
identifies **Copyright (c) 2017 Neal Grantham**. Its full text is preserved in
[`simulation/src/mimix/LICENSE`](simulation/src/mimix/LICENSE). The bundled
Julia package's `Project.toml` also declares MIT and version 1.0.0.

The upstream commit used for the supplied bundle and a complete comparison of
local changes remain to be established. `julia_env/` is a local dependency
environment; installed dependencies retain their own terms.

## Data and derived inputs

Retain source citations, applicable licenses or waivers, and descriptions of
processing for redistributed datasets. An article's license does not override
separate terms for datasets or credited material.

| Local material | Source and license evidence | Remaining work |
|---|---|---|
| Prasse: `2015_prasse_appliedEnvMicrobiology_wetlands/Data/TFW_comdat.RDS`, `TFW_metadata.csv` | Matching filenames occur in [QSeq_Model/Data](https://github.com/Djeppschmidt/QSeq_Model/tree/main/Data), linked by [Epp Schmidt et al. (2023)](https://doi.org/10.1007/s00248-023-02273-z). The README describes derived community datasets. No explicit reuse license was found in QSeq_Model. | Verify file matches and processing; establish redistribution terms with the provider, or supply upstream download instructions. The separate QSeq package's MIT license does not establish terms for these data. |
| Stämmler: `2016_stammler_microbiome_micehuman/40168_2016_175_MOESM{1,7,9,11}_ESM.*` | [Stämmler et al. (2016)](https://doi.org/10.1186/s40168-016-0175-0). The [article notice](https://pmc.ncbi.nlm.nih.gov/articles/PMC4915049/) specifies [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) for the article and [CC0](https://creativecommons.org/publicdomain/zero/1.0/) for its data unless otherwise stated. | Check individual supplement notices and identify unchanged files and local transformations. Cite the original study. |
| Zemb: `2020_zemb_microOpen_spike/zemb_counts.csv`, `zemb_metadata.csv`, `zemb_qPCR.csv` | [Zemb et al. (2020)](https://doi.org/10.1002/mbo3.977) links [dataset 10.15454/YWMMH4](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi%3A10.15454%2FYWMMH4). Dataset V1 lists `Dada2Table.tab` under Etalab Open License 2.0, separately from the article's license. | Match counts to that dataset and document conversion. Trace metadata and qPCR CSVs separately; the OTU dataset license does not establish their terms. Local analysis preprocessing is described in `zemb_parse.R`. |
| Jin: `2022_Jin_natureComm_technicalReplicates/Supplementary Data {3,4,5,7,8,9}.xlsx` | [Jin et al. (2022)](https://doi.org/10.1038/s41467-022-28426-1). Publisher specifies [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) with exceptions for separately credited material, and links supplementary/source data. | Preserve attribution, license link, and change descriptions. Check individual notices and file matches. Publisher lists Supplementary Data 1–8 plus Source Data; local Data 9 needs confirmation. |
| Vandeputte: `Vandeputte/data/41586_2017_BFnature24460_MOESM10_ESM.csv`, `OTU_nochim.csv`, `taxa_assignments.csv` | [Vandeputte et al. (2017)](https://doi.org/10.1038/nature24460). Publisher identifies MOESM10 as Source data to Figure 4. No explicit open redistribution license was found on the article page. | Establish source-data terms and trace the processed OTU/taxonomy files. Use upstream downloads or obtain permission where redistribution terms cannot be established. |
| SELEX: counts loaded with `data("selex")` from ALDEx2; bundled `selex/selex-truth.txt` | [ALDEx2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html) attributes the experiment to [McMurrough et al. (2014)](https://doi.org/10.1073/pnas.1322352111). Counts are loaded from the separately installed package. | Identify who assembled the truth file, its source validation measurements, transformations, and terms. Documentation for the counts does not establish provenance of the truth file. |
| Legacy fit: `2016_stammler_microbiome_micehuman/stammler_mrckf_out_filter5.RDS` | Supplied derived analysis input. Original instructions say it is unnecessary for manuscript figures. | Confirm who generated it and from which data and code. |

### Jin publisher download mapping

Publisher download base:
`https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-28426-1/MediaObjects/`.

| Local workbook | Publisher filename |
|---|---|
| Supplementary Data 3.xlsx | `41467_2022_28426_MOESM5_ESM.xlsx` |
| Supplementary Data 4.xlsx | `41467_2022_28426_MOESM6_ESM.xlsx` |
| Supplementary Data 5.xlsx | `41467_2022_28426_MOESM7_ESM.xlsx` |
| Supplementary Data 7.xlsx | `41467_2022_28426_MOESM9_ESM.xlsx` |
| Supplementary Data 8.xlsx | `41467_2022_28426_MOESM10_ESM.xlsx` |
| Supplementary Data 9.xlsx | Candidate: `41467_2022_28426_MOESM12_ESM.xlsx` (Source Data); unverified |

Analysis scripts and their outputs are mapped in the root [README](README.md). External R and Julia dependencies are installed
separately and retain their own licenses. Where terms remain unresolved,
attribution alone does not establish redistribution permission.
