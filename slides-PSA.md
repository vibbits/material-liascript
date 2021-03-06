<!--
author:   Alexander Botzki
email:    Alexander.Botzki@vib.be
version:  0.1.1
language: en
narrator: US English Female

comment:  slides Protein Structure Analysis

logo: img/Logo.png

link:     https://cdnjs.cloudflare.com/ajax/libs/animate.css/3.7.2/animate.min.css
link:     https://raw.githubusercontent.com/vibbits/material-liascript/master/img/org.css

debug: true

-->

# Protein Structure Analysis

- Sequences, structures and databases
- Experimental methods (X-rays, electrons and NMR)
- Finding and visualising structures from the  Protein Data Bank
- Comparing structures
- Modelling mutations
- Creating homology models

## Sequences and Structures

!?[What is a protein](https://youtu.be/wvTv8TqWC48)<!-- width="890px" height="488px" -->
[^1](https://pdb101.rcsb.org/learn/videos/what-is-a-protein-video)

## Amino acids and peptide structure

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/amino-acids.png)<!-- width="100%" -->

## The Structure-Function Connection

- Folded proteins provide a well-defined 3D arrangement of functional groups, creating microenvironments and active sites.
- Structural changes are often involved in functional  mechanisms  (motor proteins, ...)

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/hemoglobin.png)<!-- height="600px" -->

## Databases

                           {{0-1}}
*********************

[**Uniprot**](https://www.uniprot.org/) approx. 1100,000 sequences

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/uniprot-logo.png)<!-- width="50%" -->

- mainly determined by large-scale DNA sequencing of individual genes or whole genomes

- increasingly automated annotation

*******************

                            {{1}}
*******************

[**Protein Data Bank**](http://www.wwpdb.org) approx. 141,000 experimentally determined structures

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/pdb-logo.png)<!-- width="50%" -->

- Mainly determined by X-ray crystallography and high-resolution NMR spectroscopy
- Automation is increasing, but there are still significant limitations in the rate of solving new structures

*******************

## X-Ray Crystallography - Technical setup


                     {{0-1}}
*********************
![](https://material.bits.vib.be/topics/protein-structure-analysis/images/xray-tech-setup.png)<!-- width="100%" -->

*******************

                       --{{0}}--
In a crystal, a large number  of macromolecules are packed together in a regular grid, with consistent  orientations and relative distances.
When exposed  to an X-ray beam, this  arrangement gives rise to diffracted rays in specific directions, resulting in discrete spots on the planar detector. By rotating the crystal, a series  of images is obtained. From these, the intensities of all the diffracted rays of the crystal can be derived.

                      {{1}}
*******************


![](https://material.bits.vib.be/topics/protein-structure-analysis/images/diffraction-pattern.png)<!-- width="50%" -->

*******************

## X-Ray Crystallography - from diffraction pattern to atomic model

![](img/X-ray-crystallography-1.png)<!-- width="100%" -->

|Diffraction Spot Intensities and Phases| | Electron density|
|:---------------------------------------:|:--------:|:----------------:|
| $F_{obs}(h,k,l)$ and $\phi_{obs}(h,k,l)$ | $$R_{cryst} = \frac{\sum_{h,k,l}F_{obs}-F_{calc}}{\sum_{h,k,l}F_{obs}}$$ |  $\rho(x,y,z)$ |

## The Protein Databank

                                  {{0-1}}
*********************
[http://www.wwpdb.org](http://www.wwpdb.org)

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/wwpdb-welcome-page.png)

*******************

 {{1}}
*******************
- contains structures of  proteins, nucleic acids  and complexes,  determined by X-ray  crystallography, NMR  spectroscopy
- No purely theoretical  or ab initio models  (since 2006)
- Also stores supporting  experimental data
- Full deposition now  required by all peer-reviewed journals

*******************

## Exercise 1: Search the PDB

[Link to exercise 1](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/explore-pdb/tutorial.html)

- Use the UniProt site to search for “dnak”.  
- Use the PDB site to search for “dnak”.
- Compare the UniProt and PDB result lists.
- Use the sequence search function to see if there are structures with sequences similar to that of the  DnaK C-terminal domain.
- Look at the summary pages of a number of  structures and note some interesting properties.
- Select a number of structures and create a report page.

## PDB File Format

                          {{0-1}}
******************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/pdb-file-format.png)<!-- height="600px" -->

******************

                 --{{0}}--
The Protein Data Bank (PDB) format provides a standard representation for macromolecular structure data derived from X-ray diffraction and NMR studies. This representation was created in the 1970's and a large amount of software using it has been written.


                 {{1-2}}
************************************************

Atomic descriptions in classical PDB format

```
ATOM     21  N  ALEU A 392      64.898 118.959  40.827  0.60 24.19           N  
ATOM     22  N  BLEU A 392      64.981 118.958  40.769  0.45 24.73           N  
ATOM     23  CA ALEU A 392      64.818 117.251  41.235  0.60 22.53           C  
ATOM     24  CA BLEU A 392      64.908 117.527  41.044  0.45 23.55           C  
ATOM     25  C  ALEU A 392      63.999 117.260  40.013  0.60 21.73           C  
ATOM     26  C  BLEU A 392      64.972 116.794  39.706  0.45 22.95           C  
ATOM     27  O  ALEU A 392      64.789 117.304  38.508  0.60 19.91           O  
ATOM     28  O  BLEU A 392      64.270 117.164  38.760  0.45 21.80           O  
ATOM     29  CB ALEU A 392      62.994 117.229  41.239  0.60 22.01           C  
ATOM     30  CB BLEU A 392      63.588 117.170  41.739  0.45 23.26           C  
ATOM     31  CG ALEU A 392      62.020 116.866  42.907  0.60 22.43           C  
ATOM     32  CG BLEU A 392      63.411 117.414  43.240  0.45 24.38           C  
ATOM     33  CD1ALEU A 392      63.906 117.128  42.987  0.60 22.36           C  
ATOM     34  CD1BLEU A 392      62.004 117.019  43.655  0.45 24.55           C  
ATOM     35  CD2ALEU A 392      62.660 117.298  44.405  0.60 22.73           C  
ATOM     36  CD2BLEU A 392      64.432 116.611  44.023  0.45 24.44           C  
ATOM     37  N  AASP A 393      65.422 116.659  39.809  0.60 21.96           N  
ATOM     38  N  BASP A 393      65.827 115.776  39.625  0.45 22.93           N  
ATOM     39  CA AASP A 393      65.972 115.093  38.433  0.60 22.44           C  
ATOM     40  CA BASP A 393      65.967 114.977  38.409  0.45 23.84           C  
ATOM     41  C  AASP A 393      64.830 113.969  38.617  0.60 21.62           C  
ATOM     42  C  BASP A 393      64.900 113.889  38.543  0.45 22.45           C  
ATOM     43  O  AASP A 393      64.810 113.091  39.785  0.60 21.11           O  
ATOM     44  O  BASP A 393      64.962 113.054  39.456  0.45 21.89           O  
ATOM     45  CB AASP A 393      66.589 114.675  39.315  0.60 24.09           C  
ATOM     46  CB BASP A 393      67.374 114.360  38.335  0.45 26.75           C  
ATOM     47  CG AASP A 393      66.965 114.032  38.106  0.60 25.18           C  
ATOM     48  CG BASP A 393      67.888 114.196  36.899  0.45 29.40           C  
ATOM     49  OD1AASP A 393      66.303 115.375  36.780  0.60 25.41           O  
ATOM     50  OD1BASP A 393      67.190 114.588  35.938  0.45 31.28           O  
ATOM     51  OD2AASP A 393      68.508 115.272  37.824  0.60 25.36           O  
ATOM     52  OD2BASP A 393      69.014 113.671  36.724  0.45 32.22           O  
```



************************************************

               --{{1}}--
The original short term goal of the working group was to fulfill the mandate set by the IUCr: to define mmCIF data names that needed to be included in the CIF dictionary in order to adequately describe the macromolecular crystallographic experiment and its results. In January 1997, the mmCIF dictionary was completed and submitted to COMCIFS for review, and version 1.0 was released in June 1997.

                 {{2-3}}
************************************************

Atomic descriptions in mmCIF format [^1](http://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html)

```
ATOM   21   N N   A LEU A 1 4   ? 64.898 118.959 40.827 0.60 24.19 ? 392  LEU A N   1
ATOM   22   N N   B LEU A 1 4   ? 64.981 118.958 40.769 0.45 24.73 ? 392  LEU A N   1
ATOM   23   C CA  A LEU A 1 4   ? 64.818 117.251 41.235 0.60 22.53 ? 392  LEU A CA  1
ATOM   24   C CA  B LEU A 1 4   ? 64.908 117.527 41.044 0.45 23.55 ? 392  LEU A CA  1
ATOM   25   C C   A LEU A 1 4   ? 63.999 117.260 40.013 0.60 21.73 ? 392  LEU A C   1
ATOM   26   C C   B LEU A 1 4   ? 64.972 116.794 39.706 0.45 22.95 ? 392  LEU A C   1
ATOM   27   O O   A LEU A 1 4   ? 64.789 117.304 38.508 0.60 19.91 ? 392  LEU A O   1
ATOM   28   O O   B LEU A 1 4   ? 64.270 117.164 38.760 0.45 21.80 ? 392  LEU A O   1
ATOM   29   C CB  A LEU A 1 4   ? 62.994 117.229 41.239 0.60 22.01 ? 392  LEU A CB  1
ATOM   30   C CB  B LEU A 1 4   ? 63.588 117.170 41.739 0.45 23.26 ? 392  LEU A CB  1
ATOM   31   C CG  A LEU A 1 4   ? 62.020 116.866 42.907 0.60 22.43 ? 392  LEU A CG  1
ATOM   32   C CG  B LEU A 1 4   ? 63.411 117.414 43.240 0.45 24.38 ? 392  LEU A CG  1
ATOM   33   C CD1 A LEU A 1 4   ? 63.906 117.128 42.987 0.60 22.36 ? 392  LEU A CD1 1
ATOM   34   C CD1 B LEU A 1 4   ? 62.004 117.019 43.655 0.45 24.55 ? 392  LEU A CD1 1
ATOM   35   C CD2 A LEU A 1 4   ? 62.660 117.298 44.405 0.60 22.73 ? 392  LEU A CD2 1
ATOM   36   C CD2 B LEU A 1 4   ? 64.432 116.611 44.023 0.45 24.44 ? 392  LEU A CD2 1
ATOM   37   N N   A ASP A 1 5   ? 65.422 116.659 39.809 0.60 21.96 ? 393  ASP A N   1
ATOM   38   N N   B ASP A 1 5   ? 65.827 115.776 39.625 0.45 22.93 ? 393  ASP A N   1
ATOM   39   C CA  A ASP A 1 5   ? 65.972 115.093 38.433 0.60 22.44 ? 393  ASP A CA  1
ATOM   40   C CA  B ASP A 1 5   ? 65.967 114.977 38.409 0.45 23.84 ? 393  ASP A CA  1
ATOM   41   C C   A ASP A 1 5   ? 64.830 113.969 38.617 0.60 21.62 ? 393  ASP A C   1
ATOM   42   C C   B ASP A 1 5   ? 64.900 113.889 38.543 0.45 22.45 ? 393  ASP A C   1
ATOM   43   O O   A ASP A 1 5   ? 64.810 113.091 39.785 0.60 21.11 ? 393  ASP A O   1
ATOM   44   O O   B ASP A 1 5   ? 64.962 113.054 39.456 0.45 21.89 ? 393  ASP A O   1
ATOM   45   C CB  A ASP A 1 5   ? 66.589 114.675 39.315 0.60 24.09 ? 393  ASP A CB  1
ATOM   46   C CB  B ASP A 1 5   ? 67.374 114.360 38.335 0.45 26.75 ? 393  ASP A CB  1
ATOM   47   C CG  A ASP A 1 5   ? 66.965 114.032 38.106 0.60 25.18 ? 393  ASP A CG  1
ATOM   48   C CG  B ASP A 1 5   ? 67.888 114.196 36.899 0.45 29.40 ? 393  ASP A CG  1
ATOM   49   O OD1 A ASP A 1 5   ? 66.303 115.375 36.780 0.60 25.41 ? 393  ASP A OD1 1
ATOM   50   O OD1 B ASP A 1 5   ? 67.190 114.588 35.938 0.45 31.28 ? 393  ASP A OD1 1
ATOM   51   O OD2 A ASP A 1 5   ? 68.508 115.272 37.824 0.60 25.36 ? 393  ASP A OD2 1
```

************************************************

                 {{3}}
************************************************

mmCIF file looks like a paired collection of data item names and values.

```
data_1DKX
#
_entry.id   1DKX
#
_audit_conform.dict_name       mmcif_pdbx.dic
_audit_conform.dict_version    5.279
_audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic
#
loop_
_database_2.database_id
_database_2.database_code
PDB   1DKX         
WWPDB D_1000172826
#
_pdbx_database_status.status_code                     REL
_pdbx_database_status.entry_id                        1DKX
_pdbx_database_status.recvd_initial_deposition_date   1996-06-03
_pdbx_database_status.deposit_site                    ?
_pdbx_database_status.process_site                    ?
_pdbx_database_status.SG_entry                        .
_pdbx_database_status.pdb_format_compatible           Y
_pdbx_database_status.status_code_mr                  ?
_pdbx_database_status.status_code_sf                  ?
_pdbx_database_status.status_code_cs                  ?
#
loop_
_audit_author.name
_audit_author.pdbx_ordinal
'Zhu, X.'           1
'Zhao, X.'          2
'Burkholder, W.F.'  3
'Gragerov, A.'      4
'Ogata, C.M.'       5
'Gottesman, M.E.'   6
'Hendrickson, W.A.' 7
#
_citation.id                        primary
_citation.title                     'Structural analysis of substrate binding by the molecular chaperone DnaK.'
_citation.journal_abbrev            Science
_citation.journal_volume            272
_citation.page_first                1606
_citation.page_last                 1614
_citation.year                      1996
_citation.journal_id_ASTM           SCIEAS
_citation.country                   US
_citation.journal_id_ISSN           0036-8075
_citation.journal_id_CSD            0038
_citation.book_publisher            ?
_citation.pdbx_database_id_PubMed   8658133
_citation.pdbx_database_id_DOI      ?
#
```

************************************************

## Occupancy

                 {{0}}
************************************************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/occupancy.png)<!-- width="80%" -->
[^1](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/dealing-with-coordinates)

************************************************

                 --{{0}}--
Macromolecular crystals are composed of many individual molecules packed into a symmetrical arrangement. In some crystals, there are slight differences between each of these molecules. For instance, a sidechain on the surface may wag back and forth between several conformations, or a substrate may bind in two orientations in an active site, or a metal ion may be bound to only a few of the molecules. When researchers build the atomic model of these portions, they can use the occupancy to estimate the amount of each conformation that is observed in the crystal. For most atoms, the occupancy is given a value of 1, indicating that the atom is found in all of the molecules in the same place in the crystal. However, if a metal ion binds to only half of the molecules in the crystal, the researcher will see a weak image of the ion in the electron density map, and can assign an occupancy of 0.5 in the PDB structure file for this atom. Occupancies are also commonly used to identify side chains or ligands that are observed in multiple conformations. The occupancy value is used to indicate the fraction of molecules that have each of the conformations. Two (or more) atom records are included for each atom, with occupancies like 0.5 and 0.5, or 0.4 and 0.6, or other fractional occupancies that sum to a total of 1.

## Related Web sites

- Nucleic Acid Database: DNA and RNA structures

  [http://ndbserver.rutgers.edu/](http://ndbserver.rutgers.edu/)

- PDB-REDO: automatically re-refined deposited structures, using the latest methods

  [http://www.cmbi.ru.nl/pdb_redo/](http://www.cmbi.ru.nl/pdb_redo)

- EBI: many online tools for structure analysis

  [http://www.ebi.ac.uk/Tools/structure/](http://www.ebi.ac.uk/Tools/structure/)

- Replaced Electron Density Server: convenient overview of quality parameters for crystal structures

  [http://www.ebi.ac.uk/pdbe/litemol](http://www.ebi.ac.uk/pdbe/litemol)

## High-Resolution NMR Spectrometry - Principles


                 {{0}}
************************************************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/nmr-peaks-to-structure.png)<!-- height="307px"-->

************************************************

                              --{{0}}--
Many atomic nuclei, including the ubiquitous hydrogen nuclei,  resonate at specific radio frequencies when placed in a  strong, uniform magnetic field. The chemical environment of each individual atom slightly modulates its exact resonance  frequency. In macromolecules with thousands of atoms, many different  effects combine to generate an extremely complicated pattern of chemical shifts, which therefore more or less uniquely  identify each atom. Multidimensional spectra allow these  frequencies to be assigned to specific atoms.


## High-Resolution NMR Spectroscopy - Nuclear Overhauser Effect

                 {{0}}
************************************************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/nmr-noe.jpg)<!-- height="600px" -->

************************************************

                          --{{0}}--
When two atoms are near each other in 3D space, they can exchange magnetization, giving rise to crosspeaks at the  intersection of their respective frequencies. This nuclear Overhauser effect (NOE) is used to effectively measure the distances between pairs of atoms, at least qualitatively.

## High-Resolution NMR Spectroscopy

                 {{0}}
************************************************

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/nmr-model-example.png)<!-- height="600px" -->

************************************************

                --{{0}}--
After identification of the atoms by means of their unique chemical shifts, distance restraints are derived from the Overhauser crosspeaks. An extended model of the protein is generated, and then condensed into a shape that is consistent with as many of distance restraints as possible.

## Other methods

- Electron microscopy (and especially Cryo-electron microscopy): Electron crystallography and single particle reconstruction
- Small-angle X-ray and neutron scattering (SAXS and SANS)

![](https://material.bits.vib.be/topics/protein-structure-analysis/images/saxs.png)<!-- height="600px" -->

## Related Web sites

- BioMagResBank: experimental data for NMR-  derived structures (lists of distance restraints  and other experimentally derived properties)

  [http://www.bmrb.wisc.edu/](http://www.bmrb.wisc.edu/)

- BioIsis: database of SAXS-derived structures

  [http://www.bioisis.net](http://www.bioisis.net)

- EMBL database of SAXS-derived structures

  [http://www.sasbdb.org](http://www.sasbdb.org)

- EM Databank for cryo-EM structures

  [http://www.emdatabank.org](http://www.emdatabank.org.left)

## Assessing Structure Quality

General geometric properties (bond lengths and angles, Ramachandran distribution, …):
[MolProbity](http://molprobity.biochem.duke.edu/)

**Crystal Structures**

- Diffraction data resolution and completeness (PDB)
- Final $ R_{cryst} $ and $ R_{free} $ factors (PDB)
- Local correspondence to electron density [EDS](https://www.ebi.ac.uk/pdbe/eds)

**NMR Structures**

- Number of distance restraints and additional experimental data sources [BMRB](http://www.bmrb.wisc.edu/)
- Local restraint density [on-line NMR constraint analyser](http://molsim.sci.univr.it/bioinfo/tools/constraint/index.html)

**Other techniques**

- Difficult to generalise: carefully read and consider the  published protocol

## Molecular Graphics Software

- [PyMOL](http://www.pymol.org/): high-quality output,  good examples on [Wiki](http://www.pymolwiki.org/)
- [Chimera](http://www.cgl.ucsf.edu/chimera/): good  documentation on website
- [VMD](http://www.ks.uiuc.edu/Research/vmd/):  excellent for the analysis of MD trajectories
- [Yasara](http://www.yasara.org)  
- [SwissPDBViewer](http://spdbv.vital-it.ch/)

## YASARA

- Yasara View is freely available and provides  basic visualisation functions.

- Yasara Model, Dynamics and Structure provide additional functionality.

- An add-on module for NMR structures is available.

- The program can be combined with the WHAT-IF program for structure validation, and with FoldX for energy calculations.

## Exercise 2: Show a structure

[Link to exercise 2](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/visualise-structures/tutorial.html)

- Load PDB entry 1TRZ using the File menu.

- Compare the default representations (F1-F8)  and use the various controls to change the view  of the protein.

- Explore the Edit and View menus to change  various aspects of the graphical representation  of the structure.

- Examine the hydrogen bonding patterns.  Display a molecular surface.

- Create an interesting view of the molecule and  save it as a high-resolution image.

## Protein folds are the structures of domains

- Similarities in assembly of secondary structure elements

- So not based on sequence like motifs but on 3D structure

- Folds represent the shapes of protein domains

* SCOP (http://scop.mrc-lmb.cam.ac.uk/scop/)
  * Class, Fold, Superfamily, Family
  The human hemoglobin protein (1HHO-A), is a domain classfied by SCOP as the following.
  * Class: all helical proteins; fold: globin-like; superfamily: globin-like; and family: globins.

* CATH (http://www.cathdb.info/)
  * class (equivalent to SCOP class), architecture, topology (equivalent to SCOP fold), and homologous superfamily (equivalent to SCOP superfamily).
  1HHO-A is a domain classied by CATH as the following.
  * Class: mainly helical; architecture: orthogonal bundle; topology: globin-like; and superfamily: globins.

## Predict protein structures by fold recognition

1. Search SCOP/CATH for protein with same fold and known 3D structure
2. Align each amino acid of query sequence to a position in the template structure
3. Evaluate how well the sequence fits the fold and select best-fit fold
4. Build structural model of query based on alignment with selected fold

- Phyre (http://www.sbg.bio.ic.ac.uk/phyre2/html/page.cgi?id=index)
- HHpred (http://toolkit.lmb.uni-muenchen.de/hhpred)
- DescFold (http://202.112.170.199/DescFold/)

Works because:

- Number of different folds in nature is fairly small (approximately 1300)
- 90% of new submissions in PDB have similar folds to those already in PDB
- Not always accurate [^1](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.25823)

[^2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3003448/#!po=87.5000)

## Guidelines to improve fold recognition results

- Run as many methods as you can
- Run each method on many sequences from your homologous protein family
- After all of these runs, build up a consensus picture of the likely fold
- Compare function of your protein to function of the proteins with the likely fold
- Compare secondary structure of your protein to that of the likely fold

## Similarity searches based on 3D structure

Similarity on structural level: aligning 3D structures

Structure of query protein is known and aligned to PDB structures

- DALI(http://ekhidna.biocenter.helsinki.fi/dali_server/)
- VAST+(https://www.ncbi.nlm.nih.gov/Structure/vastplus/vastplus.cgi)

Compare proteins with low sequence similarity: similar structure implies homology -> same function

Can help to find active sites

[Additional information](https://web.stanford.edu/class/cs279/lectures/lecture5.pdf)

## Structure Alignment

Structures can be structurally aligned to achieve the best possible match between corresponding atoms. It is  possible to consider all atoms, or only a subset  (such as the Cα atoms).

When the structures are very similar, determining  which atoms should be matched to each other is trivial. When there are larger differences, it takes more preparation.

Different algorithms use different combinations of sequence- and secondary of tertiary structure-based information.

- [DALI](http://ekhidna.biocenter.helsinki.fi/dali_server/) [^1](Holm and Sander, 1993 and Holm 2019 https://doi.org/10.1093/bioinformatics/btz536)
- [TM-Align](https://zhanglab.ccmb.med.umich.edu/TM-align/)
- [DeepAlign](https://www.nature.com/articles/srep01448)
- [LGA](https://dx.doi.org/10.1093%2Fnar%2Fgkg571)
- [CE](https://doi.org/10.1093/protein/11.9.739)
- [FatCat](https://doi.org/10.1093/bioinformatics/btg1086)
- [MUSTANG](http://lcb.infotech.monash.edu.au/mustang/mustang_psfb-final.pdf)
- [MMLigner](https://doi.org/10.1093/bioinformatics/btw757)

## Exercise 3: Compare Structures

[Link to exercise 3](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/compare-structures/tutorial.html)

- Download the five provided PDB files and open  them in Yasara.

- Use the `Analyze|Align|Objects` with MUSTANG function to align the four last objects with the first one.

- Use the space bar to open the text console and  see the reported root mean square deviations as well as the number of residues matched.

$$ rmsd = \sqrt{\frac{1}{N}\sum_{i=1}^{N}R_{i}^{2}} $$

- Color all structures by B factor and compare the  distribution to the local variability of the structures.

## PDB File Surprises

- Missing atoms, residues and loops  
- Multiple molecules in the asymmetric unit
- Incomplete oligomers due to coincident crystal and  oligomer symmetry
- Cα-only models
- alternate conformations
- multiple models, especially for NMR ensembles
- use of B factors to represent other properties  
- other non-standard extensions (PDBQT, ...)

## Force Fields

- Energy terms representing physical interactions

  - Covalent bond lengths and angles
  - Dihedral angles and van der Waals forces (steric effects)
  - Electrostatic interactions and hydrogen bonds
  - …

- Energy can be minimized, and forces can literally be derived from the potential function.

- Consistency and careful consideration of the properties to be simulated are essential.

## Force Field Terms

Each energy term has a functional form, which includes one or more parameters:

- Covalent bond energy term

  $ E_{bond} = k(r-r_{0})^{2} $ with $r_0$ being the equilibrium length

- Van der Waals contact energy term

  $ E_{VdW} = C_6r^{-6} + C_{12}r^{-12} $ with $C_6$ attractive, $C_{12}$ repulsive term

The parameters are collectively optimized to  reproduce a chosen set of experimentally observed parameters.

A given force field should be used as a  consistent system, and can only be used to predict properties that are covered by the  training set.

## FoldX

Is designed for quantitative modelling of the contribution of structural interactions to the stability of proteins and protein complexes. It also supports  protein/DNA complexes.

The force field describes the different interactions in a protein structure or complex in analytical terms. It has  been calibrated using a set of experimentally determined  stabilities.

Applications include the optimisation of structures, the calculation of the stability of complexes, and predicting  the effect of amino-acid or base-pair mutations on these  properties.

![FoldX Force Field](img/foldx-formula.PNG)

[FoldX](https://dx.doi.org/10.1093%2Fnar%2Fgki387)

## The FoldX Plugin for YASARA

In order to make FoldX more accessible and  integrate its functions into Yasara, dr. Joost van  Durme (SWITCH laboratory) made a Yasara plugin  module that can apply FoldX functions to structures  that are loaded as Yasara objects.

This greatly simplifies the use of FoldX, and allows  for a quick visual analysis of the resulting changes  in the structures.

More information can be found at [wiki](http://foldxyasara.switchlab.org/index.php/) [FoldX](http://foldxsuite.crg.eu/)

## Exercise 4a: Repair a PDB File

[Link to exercise 4](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/mutate-structure/tutorial.html)

- Load the 1CRN PDB file.
- Use the “Repair object” option in the `Analysis|FoldX` menu to activate the corresponding FoldX function.
- Select the object to repair.

This exports the object as a temporary PDB file,  starts FoldX with the appropriate options, and loads the repaired PDB file as a new object in Yasara.

- Compare the original and repaired objects.  
- Describe the changes that were introduced.

## Exercise 4b: Model a Mutation

[Link to exercise 4](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/mutate-structure/tutorial.html)

- Load the 2AC0.sce Yasara scene file.
- Set an appropriate structure representation.
- Locate residue Ala159 using the sequence view, and right-click to access the `FoldX|Mutate` residue function. Change it to a Trp residue.
- Look at the effect of the substitution on the structure, and use the space bar to open the text  console and read the output of the FoldX calculation.
- Mutate Arg273 to an alanine side chain. Discuss  the effects of this substitution.

## Homology Modelling

- When a structure is available for a protein with a  similar sequences, it is possible to predict the  structure of a new sequence with varying degrees of  confidence.
- Use PSI-BLAST/DELTA-BLAST to detect sequences with similar structures.
- All homology modelling procedures start from an  alignment of the template and target sequences.  The quality of this alignment will have a major  impact on the resulting model.
- Available applications include stand-alone programs  (Modeller, FoldX, …) and web-based services (such as SwissModel).

## Exercise 5: Make a Homology Model using Swiss Model

[Link to exercise 5](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/homology-modeling/tutorial.html)


## Exercise 6: Study Protein-Ligand Interactions

[Link to exercise 6](https://material.bits.vib.be/topics/protein-structure-analysis/tutorials/protein-ligand-interaction/tutorial.html)

## On-Line References

- [Crystallography 101 (Bernhard Rupp)](http://www.ruppweb.org/Xray/101index.html)

- [Protein NMR, A Practical Guide (Vicky Higman)](http://www.protein-nmr.org.uk/)

- [Model validation course](http://xray.bmc.uu.se/gerard/embo2001/modval/index.html)

- [Assessing model quality](http://spdbv.vital-it.ch/TheMolecularLevel/ModQual/)

- [Lectures by Burkhard Rost on protein structure prediction](https://www.youtube.com/channel/UCU6j8BG4RbEtTgyIZJ6Vpow)
