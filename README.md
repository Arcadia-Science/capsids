# Searching for Capsids in Venom Proteomes

This README documents searching for viral capsid proteins sets of venom transcriptomes. To use the scripts and programs described in the README, create an environment with:
```
conda env create -n capsid_searches -f environment.yml
```

The most updated set of results are in `2023-06-01-updated-vog-results` where VOG capsid HMMs were run against the entire venom/tick database and 
additional tick proteins from the TSA. The specific capsid VOGs from release 217 are listed in `metadata/Updated_capsid_VOGs_v217.txt`. 

## Searching for viral capsid proteins
Run capsid HMMs against sets of venomous species proteins with `python3 scripts/run_parse_hmms.py dbs/capsid_hmm.hmm proteins/all_proteins/all_venom_proteins.fasta --evalue 1e0` with no evalue cutoff. The results are in `results/2023-04-11_capsid_v1/` where I experimented with some different evalue cutoffs before settling with no cutoffs.

To combine multiple FASTA files into one FASTA file where the filename is carried through to the header (important if you want to retain the accession or species number for example), use the script `reorganize_fastas.py`.

## Searching for a specific viral capsid protein from Bracovirus
We were specifically interested in a Bracovirus capsid protein that wasn't being captured from the previous HMM search, and found that for whatever reason this particular capsid wasn't included in the HMM above. We found the protein sequence for this capsid, BLASTed it in NCBI, and downloaded the hits and turned it into an HMM.

Then we created an alignment of the multi-FASTA file of BLAST results to the capsid protein with `muscle -align vp39_psiBLAST_hits.fasta -output vp39_psiBLAST_hits.aln` then build an HMM profile with `hmmbuild vp39_capsid.hmm vp39_psiBLAST_hits.aln`. This was used to run `hmmsearch` against the venom proteins using the `run_parse_hmms.py` script with: `python3 ../scripts/run_parse_hmms.py vp39_capsid.hmm ../proteins/all_proteins/all_venom_proteins.fasta --evalue 1e0 --output_file ../results/vp39_capsid_summary.txt`. The files for creating this vp39 HMM and the HMM are in `hmms`.

## BLAST searches for _Amblyomma americanum_
BLAST should already be installed in the environment you created above. First download the tick genome:
```
wget -O Arcadia_Amblyomma_americanum_asm001_purged_cleanedup1.fasta https://zenodo.org/record/7783368/files/Arcadia_Amblyomma_americanum_asm001_purged_cleanedup1.fasta\?download\=1
```

Then extract protein sequences of interest from the pangenome that was created by the Arcadia-Science/rehgt pipeline:
```
for acc in EEC17452.1 KAG0427517.1 KAG0443635.1 EEC12039.1 KAG0444350.1 KAG0420414.1 XP_042148722.1
do
grep -A 1 ${acc} ../../outputs/genus_pangenome_clustered/Ixodes_cds_rep_seq.fasta >> query_cds.fasta
done
```

Make a BLAST db of the genome and run blast:
```
makeblastdb -in Arcadia_Amblyomma_americanum_asm001_purged_cleanedup1.fasta -dbtype nucl -out Arcadia_Amblyomma_americanum_asm001_purged_cleanedup1

 blastn -task megablast -db Arcadia_Amblyomma_americanum_asm001_purged_cleanedup1 -query query_cds.fasta -dust no -max_target_seqs 100 -outfmt 6 -out outputfile.txt
```

And do the same for proteins:
```
for acc in EEC17452.1 KAG0427517.1 KAG0443635.1 EEC12039.1 KAG0444350.1 KAG0420414.1 XP_042148722.1
do
grep -A 1 ${acc} ../../outputs/genus_pangenome_clustered/Ixodes_aa_rep_seq.fasta >> query_aa.fasta
done

tblastn -db Arcadia_Amblyomma_americanum_asm001_purged_cleanedup1 -query query_aa.fasta -max_target_seqs 100 -outfmt 6 -out outputfile_tblastn.txt
```

The queries and results files are in `results/2023-04-20-blast-results`.

## Structural clustering to identify capsids in *Ornithodoros turicata*

In order to compare structures of tick proteins to viral capsid proteins, we ran [ProteinCartography](https://github.com/Arcadia-Science/ProteinCartography/tree/v0.4.2) using v0.4.2. 

First, we cloned the ProteinCartography repository and followed the [README.md](https://github.com/Arcadia-Science/ProteinCartography/blob/main/README.md) to install the conda environments. We used the notebook in this repository, `notebooks/1_prep_metadata.ipynb` to fetch UniProt metadata for the *Ornithodoros* proteins in `metadata/ornithodoros.txt`, which we created by searching [UniProt](https://www.uniprot.org) for all *Ornithodoros* proteins. This notebook should be ran using the [`envs/web_apis.yml`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.2/envs/web_apis.yml) environment from the ProteinCartography repository.

We then used the notebook in this repository, `notebooks/2_get_alphafold_structures.ipynb` to fetch predicted structures from the [AlphaFold database](https://alphafold.ebi.ac.uk) using the list of proteins in `metadata/ornithodoros.txt`. This notebook should also be ran using the [`envs/web_apis.yml`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.2/envs/web_apis.yml) environment from the ProteinCartography repository. We folded all of the proteins from the [Virus Orthologous Groups Database](https://vogdb.org) (VOG database) that contained "capsid" in the name of the group using [ESMFold](https://www.science.org/doi/10.1126/science.ade2574).

Finally, we ran ProteinCartography using the `config_ff.yml` that can be found in the notebooks folder of this repository. The pipeline was initiated from inside the ProteinCartography repository with the `envs/cartography_tidy` environment activated using the following command: 
```
snakemake --snakefile Snakefile_ff --configfile ../capsids/metadata/config_ff.yml --use-conda --cores 4
```
We then used the notebook in this repository, `notebooks/3_plotting_overlays.ipynb` to create custom plotting overlays for the UMAP. This notebook should be ran using the [`envs/plotting.yml`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.2/envs/plotting.yml) environment from the ProteinCartography repository. 

The results from this analysis can be found in the [Zenodo](10.5281/zenodo.12796464) repo.
