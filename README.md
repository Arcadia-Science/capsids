# Searching for Capsids in Venom Proteomes

This README documents searching for viral capsid proteins sets of venom transcriptomes. To use the scripts and programs described in the README, create an environment with:
```
conda env create -n capsid_searches -f environment.yml
```

## Searching for viral capsid proteins
Run capsid HMMs against sets of venomous species proteins with `python3 scripts/run_parse_hmms.py dbs/capsid_hmm.hmm proteins/all_proteins/all_venom_proteins.fasta --evalue 1e0` with no evalue cutoff so that I could just hand all of the results to the strike team. The results are in `results/2023-04-11_capsid_v1/` where I experimented with some different evalue cutoffs before settling with no cutoff based from input from strike/pilot people.

To combine multiple FASTA files into one FASTA file where the filename is carried through to the header (important if you want to retain the accession or species number for example), use the script `reorganize_fastas.py`.

## Searching for a specific viral capsid protein from Bracovirus
We were specifically interested in a Bracovirus capsid protein that wasn't being captured from the previous HMM search, and found that for whatever reason this particular capsid wasn't included in the HMM above. We found the protein sequence for this capsid, BLASTed it in NCBI, and downloaded the hits and turned it into an HMM.

Then we created an alignment of the multi-FASTA file of BLAST results to the capsid protein with `muscle -align vp39_psiBLAST_hits.fasta -output vp39_psiBLAST_hits.aln` then build an HMM profile with `hmmbuild vp39_capsid.hmm vp39_psiBLAST_hits.aln`. This was used to run `hmmsearch` against the venom proteins using the `run_parse_hmms.py` script with: `python3 ../scripts/run_parse_hmms.py vp39_capsid.hmm ../proteins/all_proteins/all_venom_proteins.fasta --evalue 1e0 --output_file ../results/vp39_capsid_summary.txt`

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
