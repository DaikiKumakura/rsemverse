Output of RSEM 1.3.3 itself, from a toy data set (a 6 kb random genome, 4 genes, 2150 simulated 50 nt reads written as
transcript alignments), run with `--calc-pme --calc-ci --append-names --seed 1`. The two samples are the same reads.
Made on November 23, 2026 with `rsem-prepare-reference --gtf` and `rsem-calculate-expression --alignments`.
Gene GC has a transcript of 60 nt, so its effective length is 0 although its expected count is 50.
The .model files are cut after the sections that rsem_model_qc() reads (the sequencing error profiles that follow are omitted). sampC was run with --estimate-rspd.
