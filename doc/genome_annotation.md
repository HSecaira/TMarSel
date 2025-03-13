# Genome/MAG annotation

Below are the recommendations for formatting the output of KEGG and EggNOG annotation to contain only three columns:

| orf | bit_score | gene_family |
| --- | --- | --- |
| G000006605_1 | 1130.0 | COG0593 |
| G000006605_1 | 644.8 | K02313 |
| ... | ... | ... |

## KEGG annotation with KofamScan. See [documentation](https://github.com/takaram/kofam_scan).

Execute KofamScan

```bash
exec_annotation -f detail-tsv --no-report-unannotated -o $output $query
```

In the output from KofamScan, the `orf`, `bit_score`, and `gene_family` are in columns 2, 5, and 3. Also, the best hit of every gene family is highlighted with an asterisk. The best hit refers to an score higher than the threshold defined for the gene family in KEGG. Therefore, we want to extract only the columns of interest and the gene families with significant thresholds.

```bash
awk -F'\t' '/^\*/ {print $2, $5, $3}' OFS='\t' $output > new_dir/$output
```

## EggNOG annotation with EggNOG-mapper. See [documentation](https://github.com/eggnogdb/eggnog-mapper)

Execute EggNOG-mapper

```bash
emapper.py -i $input -o $output
```

In the output from EggNOG-mapper, the `orf`, `bit_score`, and `gene_family` are in columns 1, 4, and 5. Column 5 contains gene families across taxonomic scales. We will use only the most basal assignment.

```bash
awk -F'\t' '!/^#/ {split($5, a, "|"); out=""; for (i in a) {split(a[i], b, "@"); if (i == 1 || out == "") out = b[1]; else out = out "|" b[1]} print $1, $4, out}' OFS='\t' $output > new_dir/$output
```



