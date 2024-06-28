Thanks for trying out LAM!

### Input file requirements

LAM requires a Limma output file, from either 450K or EPIC type microarrays.
The CSV file requires that the probe names have their own column heading.
It also needs a column of t values, with the header matching "t".
Check out the [Demo CSV file](https://raw.githubusercontent.com/markziemann/gmea_app/main/example_data/DMPs_sample.csv) which corresponds to a 450K array.

LAM also requires a GMT file to work, which is the standard format for gene set
libraries.
You can find many useful GMT files at the [MSigDB website](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp).
We also provide a [Demo GMT file](https://raw.githubusercontent.com/markziemann/gmea_app/main/example_data/mysets.gmt) for testing.

After uploading, click on the CSV and GMT checks to see if they look okay.

### Options

Prioritisation of enrichment results can be done by "Effect" or "Significance".
Typically "Effect" is better because we want to identify pathways with big
biological changes.
But sometimes, very few pathways are found, which is when we might want to
prioritise by significance.

Minimum gene set size is 5 by default, but some folks like to use 10 instead.

Array chip can be 450K or EPIC chip.

### Running enrichment and getting results

After uploading is complete and parameters are set to your liking, click on
"Enrichment Result" tab which will kick off the analysis.

But be patient, because it may take a few minutes to calculate those enrichment
results, especially if you are dealing with EPIC arrays and have many gene sets.
Once the enrichment table appears, you can then click on the "Enrichment Plot"
tab.

If that goes well, then you can hit the "Download Report" button.
It might also take a few minutes, but it is a worthwhile step because it will
give you a record of all the pathways found.