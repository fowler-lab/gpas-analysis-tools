# gpas-build-tables
CLI to build CRyPTIC data tables from the outputs of GPAS

Philip Fowler, 11 March 2026


Notes

1. `build-tables` should not modify the data fields and be able to cope with either JSON or CSV files downloaded from Pathogena. The glob should be recursive to deal with a sharded file system
2. `shard-files` takes a delimiter (for CRyPTIC this is `.`) and moves the output files into a sharded file system

Issues

1. Why does `variants.parquet` fail to run i.e. is `Killed`? Is it a column type?
2. Why does `genomes` have more rows than it should? -> because of a carriage return in `pipeline_build`.
3. There are multiple ENA run accessions for some UNIQUEIDs; how do I know which was used?
4. Why does mykrobe report a lineage but then record zero median depth for some samples?
5. Why are some samples missing a `main_report`? Example is `site.10.subj.YA00040368.lab.YA00040368.iso.1` / `98bc5c23-d219-43bb-9aab-e8df1c6a0f7e` which I can download via mapping but isn't in the folder, suggesting it failed. Curiously it is the last file in the mapping csv.
