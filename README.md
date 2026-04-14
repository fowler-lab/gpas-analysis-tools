# gpas-build-tables
CLI to build CRyPTIC data tables from the outputs of GPAS

Philip Fowler, 11 March 2026


Notes

1. `build-tables` should not modify the data fields and be able to cope with either JSON or CSV files downloaded from GPAS. The glob should be recursive to deal with a sharded file system
2. `correct-tables`. This is a separate command that takes the output of `build-tables` and joins `mutations` to `variants` to pick up some useful genetics metrics.
3. `shard-files` takes a delimiter (for CRyPTIC this is `.`) and moves the output files into a sharded file system (not yet implemented)

Issues

1. There are multiple ENA run accessions for some UNIQUEIDs; how do I know which was used?
2. Why does mykrobe report a lineage but then record zero median depth for some samples?
3. Why are some samples missing a `main_report`? Example is `site.10.subj.YA00040368.lab.YA00040368.iso.1` / `98bc5c23-d219-43bb-9aab-e8df1c6a0f7e` which I can download via mapping but isn't in the folder, suggesting it failed. Curiously it is the last file in the mapping csv.

## Setting up the GPAS CLI

The below sets up the GPAS CLI using a branch that contains a fix for downloading files using regexs.

```
git clone git@github.com:GlobalPathogenAnalysisService/client.git
cd client
git checkout fix/regex-download
python -m venv env
source env/bin/activate
pip install . "hostile==1.1.0"
```
Now we need to make sure our CLI communicates with the `sp3dev` GPAS instance. Issue the below in the same terminal.


```
export GPAS_UPLOAD_HOST="api.upload-sp3dev.gpas.global"
export GPAS_HOST="portal-sp3dev.gpas.global"
export GPAS_APP_HOST="app-sp3dev.gpas.global"
export AUTH_AUDIENCE="https://portal-sp3dev.gpas.global"
export AUTH_DOMAIN="eit-pathogena-dev.uk.auth0.com"
export AUTH_ISSUER="https://account.portal-sp3dev.gpas.global"
```

You also need to export the `AUTH_CLIENT_ID` but as that is a secret I won't put it here. Ask someone! Now this should let you login and authenticate using your `sp3dev` user/password and also MFA.

```
gpas auth
```

Now you have a token that is valid for 24 hours and can be used to download files from GPAS. You can check the token is valid by running `gpas auth --check-expiry`.


## NTMs test-suite

### Creating the data tables

First let's setup `gpas-analysis-tools`. Can do this in the same `env` environment as above (i.e. stay in the same terminal!).

```
cd ..
git clone https://github.com/fowler-lab/gpas-analysis-tools
cd gpas-analysis-tools
pip install -e .
```

Check it is working

```
gat -h
```

### Downloading the output files

Using the same terminal (as you've now got the `gpas` CLI setup), we can download the output files for the NTM test suite. This is a set of six samples that have already been run through the GPAS pipeline. In `test-suite/ntms/` folder you'll find `TEST_MASTER_TABLE.csv` which is setup to resemble a mapping file.  First let's get the `main_report.json` output files.

```
cd test-suite/ntms/
gpas download --filenames main_report.json --output-dir data/main_report/ --rename TEST_MASTER_TABLE.csv
gpas download --filenames '.*effects.csv' --output-dir data/effects/ --rename TEST_MASTER_TABLE.csv
gpas download --filenames '.*predictions.csv' --output-dir data/predictions/ --rename TEST_MASTER_TABLE.csv
gpas download --filenames '.*mutations.csv' --output-dir data/mutations/ --rename TEST_MASTER_TABLE.csv
gpas download --filenames ".*variants.csv" --output-dir data/variants/ --rename TEST_MASTER_TABLE.csv
```

Each of these commands will download the specified output files (if created) for each of the samples. One of the TB samples has been purposefully included as it has a `main_report.json` but has insufficient reads to go to mapping so has no subsequent CSV files. Likewise one of the NTM samples is a mixture and so has two mutations and two variants CSV files. Obviously on a very large dataset these download commands would take some time.

### Using `gat` to create the data tables


Now we can use `gat` to discover and aggregate the files we've just downloaded. The important point is that what started as the mapping file (i.e. `TEST_MASTER_TABLE.csv`) now becomes the "lookup" master table that tells us which samples have worked (and likewise the `SPECIES` table tells us which species have mutations and variants CSV files). We can therefore use these tables to track which samples are failing and at which point so they can be re-run. The nomenclature we've used is any column destined for the final data tables is named in UPPERCASE whilst any column used to track progress is in lowercase. Because the table is **changed** by the `gat` command, you have to be careful what you are doing and it is strongly recommended to keep snapshots.

This is a bit hacky but `gpas` expects a column called `sample_name` whilst `gat` expects it to be called `RUN_ACCESSION` so first let's *manually* rename it in `TEST_MASTER_TABLE.csv` and then we can use `gat` to discover the files and build the tables. It is important that you first discover the `main_report.json` files as this creates the `SPECIES` table which we need for tracking in the subsequent steps.

```
gat build-tables --lookup-table TEST_MASTER_TABLE_renamed.csv --source-files data/ --output tables/ --tablename species
```

That has written two tables to `tables/`. One called `SPECIES` and the other is a modified copy of `TEST_MASTER_TABLE_renamed`. Both are written as CSV and Parquet files. Now we can issue

```
gat build-tables --lookup-table tables/TEST_MASTER_TABLE_renamed.csv --source-files data/ --output tables/ --tablename predictions
gat build-tables --lookup-table tables/TEST_MASTER_TABLE_renamed.csv --source-files data/ --output tables/ --tablename effects
```

These are both quick. The next two will take much longer using the defaults. To speed things up I've used `pandarallel` to parallelise the pandas `apply` function that is taking the time. Also to avoid memory issues the `chunks` argument breaks up the dataframe into that many chunks and processes them sequentially. The two interact so for `pandarallel` to run efficiently you want the smallest number of chunks as the time taken to load everything into memory slows it down, but if the chunks get too large then you get memory issues. The number of `chunks` defaults to 1 whilst `pandarallel` will grab all the CPU cores it can see. (Aside: I did try making the number of cores an argument but the way `pandarallel` is written makes this difficult -- you can of course change it but you need to modify the `pandarallel.initialize()` call in the code).

```
gat build-tables --lookup-table tables/TEST_MASTER_TABLE_renamed.csv --source-files data/ --output tables/ --tablename mutations
gat build-tables --lookup-table tables/TEST_MASTER_TABLE_renamed.csv --source-files data/ --output tables/ --tablename variants
```

The last thing to do is to join `MUTATIONS` and `VARIANTS` so we can populate `COVERAGE` and `FRS`.

```
gat correct-tables --input-dir tables/ --output-name tables/MUTATIONS-CORRECTED.parquet
```

Now we have a complete set of tables!

### Exploring the tables

There is a short exploration in `walkthrough-tables.ipynb` that talks about a few things. Load that notebook using e.g. VS Code.
