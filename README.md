# gpas-build-tables
CLI to build CRyPTIC data tables from the outputs of GPAS

Philip Fowler, 11 March 2026


Notes

1. `build-tables` should not modify the data fields and be able to cope with either JSON or CSV files downloaded from GPAS. The glob should be recursive to deal with a sharded file system
2. `correct-tables`. This is a separate command that takes the output of `build-tables` and joins `mutations` to `variants` to pick up some useful genetics metrics.
3. `shard-files` takes a delimiter (for CRyPTIC this is `.`) and moves the output files into a sharded file system (not yet implemented)

Issues

1. Why does `variants.parquet` fail to run i.e. is `Killed`? Is it a column type?
2. Why does `genomes` have more rows than it should? -> because of a carriage return in `pipeline_build`.
3. There are multiple ENA run accessions for some UNIQUEIDs; how do I know which was used?
4. Why does mykrobe report a lineage but then record zero median depth for some samples?
5. Why are some samples missing a `main_report`? Example is `site.10.subj.YA00040368.lab.YA00040368.iso.1` / `98bc5c23-d219-43bb-9aab-e8df1c6a0f7e` which I can download via mapping but isn't in the folder, suggesting it failed. Curiously it is the last file in the mapping csv.

## Setting up the GPAS CLI

The below sets up the GPAS CLI using a branch that contains a fix for downloading files using regexs.

```
git clone git@github.com:GlobalPathogenAnalysisService/client.git
cd client
git checkout fix/regex-download
python -m virtualenv env
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



First we need to instantiate the `gpas` conda environment.


