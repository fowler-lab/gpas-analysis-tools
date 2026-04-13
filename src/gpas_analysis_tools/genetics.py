import pandas, numpy, glob, json

import pyarrow.parquet as pq

from tqdm import tqdm

tqdm.pandas()

from pandarallel import pandarallel

pandarallel.initialize(progress_bar=False, nb_workers=16)

def parse_variants(row):

    variant = row.variant
    is_null = False
    is_minor = False
    minor_variant = None
    minor_reads = 0
    coverage = 0

    a = json.loads(row["vcf_evidence"])

    # From Jeremy:
    #  "Just heard back from Martin, because of minos mapping
    #   against a graph rather than genome, COV_TOTAL is reads
    #   mapping to the site (so reads can map to multiple sites
    #   in an allele), whereas COV is unique reads mapping to
    #   the allele, which is why COV_TOTAL can be bigger. So
    #   I think sum(COV) is what we should be using here"

    # if "COV_TOTAL" in a.keys():
    #     coverage = a["COV_TOTAL"]
    if "COV" in a.keys():
        if isinstance(a["COV"], list):
            coverage = sum(a["COV"])
        elif isinstance(a["COV"], int):
            coverage = a["COV"]
        else:
            coverage = -1

    # FRS can be the wrong way round so ignore
    # if "FRS" in a.keys():
    #     frs = a["FRS"]
    #     if frs == ".":
    #         frs = None

    if variant[-1] == "x":
        is_null = True
    elif ":" in variant:
        is_minor = True
        cols = variant.split(":")
        minor_variant = cols[0]
        minor_reads = int(cols[1])
        if "ins" in variant or "del" in variant:
            variant = minor_variant.split("_")[0] + "_minorindel"
        else:
            variant = minor_variant[:-1] + "z"

    return pandas.Series(
        [variant, is_null, is_minor, minor_variant, minor_reads, coverage]
    )


def parse_mutations(row):

    mutation = row.mutation
    is_null = False
    is_minor = False
    minor_mutation = None
    minor_reads = None

    if ":" in mutation:
        is_minor = True
        cols = mutation.split(":")
        minor_mutation = cols[0]
        minor_reads = int(cols[1])
        if "ins" in mutation or "del" in mutation:
            mutation = mutation.split("_")[0] + "_minorindel"
        else:
            if minor_mutation[0] in ["a", "t", "c", "g"]:
                mutation = minor_mutation[:-1] + "z"
            else:
                mutation = minor_mutation[:-1] + "Z"
    elif mutation[-1] in ["X", "x"]:
        is_null = True

    return pandas.Series([mutation, is_null, is_minor, minor_mutation, minor_reads])


def build_genetics_table(filename, data_path, tables_path, master_table, max_samples, chunks, uppercase = True):

    tables = []
    n_samples = 0

    species_table = pandas.read_parquet(tables_path / "SPECIES.parquet")

    n_files = sum(1 for i in (data_path).rglob("*" + filename + "*.csv"))

    for i in tqdm((data_path).rglob("*" + filename + "*.csv"), total=n_files):

        n_samples += 1
        if max_samples is not None and n_samples > max_samples:
            break

        uid = i.stem.split("_")[0]
        if filename in uid:
            uid = uid.split('.'+filename)[0]

        # let's check the uid is at least in the master_table!
        assert uid in master_table.index, f"UID {uid} not found in master table"
        
        if not master_table.at[uid, "has_main_report"]:
            print(f"UID {uid} does not have main report, skipping")
            continue

        df = pandas.read_csv(i)
        
        # check to see if the same has the 'Assembled NTM Results' block in the main report
        if master_table.at[uid, 'has_new_block_in_main_report']:
            for sn in list(species_table.at[uid, 'SPECIES_NAME']):
                if sn.replace(' ','_') in i.stem:
                    species_name = sn
        else:
            species_name = species_table.at[uid, 'SPECIES_NAME']
        df.insert(1, 'species_name', species_name)

        species_table.at[uid, "has_" + filename] = True
        tables.append(df)

    df = pandas.concat(tables)

    if filename == "effects":
        for col in [
            "species_name",
            "gene",
            "drug",
            "prediction",
            "catalogue_name",
            "prediction_values",
        ]:
            df[col] = df[col].astype("category")

        df = df.rename(columns=str.upper)
        df.set_index(
            [
                "UNIQUEID",
                "SPECIES_NAME",
                "CATALOGUE_NAME",
                "CATALOGUE_VERSION",
                "PREDICTION_VALUES",
                "DRUG",
                "GENE",
                "MUTATION",
            ],
            inplace=True,
        )

    elif filename == "predictions":
        for col in [
            "species_name",
            "drug",
            "prediction",
            "catalogue_name",
            "catalogue_values",
        ]:
            df[col] = df[col].astype("category")

        df = df.rename(columns=str.upper)
        df.set_index(
            [
                "UNIQUEID",
                "SPECIES_NAME",
                "CATALOGUE_NAME",
                "CATALOGUE_VERSION",
                "CATALOGUE_VALUES",
                "DRUG",
            ],
            inplace=True,
        )

    elif filename == "variants":
        tables = []
        counter = 0
        for df_i in tqdm(numpy.array_split(df, chunks)):
            df_i[
                [
                    "var",
                    "is_null",
                    "is_minor",
                    "minor_variant",
                    "minor_reads",
                    "coverage",
                ]
            ] = df_i.parallel_apply(parse_variants, axis=1)
            df_i.drop(columns=["variant"], inplace=True)
            df_i.rename(columns={"var": "variant"}, inplace=True)
            for col in [
                "gene",
                "species_name",
            ]:
                df_i[col] = df_i[col].astype("category")

            df_i = df_i.rename(columns=str.upper)
            df_i.set_index(["UNIQUEID", "SPECIES_NAME", "GENE", "VARIANT"], inplace=True)
            df_i.to_csv(
                str(tables_path / (filename.upper() + "_" + str(counter)))
                + ".csv"
            )
            df_i.drop(columns=["VCF_EVIDENCE"], inplace=True)
            df_i.to_parquet(
                str(tables_path / (filename.upper() + "_" + str(counter)))
                + ".parquet"
            )

            tables.append(df_i)
            counter += 1
        df = pandas.concat(tables)

        files = glob.glob(str(tables_path) + "/VARIANTS_*.parquet")
        schema = pq.ParquetFile(files[0]).schema_arrow
        with pq.ParquetWriter(
            str(tables_path) + "/VARIANTS.parquet", schema=schema
        ) as writer:
            for file in tqdm(files):
                writer.write_table(pq.read_table(file, schema=schema))

    elif filename == "mutations":
        tables = []
        print('starting to parse mutations, this may take a while...')
        for df_i in tqdm(numpy.array_split(df, chunks)):
            df_i[
                ["mut", "is_null", "is_minor", "minor_mutation", "minor_reads"]
            ] = df_i.parallel_apply(parse_mutations, axis=1)
            df_i.drop(columns=["mutation"], inplace=True)
            df_i.rename(columns={"mut": "mutation"}, inplace=True)
            for col in [
                "species_name",
                "gene",
                "ref",
                "alt",
                "amino_acid_sequence",
                "minor_mutation",
            ]:
                df_i[col] = df_i[col].astype("category")

            df_i = df_i.rename(columns=str.upper)
            df_i.set_index(["UNIQUEID", "SPECIES_NAME", "GENE", "MUTATION"], inplace=True)

            tables.append(df_i)

        df = pandas.concat(tables)

    if filename != "variants":
        df.to_csv(str(tables_path / filename.upper()) + ".csv", index=True)
        df.to_parquet(str(tables_path / filename.upper()) + ".parquet")

    with pandas.option_context("future.no_silent_downcasting", True):
        species_table = species_table.fillna(value = {
                "has_" + filename: False
                }).infer_objects(copy=False)

    return species_table