import pandas, pathlib, json, numpy, glob

from tqdm import tqdm

from collections import defaultdict

import pyarrow.parquet as pq

tqdm.pandas()


def split_species(row):
    cols = row["name"].split(" (")
    species = cols[0]
    return species


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


def correct_tables(input_dir: str = ".", output_dir: str = "output/"):

    assert input_dir != output_dir, "Input and output directories must be different!"

    input = pathlib.Path(input_dir)
    output = pathlib.Path(output_dir)

    print("reading the VARIANTS dataframe")
    variants = pandas.read_parquet(input / "VARIANTS.parquet")
    variants.reset_index(inplace=True)

    print("aggregating..")
    df = (
        variants[["UNIQUEID", "GENE", "GENE_POSITION", "COVERAGE"]]
        .groupby(["UNIQUEID", "GENE", "GENE_POSITION"], observed=True)
        .max()
    )

    print("reading the MUTATIONS dataframe")
    mutations = pandas.read_parquet(input / "MUTATIONS.parquet")
    mutations.reset_index(inplace=True)
    mutations.set_index(["UNIQUEID", "GENE", "GENE_POSITION"], inplace=True)

    print("joining..")
    mutations = mutations.join(df, how="left")

    mutations["FRS"] = mutations["MINOR_READS"] / mutations["COVERAGE"]

    print("correcting..")
    mask = (mutations.IS_MINOR) & (mutations.FRS >= 0.9)
    mutations.loc[mask, "MUTATION"] = mutations[mask]["MINOR_MUTATION"]
    mutations.loc[mask, "MINOR_MUTATION"] = None
    mutations.loc[mask, "MINOR_READS"] = numpy.nan
    mutations.loc[mask, "IS_MINOR"] = False

    assert (
        len(mutations[mutations.IS_MINOR & (mutations.FRS >= 0.9)]) == 0
    ), "Some mutations were not corrected"

    print("writing the new MUTATIONS dataframe to disc")
    mutations.reset_index(inplace=True)
    mutations.set_index(["UNIQUEID", "GENE", "MUTATION"], inplace=True)
    mutations.to_parquet(output / "MUTATIONS.parquet")


def build_tables(
    lookup_table: str = None,
    source_files: str = "data/",
    max_samples: int = None,
    output: str = None,
    uppercase: bool = True,
    named_run_accession: bool = True,
    filename: str = None,
    chunks: int = 100,
):
    master_file = pathlib.Path(lookup_table)
    master_table = pandas.read_csv(master_file)
    master_table.set_index("UNIQUEID", inplace=True)

    path = pathlib.Path(source_files)
    tables_path = pathlib.Path(output)

    assert filename in [
        "effects",
        "variants",
        "mutations",
        "predictions",
        "genomes",
    ], "can only specify one from this list"

    if filename in ["effects", "mutations", "predictions", "variants"]:

        tables = []
        n_samples = 0

        n_files = sum(1 for i in (path).rglob("*" + filename + ".csv"))

        for i in tqdm((path).rglob("*" + filename + ".csv"), total=n_files):

            df = pandas.read_csv(i)

            n_samples += 1
            if max_samples is not None and n_samples > max_samples:
                break

            if named_run_accession:
                ena_run_accession = i.stem.split("." + filename)[0]
                uid = master_table[
                    master_table.run_accession == ena_run_accession
                ].index.values[0]
                df["uniqueid"] = uid
            else:
                uid = i.stem.split("." + filename)[0]
                df["uniqueid"] = uid
            master_table.at[uid, "has_" + filename] = True
            tables.append(df)

        df = pandas.concat(tables)

        if filename == "effects":
            for col in [
                "gene",
                "drug",
                "prediction",
                "catalogue_name",
                "prediction_values",
            ]:
                df[col] = df[col].astype("category")
            if uppercase:
                df = df.rename(columns=str.upper)
                df.set_index(
                    [
                        "UNIQUEID",
                        "CATALOGUE_NAME",
                        "CATALOGUE_VERSION",
                        "PREDICTION_VALUES",
                        "DRUG",
                        "GENE",
                        "MUTATION",
                    ],
                    inplace=True,
                )
            else:
                df.set_index(
                    [
                        "uniqueid",
                        "catalogue_name",
                        "catalogue_version",
                        "prediction_values",
                        "drug",
                        "gene",
                        "mutation",
                    ],
                    inplace=True,
                )

        elif filename == "predictions":
            for col in [
                "drug",
                "prediction",
                "catalogue_name",
                "catalogue_values",
            ]:
                df[col] = df[col].astype("category")

            if uppercase:
                df = df.rename(columns=str.upper)
                df.set_index(
                    [
                        "UNIQUEID",
                        "CATALOGUE_NAME",
                        "CATALOGUE_VERSION",
                        "CATALOGUE_VALUES",
                        "DRUG",
                    ],
                    inplace=True,
                )

            else:
                df.set_index(
                    [
                        "uniqueid",
                        "catalogue_name",
                        "catalogue_version",
                        "catalogue_values",
                        "drug",
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
                ] = df_i.apply(parse_variants, axis=1)
                df_i.drop(columns=["variant"], inplace=True)
                df_i.rename(columns={"var": "variant"}, inplace=True)
                for col in [
                    "gene",
                ]:
                    df_i[col] = df_i[col].astype("category")

                if uppercase:
                    df_i = df_i.rename(columns=str.upper)
                    df_i.set_index(["UNIQUEID", "GENE", "VARIANT"], inplace=True)
                    df_i.to_csv(
                        str(tables_path / (filename.upper() + "_" + str(counter)))
                        + ".csv"
                    )
                    df_i.drop(columns=["VCF_EVIDENCE"], inplace=True)
                    df_i.to_parquet(
                        str(tables_path / (filename.upper() + "_" + str(counter)))
                        + ".parquet"
                    )
                else:
                    df_i.set_index(["uniqueid", "gene", "variant"], inplace=True)
                    df_i.to_csv(
                        str(tables_path / (filename + "_" + str(counter))) + ".csv"
                    )
                    df_i.drop(columns=["vcf_evidence"], inplace=True)
                    df_i.to_parquet(
                        str(tables_path / (filename + "_" + str(counter))) + ".parquet"
                    )

                tables.append(df_i)
                counter += 1
            df = pandas.concat(tables)
            if uppercase:
                files = glob.glob(str(tables_path) + "/VARIANTS_*.parquet")
                schema = pq.ParquetFile(files[0]).schema_arrow
                with pq.ParquetWriter(
                    str(tables_path) + "/VARIANTS.parquet", schema=schema
                ) as writer:
                    for file in tqdm(files):
                        writer.write_table(pq.read_table(file, schema=schema))
            else:
                files = glob.glob(str(tables_path) + "/variants_*.parquet")
                schema = pq.ParquetFile(files[0]).schema_arrow
                with pq.ParquetWriter(
                    str(tables_path) + "/variants.parquet", schema=schema
                ) as writer:
                    for file in tqdm(files):
                        writer.write_table(pq.read_table(file, schema=schema))

        elif filename == "mutations":
            tables = []
            for df_i in tqdm(numpy.array_split(df, chunks)):
                df_i[
                    ["mut", "is_null", "is_minor", "minor_mutation", "minor_reads"]
                ] = df_i.apply(parse_mutations, axis=1)
                df_i.drop(columns=["mutation"], inplace=True)
                df_i.rename(columns={"mut": "mutation"}, inplace=True)
                for col in [
                    "gene",
                    "ref",
                    "alt",
                    "amino_acid_sequence",
                    "minor_mutation",
                ]:
                    df_i[col] = df_i[col].astype("category")

                if uppercase:
                    df_i = df_i.rename(columns=str.upper)
                    df_i.set_index(["UNIQUEID", "GENE", "MUTATION"], inplace=True)
                else:
                    df_i.set_index(["uniqueid", "gene", "mutation"], inplace=True)
                tables.append(df_i)
            df = pandas.concat(tables)

        if filename != "variants":
            if uppercase:
                df.to_csv(str(tables_path / filename.upper()) + ".csv", index=True)
                df.to_parquet(str(tables_path / filename.upper()) + ".parquet")
            else:
                df.to_csv(str(tables_path / filename) + ".csv", index=True)
                df.to_parquet(str(tables_path / filename) + ".parquet")

    successful_genome = 0
    too_few_reads_genome = 0
    too_few_reads_id = 0
    if filename == "genomes":

        tables = []
        n_samples = 0

        for i in tqdm((path).rglob("*main_report.json")):

            if named_run_accession:
                ena_run_accession = i.stem.split(".main_report")[0]
                uid = master_table[
                    master_table.run_accession == ena_run_accession
                ].index.values[0]
            else:
                uid = i.stem.split(".main_report")[0]

            master_table.at[uid, "has_genome"] = True

            row = [uid]
            f = open(i)
            data = json.load(f)
            n_samples += 1
            if max_samples is not None and n_samples > max_samples:
                break

            if (
                data["Pipeline Outcome"]
                == "Mycobacterial species identified. Reads mapped to M. tuberculosis (H37Rv v3) too low to proceed to M. tuberculosis complex genome assembly."
            ):
                too_few_reads_genome += 1
                master_table.at[uid, "status"] = "cannot assemble"
                continue
            elif (
                data["Pipeline Outcome"]
                == "Number of Mycobacterial reads is too low to proceed to Mycobacterial species identification."
            ):
                too_few_reads_id += 1
                master_table.at[uid, "status"] = "cannot speciate"
                continue
            elif (
                data["Pipeline Outcome"]
                == "Sufficient reads mapped to M. tuberculosis (H37Rv v3) for genome assembly, resistance prediction and relatedness assessment."
            ):
                master_table.at[uid, "status"] = "complete"
                successful_genome += 1
            else:
                print(uid, data["Pipeline Outcome"])
                continue

            row.append(data["Organism Identification"]["Mycobacterium Reads"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Name"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Num Reads"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Coverage"])
            row.append(data["Mycobacterium Results"]["Summary"][0]["Depth"])
            lineage_results = data["Mycobacterium Results"]["Lineage"]
            n_lineages = len(lineage_results)
            if n_lineages == 1:
                if lineage_results[0]["Name"][:7] == "lineage":
                    if lineage_results[0]["Name"][7] in [
                        "B",
                        "C",
                    ]:
                        lineage = lineage_results[0]["Name"]
                    else:
                        lineage = lineage_results[0]["Name"][:8]
                    sublineage = lineage_results[0]["Name"]
                    # lineage_cov = lineage_results[0]["Coverage"]
                    # lineage_depth = lineage_results[0]["Median Depth"]
                else:
                    lineage = lineage_results[0]["Name"]
                    sublineage = ""
                    # lineage_cov = lineage_results[0]["Coverage"]
                    # lineage_depth = lineage_results[0]["Median Depth"]

            else:
                lineage = "mixed"
                sublineage = ""
                # lineage_cov = None
                # lineage_depth = None
                for i in lineage_results:
                    sublineage += i["Name"] + "/"
                sublineage = sublineage[:-1]
            row.append(n_lineages)
            row.append(lineage)
            row.append(sublineage)
            # row.append(lineage_cov)
            # row.append(lineage_depth)

            antibiogram = ""
            amr_results = data["Genomes"][0]["Resistance Prediction"][
                "Resistance Prediction Summary"
            ]
            for j in amr_results:
                for k in amr_results[j]:
                    antibiogram += amr_results[j][k]
            row.append(antibiogram)

            pipeline_build = data["Metadata"]["Pipeline build"].replace("\n", "")
            row.append(pipeline_build)

            tables.append(row)

        genomes = pandas.DataFrame(
            tables,
            columns=[
                "uniqueid",
                "mycobacterial_reads",
                "name",
                "tb_reads",
                "tb_coverage",
                "tb_depth",
                "n_lineages",
                "lineage",
                "sublineage",
                "antibiogram",
                "pipeline_build",
            ],
        )
        genomes["species"] = genomes.apply(split_species, axis=1)
        genomes.drop(columns=["name"], inplace=True)
        genomes = genomes[
            [
                "uniqueid",
                "species",
                "n_lineages",
                "lineage",
                "sublineage",
                "mycobacterial_reads",
                "tb_reads",
                "tb_coverage",
                "tb_depth",
                "antibiogram",
                "pipeline_build",
            ]
        ]
        genomes["lineage"] = genomes["lineage"].astype("category")
        genomes["sublineage"] = genomes["sublineage"].astype("category")
        genomes["antibiogram"] = genomes["antibiogram"].astype("category")
        genomes["pipeline_build"] = genomes["pipeline_build"].astype("category")
        genomes.fillna(value={"status": "not processed"}, inplace=True)

        if uppercase:
            genomes = genomes.rename(columns=str.upper)
            genomes.set_index("UNIQUEID", inplace=True)
            genomes.to_csv(tables_path / "GENOMES.csv")
            genomes.to_parquet(tables_path / "GENOMES.parquet")
        else:
            genomes.set_index("uniqueid", inplace=True)
            genomes.to_csv(tables_path / "genomes.csv")
            genomes.to_parquet(tables_path / "genomes.parquet")

        total_genomes = successful_genome + too_few_reads_genome + too_few_reads_id
        print(f"{total_genomes} samples were processed.")
        print(
            f"{successful_genome} successfully reached a genome, {too_few_reads_id} samples had too few reads for identification and {too_few_reads_genome} samples had too few reads for genome assembly"
        )
    master_output = master_file.name.split(".")[0]
    master_table.to_csv(tables_path / (master_output + ".csv"), index=True)
    master_table.to_parquet(tables_path / (master_output + ".parquet"))


def main():
    import defopt

    defopt.run(
        [build_tables, correct_tables],
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
