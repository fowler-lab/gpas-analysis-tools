import pandas, pathlib, json, numpy, glob

from tqdm import tqdm

from collections import defaultdict

from .species import build_species_table
from .genetics import build_genetics_table

tqdm.pandas()


def split_species(row):
    cols = row["name"].split(" (")
    species = cols[0]
    return species




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
    named_run_accession: bool = False,
    filename: str = None,
    chunks: int = 100,
):
    master_file = pathlib.Path(lookup_table)
    master_table = pandas.read_csv(master_file)
    master_table.set_index("RUN_ACCESSION", inplace=True)

    data_path = pathlib.Path(source_files)
    tables_path = pathlib.Path(output)

    assert filename in [
        "effects",
        "variants",
        "mutations",
        "predictions",
        "species",
    ], "can only specify one from this list"

    if filename == "species":
        
        master_table = build_species_table(data_path, tables_path, master_table, max_samples)

        master_table.to_csv(tables_path / (master_file.stem + ".csv"), index=True)
        master_table.to_parquet(tables_path / (master_file.stem + ".parquet"), index=True)

    if filename in ["effects", "mutations", "predictions", "variants"]:

        species_table = build_genetics_table(filename, data_path, tables_path, master_table, max_samples, chunks, False)

        # species_table = pandas.read_parquet(tables_path / "SPECIES.parquet")

            
        species_table.to_csv(tables_path / "SPECIES.csv", index=True)
        species_table.to_parquet(tables_path / "SPECIES.parquet")


def main():
    import defopt

    defopt.run(
        [build_tables, correct_tables],
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    print("boo")
    main()
