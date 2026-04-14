import pandas
import json

from tqdm import tqdm

tqdm.pandas()


def tidy_uid(row):
    if "main" in row.RUN_ACCESSION:
        return row.RUN_ACCESSION.split(".main")[0]
    else:
        return row.RUN_ACCESSION


def tidy_species(row):
    if "tuberculosis" in row.MAPPED_SPECIES:
        species = "Mycobacterium tuberculosis"
        subspecies = row.MAPPED_SPECIES.split(" (")[0].replace("M. ", "Mycobacterium ")
        lineage = row.LINEAGE
        sublineage = row.SUBLINEAGE
    elif isinstance(row.LINEAGE, str) and "subsp." in row.LINEAGE:
        species = row.MAPPED_SPECIES  # + ' complex'
        subspecies = row.LINEAGE
        lineage = None
        sublineage = None
    else:
        species = row.MAPPED_SPECIES
        subspecies = None
        lineage = row.LINEAGE
        sublineage = row.SUBLINEAGE
    return pandas.Series(
        {
            "SPECIES": species,
            "SUBSPECIES": subspecies,
            "LINEAGE": lineage,
            "SUBLINEAGE": sublineage,
        }
    )


def split_species(row):
    cols = row["name"].split(" (")
    species = cols[0]
    return species


def get_species_number(data, uid):
    # let's get a list of the mapped species names
    #  we assume here that if it doesn't have a new block there will be only one species!
    has_new_block = False
    if "Assembled NTM Results" in data:
        species_number = len(data["Assembled NTM Results"]["Assembled Species"])
        has_new_block = True
    else:
        if len(data["Mycobacterium Results"]["Species"]) != 1:
            print(
                f"Expected only one species but found {len(data['Mycobacterium Results']['Species'])} for {uid}"
            )
        species_number = 1

    return (species_number, has_new_block)


def parse_lineage(lineage_results):
    if lineage_results["Name"][:7] == "lineage":
        if lineage_results["Name"][7] in [
            "B",
            "C",
        ]:
            lineage = lineage_results["Name"]
        else:
            lineage = lineage_results["Name"][:8]
        sublineage = lineage_results["Name"]
    else:
        lineage = lineage_results["Name"].replace("_", " ")
        sublineage = ""

    return lineage, sublineage


def build_species_table(
    data_path=None,
    tables_path=None,
    master_table=None,
    max_samples=None,
):
    successful_genome = 0
    too_few_reads_genome = 0
    too_few_reads_id = 0

    tables = []
    n_samples = 0

    # discover all the main_report.json files and loop through them
    for i in tqdm((data_path).rglob("*main_report*.json")):
        uid = i.stem.split("_")[0]
        if "main" in uid:
            uid = uid.split(".main")[0]

        # let's check the uid is at least in the master_table!
        assert uid in master_table.index, f"UID {uid} not found in master table"

        # record that it has a main_report.json file
        master_table.at[uid, "has_main_report"] = True

        f = open(i)
        data = json.load(f)
        n_samples += 1
        if max_samples is not None and n_samples > max_samples:
            print("stopping")
            break

        if "Mycobacterial species identified" in data["Pipeline Outcome"]:
            too_few_reads_genome += 1
            master_table.at[uid, "status"] = "cannot assemble"
            # continue
        elif "Number of Mycobacterial reads is too low" in data["Pipeline Outcome"]:
            too_few_reads_id += 1
            master_table.at[uid, "status"] = "cannot speciate"
            continue
        elif "Sufficient reads mapped" in data["Pipeline Outcome"]:
            master_table.at[uid, "status"] = "complete"
            successful_genome += 1
        else:
            print(uid, data["Pipeline Outcome"])
            continue

        master_table.at[uid, "HUMAN_READS"] = data["Organism Identification"][
            "Human Reads"
        ]
        master_table.at[uid, "MYCOBACTERIAL_READS"] = data["Organism Identification"][
            "Mycobacterium Reads"
        ]
        master_table.at[uid, "OTHER_BACTERIAL_READS"] = data["Organism Identification"][
            "Non-Mycobacterium Bacteria Reads"
        ]
        master_table.at[uid, "UNCLASSIFIED_READS"] = data["Organism Identification"][
            "Unclassified Reads"
        ]

        pipeline_build = data["Metadata"]["Pipeline build"].replace("\n", "")
        master_table.at[uid, "GPAS_PIPELINE_BUILD"] = pipeline_build

        number_species, has_new_block = get_species_number(data, uid)
        master_table.at[uid, "has_new_block_in_main_report"] = has_new_block
        master_table.at[uid, "NUMBER_OF_SPECIES"] = number_species

        for sn in range(number_species):
            row = [uid]
            if has_new_block:
                row.append(data["Assembled NTM Results"]["Summary"][sn]["Name"])
                row.append(data["Assembled NTM Results"]["Summary"][sn]["Num Reads"])
                row.append(data["Assembled NTM Results"]["Summary"][sn]["Coverage"])
                row.append(data["Assembled NTM Results"]["Summary"][sn]["Depth"])
                if sn >= len(data["Assembled NTM Results"]["Lineage"]):
                    lineage, sublineage = "", ""
                else:
                    lineage, sublineage = parse_lineage(
                        data["Assembled NTM Results"]["Lineage"][sn]
                    )
            else:
                row.append(data["Mycobacterium Results"]["Summary"][sn]["Name"])
                row.append(data["Mycobacterium Results"]["Summary"][sn]["Num Reads"])
                row.append(data["Mycobacterium Results"]["Summary"][sn]["Coverage"])
                row.append(data["Mycobacterium Results"]["Summary"][sn]["Depth"])
                lineage, sublineage = parse_lineage(
                    data["Mycobacterium Results"]["Lineage"][sn]
                )
            row.append(lineage)
            row.append(sublineage)
            tables.append(row)

    SPECIES = pandas.DataFrame(
        tables,
        columns=[
            "RUN_ACCESSION",
            "MAPPED_SPECIES",
            "N_READS",
            "COVERAGE",
            "DEPTH",
            "LINEAGE",
            "SUBLINEAGE",
        ],
    )
    SPECIES[["SPECIES_NAME", "SUBSPECIES_NAME", "LINEAGE", "SUBLINEAGE"]] = (
        SPECIES.apply(tidy_species, axis=1)
    )
    SPECIES.drop(columns=["MAPPED_SPECIES"], inplace=True)
    SPECIES = SPECIES[
        [
            "RUN_ACCESSION",
            "SPECIES_NAME",
            "SUBSPECIES_NAME",
            "LINEAGE",
            "SUBLINEAGE",
            "N_READS",
            "COVERAGE",
            "DEPTH",
        ]
    ]
    SPECIES["RUN_ACCESSION"] = SPECIES.apply(tidy_uid, axis=1)

    for col in [
        "SPECIES_NAME",
        "SUBSPECIES_NAME",
        "LINEAGE",
        "SUBLINEAGE",
    ]:
        SPECIES[col] = SPECIES[col].astype("category")
    for col in ["N_READS"]:
        SPECIES[col] = SPECIES[col].astype("Int64")
    for col in ["COVERAGE", "DEPTH"]:
        SPECIES[col] = SPECIES[col].astype("float64")

    SPECIES.set_index("RUN_ACCESSION", inplace=True)
    SPECIES.to_csv(tables_path / "SPECIES.csv")
    SPECIES.to_parquet(tables_path / "SPECIES.parquet")

    total_genomes = successful_genome + too_few_reads_genome + too_few_reads_id
    print(f"{total_genomes} samples were processed.")
    print(
        f"{successful_genome} successfully reached a genome, {too_few_reads_id} samples had too few reads for identification and {too_few_reads_genome} samples had too few reads for genome assembly"
    )

    with pandas.option_context("future.no_silent_downcasting", True):
        master_table = master_table.fillna(
            value={
                "has_new_block_in_main_report": False,
                "NUMBER_OF_SPECIES": 0,
                "has_main_report": False,
            }
        ).infer_objects(copy=False)

    return master_table
