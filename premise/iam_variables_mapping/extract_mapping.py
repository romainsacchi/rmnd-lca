import glob

import pandas as pd
import yaml

# Initialize a dictionary to store DataFrames
dfs = {}


# Function to process each YAML file and extract data into DataFrame
def process_yaml(file_path):
    with open(file_path, "r") as file:
        data = yaml.safe_load(file)
    df_data = []

    def append_data(key, model, variable, variable_type):
        if isinstance(variable, list):
            for var in variable:
                df_data.append([key, model, var, variable_type])
        else:
            df_data.append([key, model, variable, variable_type])

    for main_key, value in data.items():
        if isinstance(value, dict):
            if "iam_aliases" in value:
                for model, variable in value["iam_aliases"].items():
                    append_data(main_key, model, variable, "Production volume")
            if "energy_use_aliases" in value:
                for model, variable in value["energy_use_aliases"].items():
                    append_data(
                        f"{main_key} (energy use)",
                        model,
                        variable,
                        "Specific Energy Use",
                    )
            if "eff_aliases" in value:
                for model, variable in value["eff_aliases"].items():
                    append_data(
                        f"{main_key} (efficiency)", model, variable, "Efficiency"
                    )
            if "electricity_use_aliases" in value:
                for model, variable in value["electricity_use_aliases"].items():
                    append_data(
                        f"{main_key} (electricity use)",
                        model,
                        variable,
                        "Electricity Use",
                    )
            if "heat_use_aliases" in value:
                for model, variable in value["heat_use_aliases"].items():
                    append_data(f"{main_key} (heat use)", model, variable, "Heat Use")
            if "land_use" in value:
                for model, variable in value["land_use"].items():
                    append_data(main_key, model, variable, "Land Use")
            if "land_use_change" in value:
                for model, variable in value["land_use_change"].items():
                    append_data(main_key, model, variable, "Land Use Change")

            for sub_key, sub_value in value.items():
                if isinstance(sub_value, dict):
                    if "iam_aliases" in sub_value:
                        for model, variable in sub_value["iam_aliases"].items():
                            append_data(
                                f"{main_key} - {sub_key}",
                                model,
                                variable,
                                "Production volume",
                            )
                    if "energy_use_aliases" in sub_value:
                        for model, variable in sub_value["energy_use_aliases"].items():
                            append_data(
                                f"{main_key} - {sub_key} (energy use)",
                                model,
                                variable,
                                "Specific Energy Use",
                            )
                    if "eff_aliases" in sub_value:
                        for model, variable in sub_value["eff_aliases"].items():
                            append_data(
                                f"{main_key} - {sub_key} (efficiency)",
                                model,
                                variable,
                                "Efficiency",
                            )
                    if "electricity_use_aliases" in sub_value:
                        for model, variable in sub_value[
                            "electricity_use_aliases"
                        ].items():
                            append_data(
                                f"{main_key} - {sub_key} (electricity use)",
                                model,
                                variable,
                                "Electricity Use",
                            )
                    if "heat_use_aliases" in sub_value:
                        for model, variable in sub_value["heat_use_aliases"].items():
                            append_data(
                                f"{main_key} - {sub_key} (heat use)",
                                model,
                                variable,
                                "Heat Use",
                            )
                    if "land use" in sub_value:
                        append_data(
                            f"{main_key} - {sub_key}",
                            "IAM",
                            sub_value["land use"],
                            "Land Use",
                        )
                    if "land_use_change" in sub_value:
                        append_data(
                            f"{main_key} - {sub_key}",
                            "IAM",
                            sub_value["land_use_change"],
                            "Land Use Change",
                        )
    df = pd.DataFrame(
        df_data, columns=["Key", "IAM Model", "Variable", "Variable Type"]
    )
    return df


# Process all *.yaml files except "constants.yaml"
for file in glob.glob("*.yaml"):
    if file in [
        "constants.yaml",
        "missing_geography_equivalences.yaml",
        "gains_regions_mapping.yaml",
        "iam_region_to_climate.yaml",
    ]:
        continue
    file_name = (
        file.split("/")[-1].replace("_variables.yaml", "").replace("_", " ").title()
    )
    dfs[file_name] = process_yaml(file)

# Create a Pandas Excel writer
excel_file = "mapping_overview.xlsx"
with pd.ExcelWriter(excel_file, engine="xlsxwriter") as writer:
    for sheet_name, df in dfs.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False)


def create_refined_pivot_v2(df):
    df = df.copy()

    # Identify and extract the terms in parentheses
    df["Energy use"] = df["Key"].str.extract(r"\((energy use)\)", expand=False)
    df["Efficiency"] = df["Key"].str.extract(r"\((efficiency)\)", expand=False)

    # Clean the 'Key' column by removing the extracted terms
    df["Key"] = df["Key"].str.replace(r"\s*\(energy use\)\s*", "", regex=True)
    df["Key"] = df["Key"].str.replace(r"\s*\(efficiency\)\s*", "", regex=True)

    # Fill NaN with empty strings for consistency
    df.fillna("", inplace=True)

    # Combine 'Energy use' and 'Efficiency' into the 'Variable Type' column for pivoting
    df["Variable Category"] = (
        df["Variable Type"] + " " + df["Energy use"] + df["Efficiency"]
    )
    df["Variable Category"] = df["Variable Category"].str.strip()

    # Create the pivot table with combined variables separated by ";"
    pivot_df = df.pivot_table(
        index=["Key"],
        columns=["IAM Model", "Variable Category"],
        values="Variable",
        aggfunc=lambda x: "; ".join(str(v) for v in x),
    )

    # Clean up the column headers to remove any trailing "|"
    pivot_df.columns = [
        col[0] + " | " + col[1].strip(" |") if col[1] else col[0]
        for col in pivot_df.columns
    ]
    pivot_df.columns = pivot_df.columns.str.replace(r"\s+\|\s+$", "", regex=True)

    return pivot_df.T


# Apply the refined function to each DataFrame and create pivot tables
refined_pivot_tables_v2 = {}
for sheet, df in dfs.items():
    if not df.empty:
        refined_pivot_tables_v2[sheet] = create_refined_pivot_v2(df)

# Create a new Excel writer for the refined pivot tables
refined_output_path_v2 = "mapping_overview2.xlsx"
with pd.ExcelWriter(refined_output_path_v2) as writer:
    for sheet, pivot_df in refined_pivot_tables_v2.items():
        if not pivot_df.empty:
            pivot_df.to_excel(writer, sheet_name=sheet)
