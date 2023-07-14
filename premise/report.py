"""
This module export a summary of scenario to an Excel file.
"""

import os
from datetime import datetime
from pathlib import Path

import openpyxl
import pandas as pd
import xarray as xr
import yaml
from openpyxl.chart import AreaChart, LineChart, Reference
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.worksheet.dimensions import ColumnDimension, DimensionHolder
from pandas.errors import EmptyDataError

from . import DATA_DIR, VARIABLES_DIR, __version__

IAM_ELEC_VARS = VARIABLES_DIR / "electricity_variables.yaml"
IAM_FUELS_VARS = VARIABLES_DIR / "fuels_variables.yaml"
IAM_BIOMASS_VARS = VARIABLES_DIR / "biomass_variables.yaml"
IAM_CEMENT_VARS = VARIABLES_DIR / "cement_variables.yaml"
IAM_STEEL_VARS = VARIABLES_DIR / "steel_variables.yaml"
IAM_DACCS_VARS = VARIABLES_DIR / "direct_air_capture_variables.yaml"
IAM_OTHER_VARS = VARIABLES_DIR / "other_variables.yaml"
IAM_CARBON_CAPTURE_VARS = VARIABLES_DIR / "carbon_capture_variables.yaml"
REPORT_METADATA_FILEPATH = DATA_DIR / "utils" / "report" / "report.yaml"
VEHICLES_MAP = DATA_DIR / "transport" / "vehicles_map.yaml"

LOG_REPORTING_FILEPATH = DATA_DIR / "utils" / "logging" / "reporting.yaml"
# directory for log files
DIR_LOGS = Path.cwd() / "export" / "logs"
DIR_LOG_REPORT = Path.cwd() / "export" / "change reports"

# if DIR_LOG_REPORT folder does not exist
# we create it
if not Path(DIR_LOG_REPORT).exists():
    Path(DIR_LOG_REPORT).mkdir(parents=True, exist_ok=True)


def get_variables(
    filepath,
):
    """
    Get the variables from a yaml file.
    :param filepath: path to the yaml file
    """
    with open(filepath, "r", encoding="utf-8") as stream:
        out = yaml.safe_load(stream)

    return list(out.keys())


def fetch_data(
    iam_data: xr.DataArray, sector: str, variable: str
) -> [xr.DataArray, None]:
    data = {
        "Population": iam_data.other_vars if hasattr(iam_data, "other_vars") else None,
        "GDP": iam_data.other_vars if hasattr(iam_data, "other_vars") else None,
        "CO2": iam_data.other_vars if hasattr(iam_data, "other_vars") else None,
        "GMST": iam_data.other_vars if hasattr(iam_data, "other_vars") else None,
        "Electricity - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Electricity (biom) - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Electricity - efficiency": iam_data.electricity_efficiencies
        if hasattr(iam_data, "electricity_efficiencies")
        else None,
        "Fuel (gasoline) - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Fuel (gasoline) - efficiency": iam_data.petrol_efficiencies
        if hasattr(iam_data, "petrol_efficiencies")
        else None,
        "Fuel (diesel) - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Fuel (diesel) - efficiency": iam_data.diesel_efficiencies
        if hasattr(iam_data, "diesel_efficiencies")
        else None,
        "Fuel (gas) - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Fuel (gas) - efficiency": iam_data.gas_efficiencies
        if hasattr(iam_data, "gas_efficiencies")
        else None,
        "Fuel (hydrogen) - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Fuel (hydrogen) - efficiency": iam_data.hydrogen_efficiencies
        if hasattr(iam_data, "hydrogen_efficiencies")
        else None,
        "Cement - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Cement - efficiency": iam_data.cement_efficiencies
        if hasattr(iam_data, "cement_efficiencies")
        else None,
        "Cement - CCS": iam_data.carbon_capture_rate
        if hasattr(iam_data, "carbon_capture_rate")
        else None,
        "Steel - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Steel - efficiency": iam_data.steel_efficiencies
        if hasattr(iam_data, "steel_efficiencies")
        else None,
        "Steel - CCS": iam_data.carbon_capture_rate
        if hasattr(iam_data, "carbon_capture_rate")
        else None,
        "Direct Air Capture - generation": iam_data.production_volumes
        if hasattr(iam_data, "production_volumes")
        else None,
        "Transport (cars)": iam_data.trsp_cars
        if hasattr(iam_data, "trsp_cars")
        else None,
        "Transport (buses)": iam_data.trsp_buses
        if hasattr(iam_data, "trsp_buses")
        else None,
        "Transport (trucks)": iam_data.trsp_trucks
        if hasattr(iam_data, "trsp_trucks")
        else None,
    }

    if data[sector] is not None:
        iam_data = data[sector]

        if any(x in sector for x in ["car", "bus", "truck"]):
            iam_data = iam_data.sum(dim=["size", "construction_year"])
            iam_data = iam_data.rename({"powertrain": "variables"}).T

        return iam_data.sel(
            variables=[v for v in variable if v in iam_data.coords["variables"].values]
        )
    else:
        return None


def generate_summary_report(scenarios: list, filename: Path) -> None:
    """
    Generate a summary report of the scenarios.
    """

    SECTORS = {
        "Population": {
            "filepath": IAM_OTHER_VARS,
            "variables": [
                "population",
            ],
        },
        "GDP": {
            "filepath": IAM_OTHER_VARS,
            "variables": [
                "gdp",
            ],
        },
        "CO2": {
            "filepath": IAM_OTHER_VARS,
            "variables": [
                "CO2",
            ],
        },
        "GMST": {
            "filepath": IAM_OTHER_VARS,
            "variables": [
                "GMST",
            ],
        },
        "Electricity - generation": {
            "filepath": IAM_ELEC_VARS,
        },
        "Electricity (biom) - generation": {
            "filepath": IAM_BIOMASS_VARS,
        },
        "Electricity - efficiency": {
            "filepath": IAM_ELEC_VARS,
        },
        "Fuel (gasoline) - generation": {
            "filepath": IAM_FUELS_VARS,
            "filter": ["gasoline", "ethanol"],
        },
        "Fuel (gasoline) - efficiency": {
            "filepath": IAM_FUELS_VARS,
            "filter": ["gasoline", "ethanol"],
        },
        "Fuel (diesel) - generation": {
            "filepath": IAM_FUELS_VARS,
            "filter": [
                "diesel",
            ],
        },
        "Fuel (diesel) - efficiency": {
            "filepath": IAM_FUELS_VARS,
            "filter": [
                "diesel",
            ],
        },
        "Fuel (gas) - generation": {
            "filepath": IAM_FUELS_VARS,
            "filter": ["natural gas", "biogas", "methane"],
        },
        "Fuel (gas) - efficiency": {
            "filepath": IAM_FUELS_VARS,
            "filter": ["natural gas", "biogas", "methane"],
        },
        "Fuel (hydrogen) - generation": {
            "filepath": IAM_FUELS_VARS,
            "filter": [
                "hydrogen",
            ],
        },
        "Fuel (hydrogen) - efficiency": {
            "filepath": IAM_FUELS_VARS,
            "filter": [
                "hydrogen",
            ],
        },
        "Cement - generation": {
            "filepath": IAM_CEMENT_VARS,
        },
        "Cement - efficiency": {
            "filepath": IAM_CEMENT_VARS,
        },
        "Cement - CCS": {
            "filepath": IAM_CARBON_CAPTURE_VARS,
            "variables": ["cement"],
        },
        "Steel - generation": {
            "filepath": IAM_STEEL_VARS,
        },
        "Steel - efficiency": {
            "filepath": IAM_STEEL_VARS,
        },
        "Steel - CCS": {
            "filepath": IAM_CARBON_CAPTURE_VARS,
            "variables": ["steel"],
        },
        "Direct Air Capture - generation": {
            "filepath": IAM_DACCS_VARS,
            "variables": ["dac_solvent", "dac_sorbent"],
        },
        "Transport (cars)": {
            "filepath": VEHICLES_MAP,
            "variables": [
                "BEV",
                "FCEV",
                "ICEV-d",
                "ICEV-g",
                "ICEV-p",
                "PHEV-d",
                "PHEV-p",
            ],
        },
        "Transport (buses)": {
            "filepath": VEHICLES_MAP,
            "variables": [
                "BEV",
                "FCEV",
                "ICEV-d",
                "ICEV-g",
                "ICEV-p",
                "PHEV-d",
                "PHEV-p",
            ],
        },
        "Transport (trucks)": {
            "filepath": VEHICLES_MAP,
            "variables": [
                "BEV",
                "FCEV",
                "ICEV-d",
                "ICEV-g",
                "ICEV-p",
                "PHEV-d",
                "PHEV-p",
            ],
        },
    }

    with open(REPORT_METADATA_FILEPATH, "r", encoding="utf-8") as stream:
        metadata = yaml.safe_load(stream)

    workbook = openpyxl.Workbook()
    workbook.remove(workbook.active)

    for sector, filepath in SECTORS.items():
        if "variables" in filepath:
            variables = filepath["variables"]
        else:
            variables = get_variables(filepath["filepath"])
            if "filter" in filepath:
                variables = [
                    x for x in variables if any(y in x for y in filepath["filter"])
                ]

        worksheet = workbook.create_sheet(sector)

        col, row = (1, 1)

        worksheet.cell(
            column=col,
            row=row,
            value=metadata[sector]["expl_text"],
        )

        scenario_list = []

        last_col_used = 0

        for scenario_idx, scenario in enumerate(scenarios):
            if (scenario["model"], scenario["pathway"]) not in scenario_list:
                iam_data = fetch_data(
                    iam_data=scenario["iam data"],
                    sector=sector,
                    variable=variables,
                )

                if iam_data is None:
                    continue

                if scenario_idx > 0:
                    col = last_col_used + metadata[sector]["offset"]

                row = 3

                worksheet.cell(
                    column=col,
                    row=row,
                    value=f"{scenario['model'].upper()} - {scenario['pathway'].upper()}",
                )
                worksheet.cell(column=col, row=row).font = Font(
                    bold=True, size=14, underline="single"
                )

                row += 2

                for region in scenario["iam data"].regions:
                    if sector == "GMST" and region != "World":
                        continue

                    worksheet.cell(column=col, row=row, value=region)

                    row += 3

                    dataframe = iam_data.sel(
                        variables=[
                            v for v in variables if v in iam_data.variables.values
                        ],
                        region=region,
                        year=[y for y in iam_data.coords["year"].values if y <= 2100],
                    )

                    if len(dataframe) > 0:
                        dataframe = dataframe.to_dataframe("val")
                        dataframe = dataframe.unstack()["val"]
                        dataframe = dataframe.T
                        dataframe = dataframe.rename_axis(index=None)

                        data = dataframe_to_rows(dataframe)

                        counter = 0
                        for _, data_row in enumerate(data, 1):
                            if data_row != [None]:
                                for c_idx, value in enumerate(data_row, 1):
                                    worksheet.cell(
                                        row=row + counter,
                                        column=col + c_idx,
                                        value=value,
                                    )
                                    last_col_used = col + c_idx
                                counter += 1

                        values = Reference(
                            worksheet,
                            min_col=col + 2,
                            min_row=row,
                            max_col=col + c_idx,
                            max_row=row + counter - 1,
                        )
                        cats = Reference(
                            worksheet,
                            min_col=col + 1,
                            min_row=row + 1,
                            max_row=row + counter - 1,
                        )

                        if "generation" in sector:
                            chart = AreaChart(grouping="stacked")
                        elif "efficiency" in sector:
                            chart = LineChart()
                        elif "CCS" in sector:
                            chart = AreaChart()
                        elif "Transport" in sector:
                            chart = AreaChart(grouping="stacked")
                        else:
                            chart = LineChart()

                        chart.add_data(values, titles_from_data=True)
                        chart.set_categories(cats)
                        chart.title = f"{region} - {sector}"
                        chart.y_axis.title = metadata[sector]["label"]
                        chart.height = 8
                        chart.width = 16
                        chart.anchor = f"{get_column_letter(col + 2)}{row + 1}"
                        worksheet.add_chart(chart)

                        row += counter + 2

                scenario_list.append((scenario["model"], scenario["pathway"]))

    workbook.save(filename)


def generate_change_report(source, version, source_type, system_model):
    """
    Generate a change report of the scenarios from the log files.
    """

    # create an Excel workbook
    workbook = openpyxl.Workbook()
    workbook.remove(workbook.active)

    log_filepaths = [
        "premise_dac",
        "premise_electricity",
        "premise_fuel",
        "premise_transport",
        "premise_steel",
        "premise_cement",
        "premise_emissions",
        "premise_external_scenarios",
    ]

    # fetch YAML file containing the reporting metadata
    with open(LOG_REPORTING_FILEPATH, "r") as f:
        metadata = yaml.load(f, Loader=yaml.FullLoader)

    # create a first tab
    # where is displayed
    # the name and version of the library
    # the date of the report
    # and the name of the source database
    worksheet = workbook.create_sheet("Change report")
    worksheet.cell(row=1, column=1, value="Library name")
    worksheet.cell(row=1, column=2, value="Library version")
    worksheet.cell(row=1, column=3, value="Report date")
    worksheet.cell(row=1, column=4, value="Source database")
    worksheet.cell(row=1, column=5, value="Source database format")
    worksheet.cell(row=1, column=6, value="Database version")
    worksheet.cell(row=1, column=7, value="Database system model")
    worksheet.cell(row=2, column=1, value="premise")
    worksheet.cell(row=2, column=2, value=".".join(map(str, __version__)))
    worksheet.cell(row=2, column=3, value=datetime.now())
    worksheet.cell(row=2, column=4, value=source)
    worksheet.cell(row=2, column=5, value=source_type)
    worksheet.cell(row=2, column=6, value=version)
    worksheet.cell(row=2, column=7, value=system_model)

    dim_holder = DimensionHolder(worksheet=worksheet)
    for col in range(worksheet.min_column, worksheet.max_column + 1):
        dim_holder[get_column_letter(col)] = ColumnDimension(
            worksheet, min=col, max=col, width=20
        )

    worksheet.column_dimensions = dim_holder

    for filepath in log_filepaths:
        fp = Path(DIR_LOGS / filepath).with_suffix(".log")
        # check if log file exists
        if not fp.is_file():
            continue
        # if file exists, check that it is not empty
        if os.stat(fp).st_size == 0:
            continue

        df = convert_log_to_excel_file(fp)

        # create a worksheet for this sector
        worksheet = workbook.create_sheet(fetch_tab_name(filepath))

        # add a description of each column
        # in each row
        for col, column in enumerate(fetch_columns(fp), 1):
            worksheet.cell(
                row=1,
                column=col,
                value=metadata[filepath]["columns"][column]["description"],
            )
            worksheet.cell(
                row=2,
                column=col,
                value=metadata[filepath]["columns"][column].get("unit"),
            )

        # add the df dataframe to the sheet
        for r in dataframe_to_rows(df, index=False):
            worksheet.append(r)

    # save the workbook in the working directory
    # the file name is change_report with the current date
    fp = Path(
        DIR_LOG_REPORT / f"change_report {datetime.now().strftime('%Y-%m-%d')}.xlsx"
    )
    workbook.save(fp)


def fetch_columns(variable):
    """
    Read reporting.yaml which return
    the columns for the variable.
    """

    with open(LOG_REPORTING_FILEPATH, "r", encoding="utf-8") as stream:
        reporting = yaml.safe_load(stream)

    return list(reporting[variable.stem]["columns"].keys())


def fetch_tab_name(variable):
    """
    Read reporting.yaml which return
    the tab name for the variable.
    """

    with open(LOG_REPORTING_FILEPATH, "r", encoding="utf-8") as stream:
        reporting = yaml.safe_load(stream)

    return reporting[variable]["tab"]


def convert_log_to_excel_file(filepath):
    """
    Read the log file premise.log in the working directory.
    Load into pandas dataframe and group the data by
    scenario and variable.
    """

    try:
        df = pd.read_csv(filepath, sep="|", header=None)
        df.columns = fetch_columns(filepath)
        return df

    except EmptyDataError:
        # return an empty dataframe
        return pd.DataFrame(columns=fetch_columns(filepath))
