Mapping
=======

Link to a new IAM model
-----------------------

Although *premise* comes with a set of scenarios from the REMIND
and IMAGE IAM models, it is possible to link it to a new IAM model.
To do so, you need to populate the .yaml mapping files under the
folder https://github.com/polca/premise/tree/master/premise/iam_variables_mapping.

For each variable in each of the .yaml files, specify the
corresponding IAM variable name as follows:

.. code-block:: yaml

    Biomass CHP:
      iam_aliases:
        remind: SE|Electricity|Biomass|++|Combined Heat and Power w/o CC
        image: Secondary Energy|Electricity|Biomass|w/o CCS|3

        new_IAM: new_IAM_variable_name <--- this is the new IAM variable name

      eff_aliases:
        remind: Tech|Electricity|Biomass|Combined Heat and Power w/o CC|Efficiency
        image: Efficiency|Electricity|Biomass|w/o CCS|3

        new_IAM: new_IAM_efficiency_variable_name <--- this is the new IAM variable name

      ecoinvent_aliases:
        fltr:
          - heat and power co-generation, wood chips
        mask:
          reference product: heat
      ecoinvent_fuel_aliases:
        fltr:
          - market for wood chips, wet, measured as dry mass

If efficiency-related variables are not available, the corresponding
technologies will simply not have their efficiency adjusted.

Additionally, add your model name to the models list as well as
the list of geographical regions as LIST_xxx_REGIONS, with xxx
being the IAM model name, in the file iam_variables_mapping/constants.yaml.

Lastly, inform premise about the geographical definitions of
the IAM model you are using.
Create a .json file listing ISO 3166-1 alpha-2 country codes
and their corresponding IAM regions, as shown below, and store it under
premise/iam_variables_mapping/topologies, under the name: iamname-topology.json.

.. code-block:: json
    {
        ...
        "REF": ["AM", "AZ", "BY", "GE", "KZ", "KG", "MD", "RU", "TJ", "TM", "UA", "UZ"],
        "CAZ": ["AU", "CA", "NZ"],
        "CHA": ["CN", "HK", "MO", "TW"],
        "IND": ["IN"],
        "JPN": ["JP"],
        "USA": ["US", "PM"],
        "World": ["GLO", "RoW"]
    }

Note that the IAM region names must be identical to the ones used in the IAM scenario files.

IAM scenario file
-----------------

The scenario file should be a comma-separated text file (i.e., csv)
with data presented in a tabular format, such as:

+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| Model  | Scenario    | Region | Variable         | Unit      | 2005        | 2010        | 2015        | 2020        | 2025        |
+========+=============+========+==================+===========+=============+=============+=============+=============+=============+
| REMIND | SSP2EU-Base | CAZ    | Emi|CO2|+|Energy | Mt CO2/yr | 1011.34074  | 976.7202877 | 993.8525168 | 957.3199102 | 945.014101  |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | CHA    | Emi|CO2|+|Energy | Mt CO2/yr | 6720.313463 | 8601.575671 | 10086.37126 | 11281.46999 | 10996.79931 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | EUR    | Emi|CO2|+|Energy | Mt CO2/yr | 4235.648974 | 3730.532814 | 3392.421123 | 3114.284044 | 2860.549231 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | IND    | Emi|CO2|+|Energy | Mt CO2/yr | 1215.466496 | 1664.185158 | 2146.940653 | 2477.459967 | 2946.357462 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | JPN    | Emi|CO2|+|Energy | Mt CO2/yr | 1457.252288 | 1415.666384 | 1345.278014 | 1181.679212 | 1060.684659 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | LAM    | Emi|CO2|+|Energy | Mt CO2/yr | 1410.609298 | 1575.558465 | 1682.930038 | 1613.4512   | 1739.260156 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | MEA    | Emi|CO2|+|Energy | Mt CO2/yr | 1782.408233 | 2254.050107 | 2607.952516 | 2793.972343 | 3064.426497 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | NEU    | Emi|CO2|+|Energy | Mt CO2/yr | 378.1710003 | 421.2277231 | 477.6241091 | 498.465216  | 500.4845903 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | OAS    | Emi|CO2|+|Energy | Mt CO2/yr | 1787.07182  | 2073.863804 | 2442.52372  | 2780.880819 | 3264.746917 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+
| REMIND | SSP2EU-Base | REF    | Emi|CO2|+|Energy | Mt CO2/yr | 2551.110779 | 2472.637216 | 2544.690495 | 2607.286302 | 2681.647657 |
+--------+-------------+--------+------------------+-----------+-------------+-------------+-------------+-------------+-------------+


The following columns must be present:

* Region
* Variable
* Unit

as well as the time steps (e..g, 2005 to 2100).
Other columns can be present, but they will be ignored.

You need to point to that file when initiating `NewDatabase`, like so:

.. code-block:: python

    ndb = NewDatabase(
        scenarios = [{"model":"remind", "pathway":"my_special_scenario", "year":2028,
                      "filepath":r"C:\filepath\to\your\scenario\folder"}],
        source_db="ecoinvent 3.6 cutoff", # <-- name of the database
        source_version="3.6", # <-- version of ecoinvent
    )

There are essentially two types of variables needed from the IAM scenario files:

- variables that relate to the production volumes of technologies. These variables are used to scale the production volumes of the corresponding activities in the ecoinvent database. For example, if the IAM scenario file contains a variable named ``Electricity|Production|Wind`` for the region ``EUR``, it will help premise calculate the share of wind power in the electricity consumption mix of the said region. Hence, the unit of such variables should refer to a production volume over time (e.g., ``GWh/year``, ``EJ/year``, etc.).
- variables that relate to the efficiency of technologies over time. These variables are used to calculate scaling factors (which are relative by default to 2020), to adjust the energy or material efficiency of the corresponding activities in the ecoinvent database. For example, if the IAM scenario file contains a variable named ``Electricity|Efficiency|Coal`` for the region ``EUR``, it will help premise adjust the amount of coal and related emissions per unit of kWh produced in the said region. Hence, the unit of such variables can be unitless, or relate to an efficiency ratio or percentage.