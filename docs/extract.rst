EXTRACT
=======

The **EXTRACT** phase consists of the following steps:

* Extraction and cleaning of the ecoinvent database
* Import and cleaning of additional inventories
* Import and cleaning of user-provided inventories (optional)
* Caching, if these database and inventories are imported for the first time
* Loading of IAM data

Current IAM scenarios
"""""""""""""""""""""

*premise* includes several Integrated Assessment Model (IAM) scenarios,
but you can also use other scenarios.
In *premise*, scenarios are defined by their Shared Socio-economic
Pathway (SSP), a climate trajectory—often represented by a Representative
Concentration Pathway (RCP)—and a year (e.g., SSP1, Base, 2035).


+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP/RCP scenario | GMST increase by 2100 | Society/economy trend                                                              | Climate policy                              | REMIND          | IMAGE      |
+==================+=======================+====================================================================================+=============================================+=================+============+
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP1-None        | 2.3-2.8 °C            | Optimistic trends for human develop. and economy, driven by sustainable practices. | None                                        | SSP1-Base       | SSP1-Base  |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP1-None        | ~2.2 °C               | Optimistic trends for human develop. and economy, driven by sustainable practices. | National Policies Implemented (NPI).        | SSP1-NPi        |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP1-None        | ~1.9 °C               | Optimistic trends for human develop. and economy, driven by sustainable practices. | Nationally Determined Contributions (NDCs). | SSP1-NDC        |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP1-RCP2.6      | ~1.7 °C               | Optimistic trends for human develop. and economy, driven by sustainable practices. | Paris Agreement objective.                  | SSP1-PkBudg1150 |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP1-RCP1.9      | ~1.3 °C               | Optimistic trends for human develop. and economy, driven by sustainable practices. | Paris Agreement objective.                  | SSP1-PkBudg500  |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP2-None        | ~3.5 °C               | Extrapolation from historical developments.                                        | None (eq. to RCP6)                          | SSP2-Base       | SSP2-Base  |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP2-None        | ~3.3 °C               | Extrapolation from historical developments.                                        | National Policies Implemented (NPI).        | SSP2-NPi        |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP2-None        | ~2.5 °C               | Extrapolation from historical developments.                                        | Nationally Determined Contributions (NDCs). | SSP2-NDC        |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP2-RCP2.6      | 1.6-1.8 °C            | Extrapolation from historical developments.                                        | Paris Agreement objective.                  | SSP2-PkBudg1150 | SSP2-RCP26 |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP2-RCP1.9      | 1.2-1.4 °C            | Extrapolation from historical developments.                                        | Paris Agreement objective.                  | SSP2-PkBudg500  | SSP2-RCP19 |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP5-None        | ~4.5 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | None                                        | SSP5-Base       |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP5-None        | ~4.0 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | National Policies Implemented (NPI).        | SSP5-NPi        |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP5-None        | ~3.0 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | Nationally Determined Contributions (NDCs). | SSP5-NDC        |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP5-RCP2.6      | ~1.7 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | Paris Agreement objective.                  | SSP5-PkBudg1150 |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+
| SSP5-RCP1.9      | ~1.0 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | Paris Agreement objective.                  | SSP5-PkBudg500  |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+


.. note::

    A summary report of the main variables of the scenarios
    selected is generated automatically after each database export.
    You can also generate it manually:

.. python::

    ndb = NewDatabase(...)
    ndb.generate_scenario_report()


Supported versions of ecoinvent
"""""""""""""""""""""""""""""""

*premise* currently works with the following ecoinvent database versions:

* v.3.5, cut-off
* v.3.6, cut-off
* v.3.7, cut-off
* v.3.7.1, cut-off
* **v.3.8, cut-off and consequential**
* **v.3.9/3.9.1, cut-off and consequential**


Supported sources of ecoinvent
""""""""""""""""""""""""""""""

*premise* can extract the ecoinvent database from:

* a brightway2_ project that contains the ecoinvent database
* ecosposld2 files, that can be downloaded from the ecoinvent_ website

.. _ecoinvent: https://ecoinvent.org
.. _brightway2: https://brightway.dev/



From a brightway2 project
-------------------------

To extract from an ecoinvent database located in a brightway2 project, simply
indicate the database name in `source_db` and its version in `source_version`:

.. code-block:: python

  from premise import *
  import brightway2 as bw

  bw.projects.set_current("my_project)

  ndb = NewDatabase(
        scenarios=[
                {"model":"remind", "pathway":"SSP2-Base", "year":2028}
            ],
        source_db="ecoinvent 3.7 cutoff", # <-- this is NEW.
        source_version="3.7.1", # <-- this is NEW
        key='xxxxxxxxxxxxxxxxxxxxxxxxx'
    )

Note that a cache of the database will be created the first time and
store in the library folder. Any subsequent creation of databases
using the same ecoinvent version will no longer require this extraction
step.

If you wish to clear that cache folder, do:

.. code-block:: python

    from premise import *

    clear_cache()

.. note::

    It is recommended to restart your notebook once
    the data has been cached for the first time, so that
    the remaining steps can be performed using the
    cached data (much faster).


From ecospold2 files
--------------------

To extract from a set of ecospold2 files, you need to point to the location of those files
in `source_file_path`, as well as indicate the database format in `source_type`:

.. code-block:: python

    from premise import *

    ndb = NewDatabase(
        scenarios = [
            {"model":"remind", "pathway":"SSP2-Base", "year":2028}
                    ],
        source_type="ecospold", # <--- this is NEW
        source_file_path=r"C:\file\path\to\ecoinvent 3.5_cutoff_ecoSpold02\datasets", # <-- this is NEW
        source_version="3.5",
    )

Import of additional inventories
""""""""""""""""""""""""""""""""

After the ecoinvent database is extracted and checked, a number of additional inventories
are imported, regardless of the year of scenario that is being considered.


Power generation
----------------

A number of  datasets relating to power generation not originally present in
ecoinvent are imported. The next sub-sections lists such datasets.

Power plants with CCS
*********************

Datasets for power generation with Carbon Capture and Storage (CCS) are imported.
They originate from Volkart_ et al. 2013, and can be consulted here: LCI_Power_generation_.
An exception to this are the inventories for biomass-based integrated gasification combined cycle power plants (BIGCCS),
which are from Briones-Hidrovo_ et al, 2020.

.. _Volkart: https://doi.org/10.1016/j.ijggc.2013.03.003
.. _Briones-Hidrovo: https://doi.org/10.1016/j.jclepro.2020.125680
.. _LCI_Power_generation: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-Carma-CCS.xlsx

The table below lists the names of the new activities (only production datasets are shown).

 ============================================================================================================= ===========
  Power generation with CCS (activities list)                                                                   location
 ============================================================================================================= ===========
  electricity production, at power plant/hard coal, IGCC, no CCS                                                RER
  electricity production, at power plant/hard coal, PC, no CCS                                                  RER
  electricity production, at power plant/hard coal, oxy, pipeline 200km, storage 1000m                          RER
  electricity production, at power plant/hard coal, oxy, pipeline 400km, storage 3000m                          RER
  electricity production, at power plant/hard coal, post, pipeline 200km, storage 1000m                         RER
  electricity production, at power plant/hard coal, post, pipeline 400km, storage 1000m                         RER
  electricity production, at power plant/hard coal, post, pipeline 400km, storage 3000m                         RER
  electricity production, at power plant/hard coal, pre, pipeline 200km, storage 1000m                          RER
  electricity production, at power plant/hard coal, pre, pipeline 400km, storage 3000m                          RER
  electricity production, at power plant/lignite, IGCC, no CCS                                                  RER
  electricity production, at power plant/lignite, PC, no CCS                                                    RER
  electricity production, at power plant/lignite, oxy, pipeline 200km, storage 1000m                            RER
  electricity production, at power plant/lignite, oxy, pipeline 400km, storage 3000m                            RER
  electricity production, at power plant/lignite, post, pipeline 200km, storage 1000m                           RER
  electricity production, at power plant/lignite, post, pipeline 400km, storage 3000m                           RER
  electricity production, at power plant/lignite, pre, pipeline 200km, storage 1000m                            RER
  electricity production, at power plant/lignite, pre, pipeline 400km, storage 3000m                            RER
  electricity production, at power plant/natural gas, ATR H2-CC, no CCS                                         RER
  electricity production, at power plant/natural gas, NGCC, no CCS/kWh                                          RER
  electricity production, at power plant/natural gas, post, pipeline 200km, storage 1000m                       RER
  electricity production, at power plant/natural gas, post, pipeline 400km, storage 1000m                       RER
  electricity production, at power plant/natural gas, post, pipeline 400km, storage 3000m                       RER
  electricity production, at power plant/natural gas, pre, pipeline 200km, storage 1000m                        RER
  electricity production, at power plant/natural gas, pre, pipeline 400km, storage 3000m                        RER
  electricity production, at wood burning power plant 20 MW, truck 25km, no CCS                                 RER
  electricity production, at wood burning power plant 20 MW, truck 25km, post, pipeline 200km, storage 1000m    RER
  electricity production, at wood burning power plant 20 MW, truck 25km, post, pipeline 400km, storage 3000m    RER
 ============================================================================================================= ===========


Natural gas
***********

Updated inventories relating to natural gas extraction and distribution
are imported to substitute some of the original ecoinvent dataset.
These datasets originate from ESU Services and come with a report_,
and can be consulted here: LCI_Oil_NG_.

.. _LCI_Oil_NG: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-ESU-oil-and-gas.xlsx

They have been adapted to a brightway2-compatible format.
These new inventories have, among other things, higher methane slip
emissions along the natural gas supply chain, especially at extraction.

.. _report: http://www.esu-services.ch/fileadmin/download/publicLCI/meili-2021-LCI%20for%20the%20oil%20and%20gas%20extraction.pdf


 ========================================================== ==============================================================
  Original dataset                                           Replaced by
 ========================================================== ==============================================================
  natural gas production (natural gas, high pressure), DE    natural gas, at production (natural gas, high pressure), DE
  natural gas production (natural gas, high pressure), DZ    natural gas, at production (natural gas, high pressure), DZ
  natural gas production (natural gas, high pressure), US    natural gas, at production (natural gas, high pressure), US
  natural gas production (natural gas, high pressure), RU    natural gas, at production (natural gas, high pressure), RU
  petroleum and gas production, GB                           natural gas, at production (natural gas, high pressure), GB
  petroleum and gas production, NG                           natural gas, at production (natural gas, high pressure), NG
  petroleum and gas production, NL                           natural gas, at production (natural gas, high pressure), NL
  petroleum and gas production, NO                           natural gas, at production (natural gas, high pressure), NO
 ========================================================== ==============================================================

The original natural gas datasets are preserved, but they do not provide input to any
other datasets in the database. The new datasets provide natural gas at high pressure to
the original supply chains, which remain unchanged.

The table below lists the names of the new activities (only high pressure datasets are shown).

 ============================= ===========
  Natural gas extraction        location
 ============================= ===========
  natural gas, at production    AZ
  natural gas, at production    RO
  natural gas, at production    LY
  natural gas, at production    SA
  natural gas, at production    IQ
  natural gas, at production    RU
  natural gas, at production    NL
  natural gas, at production    DZ
  natural gas, at production    NG
  natural gas, at production    DE
  natural gas, at production    KZ
  natural gas, at production    NO
  natural gas, at production    QA
  natural gas, at production    GB
  natural gas, at production    MX
  natural gas, at production    US
 ============================= ===========

.. note::

    This import will be removed in the future, as the original
    ecoinvent dataset will be updated (i.e., v3.9).

Photovoltaic panels
*******************

Photovoltaic panel inventories originate the IEA's Task 12 project IEA_PV_. They have been adapted
into a brightway2-friendly format. They can be consulted here: LCI_PV_.

.. _IEA_PV: https://iea-pvps.org/wp-content/uploads/2020/12/IEA-PVPS-LCI-report-2020.pdf
.. _LCI_PV: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-PV.xlsx

They consist of the following PV installation types:

 ============================================================================================ ===========
  PV installation                                                                              location
 ============================================================================================ ===========
  photovoltaic slanted-roof installation, 1.3 MWp, multi-Si, panel, mounted, on roof           CH
  photovoltaic flat-roof installation, 156 kWp, multi-Si, on roof                              CH
  photovoltaic flat-roof installation, 156 kWp, single-Si, on roof                             CH
  photovoltaic flat-roof installation, 280 kWp, multi-Si, on roof                              CH
  photovoltaic flat-roof installation, 280 kWp, single-Si, on roof                             CH
  photovoltaic flat-roof installation, 324 kWp, multi-Si, on roof                              DE
  photovoltaic slanted-roof installation, 3 kWp, CIS, laminated, integrated, on roof           CH
  photovoltaic slanted-roof installation, 3 kWp, CIS, laminated, integrated, on roof           RER
  photovoltaic slanted-roof installation, 3 kWp, CdTe, panel, mounted, on roof                 CH
  photovoltaic slanted-roof installation, 3 kWp, CdTe, panel, mounted, on roof                 RER
  photovoltaic slanted-roof installation, 3 kWp, micro-Si, laminated, integrated, on roof      RER
  photovoltaic slanted-roof installation, 3 kWp, micro-Si, panel, mounted, on roof             RER
  photovoltaic flat-roof installation, 450 kWp, single-Si, on roof                             DE
  photovoltaic open ground installation, 560 kWp, single-Si, on open ground                    CH
  photovoltaic open ground installation, 569 kWp, multi-Si, on open ground                     ES
  photovoltaic open ground installation, 570 kWp, CIS, on open ground                          RER
  photovoltaic open ground installation, 570 kWp, CdTe, on open ground                         RER
  photovoltaic open ground installation, 570 kWp, micro-Si, on open ground                     RER
  photovoltaic open ground installation, 570 kWp, multi-Si, on open ground                     ES
  photovoltaic open ground installation, 570 kWp, multi-Si, on open ground                     RER
  photovoltaic open ground installation, 570 kWp, single-Si, on open ground                    RER
  photovoltaic slanted-roof installation, 93 kWp, multi-Si, laminated, integrated, on roof     CH
  photovoltaic slanted-roof installation, 93 kWp, multi-Si, panel, mounted, on roof            CH
  photovoltaic slanted-roof installation, 93 kWp, single-Si, laminated, integrated, on roof    CH
  photovoltaic slanted-roof installation, 93 kWp, single-Si, panel, mounted, on roof           CH
 ============================================================================================ ===========


Although these datasets have a limited number of locations (CH, RER, DE, ES),
the IEA report provides country-specific load factors:

 ======================= =========== ========= ==========
  production [kWh/kWp]    roof-top    façade    central
 ======================= =========== ========= ==========
  PT                      1427        999       1513
  IL                      1695        1187      1798
  SE                      919         643       974
  FR                      968         678       1026
  TR                      1388        971       1471
  NZ                      1240        868       1315
  MY                      1332        933       1413
  CN                      971         679       1029
  TH                      1436        1005      1522
  ZA                      1634        1144      1733
  JP                      1024        717       1086
  CH                      976         683       1040
  DE                      922         645       978
  KR                      1129        790       1197
  AT                      1044        731       1111
  GR                      1323        926       1402
  IE                      796         557       844
  AU                      1240        868       1314
  IT                      1298        908       1376
  MX                      1612        1128      1709
  NL                      937         656       994
  GB                      848         593       899
  ES                      1423        996       1509
  CL                      1603        1122      1699
  HU                      1090        763       1156
  CZ                      944         661       1101
  CA                      1173        821       1243
  US                      1401        981       1485
  NO                      832         583       882
  FI                      891         624       945
  BE                      908         635       962
  DK                      971         680       1030
  LU                      908         635       962
 ======================= =========== ========= ==========


In the report, the generation potential per installation type is multiplied by the number of installations
in each country, to produce country-specific PV power mix datasets normalized to 1 kWh.
The report specifies the production-weighted PV mix for each country, but we further split it
between residential (<=3kWp) and commercial (>3kWp) installations
(as most IAMs make such distinction):

 ==================================================== ===========
  Production-weighted PV mix                           location
 ==================================================== ===========
  electricity production, photovoltaic, residential    PT
  electricity production, photovoltaic, residential    IL
  electricity production, photovoltaic, residential    SE
  electricity production, photovoltaic, residential    FR
  electricity production, photovoltaic, residential    TR
  electricity production, photovoltaic, residential    NZ
  electricity production, photovoltaic, residential    MY
  electricity production, photovoltaic, residential    CN
  electricity production, photovoltaic, residential    TH
  electricity production, photovoltaic, residential    ZA
  electricity production, photovoltaic, residential    JP
  electricity production, photovoltaic, residential    CH
  electricity production, photovoltaic, residential    DE
  electricity production, photovoltaic, residential    KR
  electricity production, photovoltaic, residential    AT
  electricity production, photovoltaic, residential    GR
  electricity production, photovoltaic, residential    IE
  electricity production, photovoltaic, residential    AU
  electricity production, photovoltaic, residential    IT
  electricity production, photovoltaic, residential    MX
  electricity production, photovoltaic, residential    NL
  electricity production, photovoltaic, residential    GB
  electricity production, photovoltaic, residential    ES
  electricity production, photovoltaic, residential    CL
  electricity production, photovoltaic, residential    HU
  electricity production, photovoltaic, residential    CZ
  electricity production, photovoltaic, residential    CA
  electricity production, photovoltaic, residential    US
  electricity production, photovoltaic, residential    NO
  electricity production, photovoltaic, residential    FI
  electricity production, photovoltaic, residential    BE
  electricity production, photovoltaic, residential    DK
  electricity production, photovoltaic, residential    LU
  electricity production, photovoltaic, commercial     PT
  electricity production, photovoltaic, commercial     IL
  electricity production, photovoltaic, commercial     SE
  electricity production, photovoltaic, commercial     FR
  electricity production, photovoltaic, commercial     TR
  electricity production, photovoltaic, commercial     NZ
  electricity production, photovoltaic, commercial     MY
  electricity production, photovoltaic, commercial     CN
  electricity production, photovoltaic, commercial     TH
  electricity production, photovoltaic, commercial     ZA
  electricity production, photovoltaic, commercial     JP
  electricity production, photovoltaic, commercial     CH
  electricity production, photovoltaic, commercial     DE
  electricity production, photovoltaic, commercial     KR
  electricity production, photovoltaic, commercial     AT
  electricity production, photovoltaic, commercial     GR
  electricity production, photovoltaic, commercial     IE
  electricity production, photovoltaic, commercial     AU
  electricity production, photovoltaic, commercial     IT
  electricity production, photovoltaic, commercial     MX
  electricity production, photovoltaic, commercial     NL
  electricity production, photovoltaic, commercial     GB
  electricity production, photovoltaic, commercial     ES
  electricity production, photovoltaic, commercial     CL
  electricity production, photovoltaic, commercial     HU
  electricity production, photovoltaic, commercial     CZ
  electricity production, photovoltaic, commercial     CA
  electricity production, photovoltaic, commercial     US
  electricity production, photovoltaic, commercial     NO
  electricity production, photovoltaic, commercial     FI
  electricity production, photovoltaic, commercial     BE
  electricity production, photovoltaic, commercial     DK
  electricity production, photovoltaic, commercial     LU
 ==================================================== ===========

Hence, inside the *residential* PV mix of Spain ("electricity production, photovoltaic, residential"),
one will find the following inputs for the production of 1kWh:

 ========================================================================================== ============== =========== ============
  name                                                                                       amount         location    unit
 ========================================================================================== ============== =========== ============
  Energy, solar, converted                                                                   3.8503                     megajoule
  Heat, waste                                                                                0.25027                    megajoule
  photovoltaic slanted-roof installation, 3 kWp, CIS, laminated, integrated, on roof         2.48441E-08    CH          unit
  photovoltaic slanted-roof installation, 3 kWp, CdTe, panel, mounted, on roof               4.99911E-07    CH          unit
  photovoltaic slanted-roof installation, 3 kWp, micro-Si, laminated, integrated, on roof    3.93869E-09    RER         unit
  photovoltaic slanted-roof installation, 3 kWp, micro-Si, panel, mounted, on roof           6.55186E-08    RER         unit
  photovoltaic facade installation, 3kWp, multi-Si, laminated, integrated, at building       2.10481E-07    RER         unit
  photovoltaic facade installation, 3kWp, multi-Si, panel, mounted, at building              2.10481E-07    RER         unit
  photovoltaic facade installation, 3kWp, single-Si, laminated, integrated, at building      1.11463E-07    RER         unit
  photovoltaic facade installation, 3kWp, single-Si, panel, mounted, at building             1.11463E-07    RER         unit
  photovoltaic flat-roof installation, 3kWp, multi-Si, on roof                               2.20794E-06    RER         unit
  photovoltaic flat-roof installation, 3kWp, single-Si, on roof                              1.17025E-06    RER         unit
  photovoltaic slanted-roof installation, 3kWp, CIS, panel, mounted, on roof                 4.12805E-07    CH          unit
  photovoltaic slanted-roof installation, 3kWp, CdTe, laminated, integrated, on roof         3.00704E-08    CH          unit
  photovoltaic slanted-roof installation, 3kWp, multi-Si, laminated, integrated, on roof     1.08693E-07    RER         unit
  photovoltaic slanted-roof installation, 3kWp, multi-Si, panel, mounted, on roof            1.81407E-06    RER         unit
  photovoltaic slanted-roof installation, 3kWp, single-Si, laminated, integrated, on roof    5.75655E-08    RER         unit
  photovoltaic slanted-roof installation, 3kWp, single-Si, panel, mounted, on roof           9.6195E-07     RER         unit
 ========================================================================================== ============== =========== ============

with, for example, 2.48E-8 units of "photovoltaic slanted-roof installation, 3 kWp, CIS, laminated, integrated, on roof"
being calculated as:

.. code-block::

    1 / (30 [years] * 1423 [kWh/kWp] * 0.32% [share of PV capacity of such type installed in Spain])

Note that commercial PV mix datasets provide electricity at high voltage, unlike residential
PV mix datasets, which supply at low voltage only.

Geothermal
**********

Heat production by means of a geothermal well are not represented in ecoinvent.
The geothermal power plant construction inventories are from Maeder_ Bachelor Thesis.

.. _Maeder: https://www.psi.ch/sites/default/files/import/ta/PublicationTab/BSc_Mattia_Maeder_2016.pdf

The co-generation unit has been removed and replaced by heat exchanger and
district heating pipes. Gross heat output of 1,483 TJ, with 80% efficiency.

The inventories can be consulted here: LCIgeothermal_.

.. _LCIgeothermal: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-geothermal.xlsx

They introduce the following datasets (only heat production datasets shown):

 =================================== ===========
  Geothermal heat production          location
 =================================== ===========
  heat production, deep geothermal    RAS
  heat production, deep geothermal    GLO
  heat production, deep geothermal    RAF
  heat production, deep geothermal    RME
  heat production, deep geothermal    RLA
  heat production, deep geothermal    RU
  heat production, deep geothermal    CA
  heat production, deep geothermal    JP
  heat production, deep geothermal    US
  heat production, deep geothermal    IN
  heat production, deep geothermal    CN
  heat production, deep geothermal    RER
 =================================== ===========


Hydrogen
--------

*premise* imports inventories for hydrogen production via the following pathways:

* Steam Methane Reforming, using natural gas
* Steam Methane Reforming, using natural gas, with Carbon Capture and Storage
* Steam Methane Reforming, using bio-methane
* Steam Methane Reforming, using bio-methane, with Carbon Capture and Storage
* Auto Thermal Reforming, using natural gas
* Auto Thermal Reforming, using natural gas, with Carbon Capture and Storage
* Auto Thermal Reforming, using bio-methane
* Auto Thermal Reforming, using bio-methane, with Carbon Capture and Storage
* Woody biomass gasification, using a fluidized bed
* Woody biomass gasification, using a fluidized bed, with Carbon Capture and Storage
* Woody biomass gasification, using an entrained flow gasifier
* Woody biomass gasification, using an entrained flow gasifier, with Carbon Capture and Storage
* Coal gasification
* Coal gasification, with Carbon Capture and Storage
* Electrolysis
* Thermochemical water splitting
* Pyrolysis

Inventories using Steam Methane Reforming are from Antonini_ et al. 2021.
They can be consulted here: LCI_SMR_.
Inventories using Auto Thermal Reforming are from Antonini_ et al. 2021.
They can be consulted here: LCI_ATR_.
Inventories using Woody biomass gasification are from Antonini2_ et al. 2021.
They can be consulted here: LCI_woody_.
Inventories using coal gasification are from Wokaun_ et al. 2015, but updated
with Li_ et al. 2022, which also provide an option with CCS.
They can be consulted here: LCI_coal_.
Inventories using electrolysis are from Niklas Gerloff_. 2021.
They can be consulted here: LCI_electrolysis_.
Inventories for thermochemical water splitting are from Zhang2_ et al. 2022.
Inventories for pyrolysis are from Al-Qahtani_ et al. 2021, completed with
data from Postels_ et al., 2016.

.. _Antonini: https://pubs.rsc.org/en/content/articlelanding/2020/se/d0se00222d
.. _Antonini2: https://pubs.rsc.org/en/Content/ArticleLanding/2021/SE/D0SE01637C
.. _Wokaun: https://www.cambridge.org/core/books/transition-to-hydrogen/43144AF26ED80E7106B675A6E83B1579
.. _Li: https://doi.org/10.1016/j.jclepro.2022.132514
.. _Gerloff: https://doi.org/10.1016/j.est.2021.102759
.. _Zhang2: https://doi.org/10.1016/j.ijhydene.2022.02.150
.. _Al-Qahtani: https://doi.org/10.1016/j.apenergy.2020.115958
.. _Postels: https://doi.org/10.1016/j.ijhydene.2016.09.167
.. _LCI_SMR: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-hydrogen-smr-atr-natgas.xlsx
.. _LCI_ATR: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-hydrogen-smr-atr-natgas.xlsx
.. _LCI_woody: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-hydrogen-wood-gasification.xlsx
.. _LCI_coal: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-hydrogen-coal-gasification.xlsx
.. _LCI_electrolysis: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-hydrogen-electrolysis.xlsx

The new datasets introduced are listed in the table below (only production datasets are shown).

 ======================================================================================================================================= ===========
  Hydrogen production                                                                                                                     location
 ======================================================================================================================================= ===========
  hydrogen production, steam methane reforming of natural gas, 25 bar                                                                     CH
  hydrogen production, steam methane reforming of natural gas, with CCS (MDEA, 98% eff.), 25 bar                                          CH
  hydrogen production, steam methane reforming, from biomethane, high and low temperature, with CCS (MDEA, 98% eff.), 26 bar              CH
  hydrogen production, steam methane reforming, from biomethane, high and low temperature, 26 bar                                         CH
  hydrogen production, auto-thermal reforming, from biomethane, 25 bar                                                                    CH
  hydrogen production, auto-thermal reforming, from biomethane, with CCS (MDEA, 98% eff.), 25 bar                                         CH
  hydrogen production, gaseous, 25 bar, from heatpipe reformer gasification of woody biomass with CCS, at gasification plant              CH
  hydrogen production, gaseous, 25 bar, from heatpipe reformer gasification of woody biomass, at gasification plant                       CH
  hydrogen production, gaseous, 25 bar, from gasification of woody biomass in entrained flow gasifier, with CCS, at gasification plant    CH
  hydrogen production, gaseous, 25 bar, from gasification of woody biomass in entrained flow gasifier, at gasification plant              CH
  hydrogen production, gaseous, 30 bar, from hard coal gasification and reforming, at coal gasification plant                             RER
  hydrogen production, gaseous, 30 bar, from PEM electrolysis, from grid electricity                                                      RER
  hydrogen production, gaseous, 20 bar, from AEC electrolysis, from grid electricity                                                      RER
  hydrogen production, gaseous, 1 bar, from SOEC electrolysis, from grid electricity                                                      RER
  hydrogen production, gaseous, 1 bar, from SOEC electrolysis, with steam input, from grid electricity                                    RER
  hydrogen production, gaseous, 25 bar, from thermochemical water splitting, at solar tower                                               RER
  hydrogen production, gaseous, 100 bar, from methane pyrolysis                                                                           RER
 ======================================================================================================================================= ===========

Hydrogen storage and distribution
*********************************

A number of datasets relating to hydrogen storage and distribution are also imported.

They are necessary to model the distribution of hydrogen:

* via re-assigned transmission and distribution CNG pipelines, in a gaseous state
* via dedicated transmission and distribution hydrogen pipelines, in a gaseous state
* as a liquid organic compound, by hydrogenation
* via truck, in a liquid state
* hydrogen refuelling station


Small and large storage solutions are also provided:
* high pressure hydrogen storage tank
* geological storage tank

These datasets originate from the work of Wulf_ et al. 2018, and can be
consulted here: LCI_H2_distr_. For re-assigned CNG pipelines, which require the hydrogen
to be mixed together with oxygen to limit metal embrittlement,
some parameters are taken from the work of Cerniauskas_ et al. 2020.

The datasets introduced are listed in the table below.

 ================================================================== ===========
  Hydrogen distribution                                              location
 ================================================================== ===========
  Hydrogen refuelling station                                        GLO
  high pressure hydrogen storage tank                                GLO
  distribution pipeline for hydrogen, dedicated hydrogen pipeline    RER
  transmission pipeline for hydrogen, dedicated hydrogen pipeline    RER
  zinc coating for hydrogen pipeline                                 RER
  hydrogenation of hydrogen                                          RER
  dehydrogenation of hydrogen                                        RER
  dibenzyltoluene production                                         RER
  solution mining for geological hydrogen storage                    RER
  geological hydrogen storage                                        RER
  hydrogen embrittlement inhibition                                  RER
  distribution pipeline for hydrogen, reassigned CNG pipeline        RER
  transmission pipeline for hydrogen, reassigned CNG pipeline        RER
 ================================================================== ===========


.. _Wulf: https://www.sciencedirect.com/science/article/pii/S095965261832170X
.. _LCI_H2_distr: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-hydrogen-distribution.xlsx
.. _Cerniauskas: https://doi.org/10.1016/j.ijhydene.2020.02.121

Biofuels
--------

Inventories for energy crops- and residues-based production of bioethanol and biodiesel
are imported, and can be consulted here: LCI_biofuels_. They include the farming of the crop,
the conversion of hte biomass to fuel, as well as its distribution. The conversion process
often leads to the production of co-products (dried distiller's grain, electricity, CO2, bagasse.).
Hence, energy, economic and system expansion partitioning approaches are available.
These inventories originate from several different sources
(Wu_ et al. 2006 (2020 update), Cozzolino_ 2018, Pereira_ et al. 2019 and Gonzalez-Garcia_ et al. 2012),
Cavalett_ & Cherubini 2022, as indicated in the table below.

.. _LCI_biofuels: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-biofuels.xlsx
.. _Cozzolino: https://www.psi.ch/sites/default/files/2019-09/Cozzolino_377125_%20Research%20Project%20Report.pdf
.. _Gonzalez-Garcia: https://doi.org/10.1016/j.scitotenv.2012.07.044
.. _Wu: http://greet.es.anl.gov/publication-2lli584z
.. _Pereira: http://task39.sites.olt.ubc.ca/files/2019/04/Task-39-GHS-models-Final-Report-Phase-1.pdf
.. _Cavalett: https://doi.org/10.1002/bbb.2395

The following datasets are introduced:

 ================================================================================== =========== =============================
  Activity                                                                           Location    Source
 ================================================================================== =========== =============================
  Farming and supply of switchgrass                                                  US          Wu et al. 2006 (2020 update)
  Farming and supply of poplar                                                       US          Wu et al. 2006 (2020 update)
  Farming and supply of willow                                                       US          Wu et al. 2006 (2020 update)
  Supply of forest residue                                                           US          Wu et al. 2006 (2020 update)
  Farming and supply of miscanthus                                                   US          Wu et al. 2006 (2020 update)
  Farming and supply of corn stover                                                  US          Wu et al. 2006 (2020 update)
  Farming and supply of sugarcane                                                    US          Wu et al. 2006 (2020 update)
  Farming and supply of Grain Sorghum                                                US          Wu et al. 2006 (2020 update)
  Farming and supply of Sweet Sorghum                                                US          Wu et al. 2006 (2020 update)
  Farming and supply of Forage Sorghum                                               US          Wu et al. 2006 (2020 update)
  Farming and supply of corn                                                         US          Wu et al. 2006 (2020 update)
  Farming and supply of sugarcane                                                    BR          Pereira et al. 2019/RED II
  Farming and supply of sugarcane straw                                              BR          Pereira et al. 2019
  Farming and supply of eucalyptus                                                   ES          Gonzalez-Garcia et al. 2012
  Farming and supply of wheat grains                                                 RER         Cozzolino 2018
  Farming and supply of wheat straw                                                  RER         Cozzolino 2018
  Farming and supply of corn                                                         RER         Cozzolino 2018/RED II
  Farming and supply of sugarbeet                                                    RER         Cozzolino 2018
  Supply of forest residue                                                           RER         Cozzolino 2018
  Supply and refining of waste cooking oil                                           RER         Cozzolino 2018
  Farming and supply of rapeseed                                                     RER         Cozzolino 2018/RED II
  Farming and supply of palm fresh fruit bunch                                       RER         Cozzolino 2018
  Farming and supply of dry algae                                                    RER         Cozzolino 2018
  Ethanol production, via fermentation, from switchgrass                             US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from poplar                                  US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from willow                                  US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from forest residue                          US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from miscanthus                              US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from corn stover                             US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from sugarcane                               US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from grain sorghum                           US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from sweet sorghum                           US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from forage sorghum                          US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from corn                                    US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from corn, with carbon capture               US          Wu et al. 2006 (2020 update)
  Ethanol production, via fermentation, from sugarcane                               BR          Pereira et al. 2019
  Ethanol production, via fermentation, from sugarcane straw                         BR          Pereira et al. 2019
  Ethanol production, via fermentation, from eucalyptus                              ES          Gonzalez-Garcia et al. 2012
  Ethanol production, via fermentation, from wheat grains                            RER         Cozzolino 2018
  Ethanol production, via fermentation, from wheat straw                             RER         Cozzolino 2018
  Ethanol production, via fermentation, from corn starch                             RER         Cozzolino 2018
  Ethanol production, via fermentation, from sugarbeet                               RER         Cozzolino 2018
  Ethanol production, via fermentation, from forest residue                          RER         Cozzolino 2018
  Ethanol production, via fermentation, from forest residues                         RER         Cavalett & Cherubini 2022
  Ethanol production, via fermentation, from forest product (non-residual)           RER         Cavalett & Cherubini 2022
  Biodiesel production, via transesterification, from used cooking oil               RER         Cozzolino 2018
  Biodiesel production, via transesterification, from rapeseed oil                   RER         Cozzolino 2018
  Biodiesel production, via transesterification, from palm oil, energy allocation    RER         Cozzolino 2018
  Biodiesel production, via transesterification, from algae, energy allocation       RER         Cozzolino 2018
  Biodiesel production, via Fischer-Tropsch, from forest residues                    RER         Cavalett & Cherubini 2022
  Biodiesel production, via Fischer-Tropsch, from forest product (non-residual)      RER         Cavalett & Cherubini 2022
  Kerosene production, via Fischer-Tropsch, from forest residues                     RER         Cavalett & Cherubini 2022
  Kerosene production, via Fischer-Tropsch, from forest product (non-residual)       RER         Cavalett & Cherubini 2022
 ================================================================================== =========== =============================

Synthetic fuels
---------------

*premise* imports inventories for the synthesis of hydrocarbon fuels
following two pathways:

* *Fischer-Tropsch*: it uses hydrogen and CO (from CO2 via a reverse water gas
  shift process) to produce "syncrude", which is distilled into diesel, kerosene,
  naphtha and lubricating oil and waxes. Inventories are from van der Giesen_ et al. 2014.
* *Methanol-to-liquids*: methanol is synthesized from hydrogen and CO2, and further
  distilled into gasoline, diesel, LGP and kerosene. Synthetic methanol inventories
  are from Hank_ et al. 2019. The methanol to fuel process specifications are from
  FVV_ 2013.
* *Electro-chemical methanation*: methane is produced from hydrogen and CO2 using
  a Sabatier methanation reactor. Inventories are from Zhang_ et al, 2019.

.. _Giesen: https://pubs.acs.org/doi/abs/10.1021/es500191g
.. _Hank: https://doi.org/10.1039/C9SE00658C
.. _FVV: https://www.fvv-net.de/fileadmin/user_upload/medien/materialien/FVV-Kraftstoffstudie_LBST_2013-10-30.pdf
.. _Zhang: https://doi.org/10.1039/C9SE00986H

In their default configuration, these fuels use hydrogen from electrolysis and CO2
from direct air capture (DAC). However, *premise* builds different configurations
(i.e., CO2 and hydrogen sources) for these fuels, for each IAM region:

 ============================================================================================================================================================================ ================== =============================
  Fuel production dataset                                                                                                                                                      location           source
 ============================================================================================================================================================================ ================== =============================
  Diesel production, synthetic, from Fischer Tropsch process, hydrogen from coal gasification, at fuelling station                                                             all IAM regions    van der Giesen et al. 2014
  Diesel production, synthetic, from Fischer Tropsch process, hydrogen from coal gasification, with CCS, at fuelling station                                                   all IAM regions    van der Giesen et al. 2014
  Diesel production, synthetic, from Fischer Tropsch process, hydrogen from electrolysis, at fuelling station                                                                  all IAM regions    van der Giesen et al. 2014
  Diesel production, synthetic, from Fischer Tropsch process, hydrogen from wood gasification, at fuelling station                                                             all IAM regions    van der Giesen et al. 2014
  Diesel production, synthetic, from Fischer Tropsch process, hydrogen from wood gasification, with CCS, at fuelling station                                                   all IAM regions    van der Giesen et al. 2014
  Diesel production, synthetic, from methanol, hydrogen from coal gasification, at fuelling station                                                                            all IAM regions    Hank et al, 2019
  Diesel production, synthetic, from methanol, hydrogen from coal gasification, with CCS, at fuelling station                                                                  all IAM regions    Hank et al, 2019
  Diesel production, synthetic, from methanol, hydrogen from electrolysis, CO2 from cement plant, at fuelling station                                                          all IAM regions    Hank et al, 2019
  Diesel production, synthetic, from methanol, hydrogen from electrolysis, CO2 from DAC, at fuelling station                                                                   all IAM regions    Hank et al, 2019
  Gasoline production, synthetic, from methanol, hydrogen from coal gasification, at fuelling station                                                                          all IAM regions    Hank et al, 2019
  Gasoline production, synthetic, from methanol, hydrogen from coal gasification, with CCS, at fuelling station                                                                all IAM regions    Hank et al, 2019
  Gasoline production, synthetic, from methanol, hydrogen from electrolysis, CO2 from cement plant, at fuelling station                                                        all IAM regions    Hank et al, 2019
  Gasoline production, synthetic, from methanol, hydrogen from electrolysis, CO2 from DAC, at fuelling station                                                                 all IAM regions    Hank et al, 2019
  Kerosene production, from methanol, hydrogen from coal gasification                                                                                                          all IAM regions    Hank et al, 2019
  Kerosene production, from methanol, hydrogen from electrolysis, CO2 from cement plant                                                                                        all IAM regions    Hank et al, 2019
  Kerosene production, from methanol, hydrogen from electrolysis, CO2 from DAC                                                                                                 all IAM regions    Hank et al, 2019
  Kerosene production, synthetic, Fischer Tropsch process, hydrogen from coal gasification                                                                                     all IAM regions    van der Giesen et al. 2014
  Kerosene production, synthetic, Fischer Tropsch process, hydrogen from coal gasification, with CCS                                                                           all IAM regions    van der Giesen et al. 2014
  Kerosene production, synthetic, Fischer Tropsch process, hydrogen from electrolysis                                                                                          all IAM regions    van der Giesen et al. 2014
  Kerosene production, synthetic, Fischer Tropsch process, hydrogen from wood gasification                                                                                     all IAM regions    van der Giesen et al. 2014
  Kerosene production, synthetic, Fischer Tropsch process, hydrogen from wood gasification, with CCS                                                                           all IAM regions    van der Giesen et al. 2014
  Lubricating oil production, synthetic, Fischer Tropsch process, hydrogen from coal gasification                                                                              all IAM regions    van der Giesen et al. 2014
  Lubricating oil production, synthetic, Fischer Tropsch process, hydrogen from electrolysis                                                                                   all IAM regions    van der Giesen et al. 2014
  Lubricating oil production, synthetic, Fischer Tropsch process, hydrogen from wood gasification                                                                              all IAM regions    van der Giesen et al. 2014
  Lubricating oil production, synthetic, Fischer Tropsch process, hydrogen from wood gasification, with CCS                                                                    all IAM regions    van der Giesen et al. 2014
  Methane, synthetic, gaseous, 5 bar, from coal-based hydrogen, at fuelling station                                                                                            all IAM regions    Zhang et al, 2019
  Methane, synthetic, gaseous, 5 bar, from electrochemical methanation (H2 from electrolysis, CO2 from DAC using heat pump heat), at fuelling station, using heat pump heat    all IAM regions    Zhang et al, 2019
  Methane, synthetic, gaseous, 5 bar, from electrochemical methanation (H2 from electrolysis, CO2 from DAC using waste heat), at fuelling station, using waste heat            all IAM regions    Zhang et al, 2019
  Methane, synthetic, gaseous, 5 bar, from electrochemical methanation, at fuelling station                                                                                    all IAM regions    Zhang et al, 2019
  Naphtha production, synthetic, Fischer Tropsch process, hydrogen from coal gasification                                                                                      all IAM regions    van der Giesen et al. 2014
  Naphtha production, synthetic, Fischer Tropsch process, hydrogen from electrolysis                                                                                           all IAM regions    van der Giesen et al. 2014
  Naphtha production, synthetic, Fischer Tropsch process, hydrogen from wood gasification                                                                                      all IAM regions    van der Giesen et al. 2014
  Naphtha production, synthetic, Fischer Tropsch process, hydrogen from wood gasification, with CCS                                                                            all IAM regions    van der Giesen et al. 2014
  Liquefied petroleum gas production, synthetic, from methanol, hydrogen from electrolysis, CO2 from DAC, at fuelling station                                                  all IAM regions    Hank et al, 2019
 ============================================================================================================================================================================ ================== =============================

In the case of wood and coal gasification-based fuels, the CO2 needed to produce methanol
or syncrude originates from the gasification process itself. This also implies
that in the methanol and/or RWGS process, a carbon balance correction is applied to reflect the
fact that a part of the CO2 from the gasification process is redirected into
the fuel production process.

If the CO2 originates from:

* a gasification process without CCS, a negative carbon correction is added to
reflect the fact that part of the CO2 has not been emitted but has ended in the fuel instead.
* the gasification process with CCS, no carbon correction is necessary, because the CO2 is stored
in the fuel instead of being stored underground, which from a carbon accounting standpoint is
similar.

Carbon Capture
--------------

Two sets of inventories for Direct Air Capture (DAC) are available in *premise*.
One for a solvent-based system, and one for a sorbent-based system. The inventories
were developed by Qiu_ and are available in the LCI_DAC_ spreadsheet. For each,
a variant including the subsequent compression, transport and storage of the
captured CO2 is also available.

They can be consulted here: LCI_DAC_.

.. _Qiu: https://doi.org/10.1038/s41467-022-31146-1
.. _LCI_DAC: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-direct-air-capture.xlsx

Additional, two datasets for carbon capture at point sources are available:
one at cement plant from Meunier_ et al, 2020, and another one at municipal solid waste incineration plant (MSWI)
from Bisinella_ et al, 2021.

.. _Meunier: https://doi.org/10.1016/j.renene.2019.07.010
.. _Bisinella: https://doi.org/10.1016/j.wasman.2021.04.046

They introduce the following datasets:

 =============================================================================================================== ===========
  Activity                                                                                                         Location
 =============================================================================================================== ===========
  carbon dioxide, captured from atmosphere, with a sorbent-based direct air capture system, 100ktCO2               RER
  carbon dioxide, captured from atmosphere and stored, with a sorbent-based direct air capture system, 100ktCO2    RER
  carbon dioxide, captured from atmosphere, with a solvent-based direct air capture system, 1MtCO2                 RER
  carbon dioxide, captured from atmosphere and stored, with a solvent-based direct air capture system, 1MtCO2      RER
  carbon dioxide, captured at municipal solid waste incineration plant, for subsequent reuse                       RER
  carbon dioxide, captured at cement production plant, for subsequent reuse                                        RER
 =============================================================================================================== ===========

Using the transformation function `update_dac()`, *premise* creates various configurations of these processes,
using different sources for heat (industrial steam heat, high-temp heat
pump heat and excess heat), which are found under the following names, for each IAM region:

 ======================================================================================================================================================= ==================
  name                                                                                                                                                      location
 ======================================================================================================================================================= ==================
  carbon dioxide, captured from atmosphere, with a solvent-based direct air capture system, 1MtCO2, with industrial steam heat, and grid electricity       all IAM regions
  carbon dioxide, captured from atmosphere, with a solvent-based direct air capture system, 1MtCO2, with heat pump heat, and grid electricity              all IAM regions
  carbon dioxide, captured from atmosphere, with a sorbent-based direct air capture system, 100ktCO2, with waste heat, and grid electricity                all IAM regions
  carbon dioxide, captured from atmosphere, with a sorbent-based direct air capture system, 100ktCO2, with industrial steam heat, and grid electricity     all IAM regions
  carbon dioxide, captured from atmosphere, with a sorbent-based direct air capture system, 100ktCO2, with heat pump heat, and grid electricity            all IAM regions
 ======================================================================================================================================================= ==================

Note that only solid sorbent DAC can use waste heat, as teh heat requirement for liquid solvent DAC
is too high (~900 C)

Li-ion batteries
----------------

NMC-111, NMC-6222 NMC-811 and NCA Lithium-ion battery inventories are originally
from Dai_ et al. 2019. They have been adapted to ecoinvent by Crenna_ et al, 2021.
LFP and LTO Lithium-ion battery inventories are from  Schmidt_ et al. 2019.

They introduce the following datasets:

 ============================== =========== ======================================
  Battery components             location    source
 ============================== =========== ======================================
  Battery BoP                    GLO         Schmidt et al. 2019
  Battery cell, NMC-111          GLO         Dai et al. 2019, Crenna et al. 2021
  Battery cell, NMC-622          GLO         Dai et al. 2019, Crenna et al. 2021
  Battery cell, NMC-811          GLO         Dai et al. 2019, Crenna et al. 2021
  Battery cell, NCA              GLO         Dai et al. 2019, Crenna et al. 2021
  Battery cell, LFP              GLO         Schmidt et al. 2019
  Battery cell, LTO              GLO         Schmidt et al. 2019
 ============================== =========== ======================================

These battery inventories are mostly used by battery electric vehicles
(also imported by *premise*), and are to be preferred to battery
inventories coming with ecoinvent (battery inventories since ecoinvent 3.8
are also from Crenna_ et al, 2021, but have been implemented with
some errors, which may be corrected in the future in ecoinvent 3.9).

These inventories can be found here: LCI_batteries_.

Graphite
--------

*premise* includes new inventories for:

* natural graphite, from Engels_ et al. 2022,
* synthetic graphite, from Surovtseva_ et al. 2022,

forming a new market for graphite, with the following datasets:

 ===================================== =========== ===========
  Activity                              Location
 ===================================== =========== ===========
  market for graphite, battery grade                1.0
  graphite, natural                     CN          0.8
  graphite, synthetic                   CN          0.2
 ===================================== =========== ===========

to represent a 80:20 split between natural and synthetic graphite,
according to Surovtseva_ et al, 2022.

These inventories can be found here: LCI_graphite_.

Cobalt
------

New inventories of cobalt are added, from the work of Dai, Kelly and Elgowainy_, 2018.
They are available under the following datasets:

=================================================================================== ===========
Activity                                                                             Location
=================================================================================== ===========
cobalt sulfate production, from copper mining, economic allocation                   CN
cobalt sulfate production, from copper mining, energy allocation                     CN
cobalt metal production, from copper mining, via electrolysis, economic allocation   CN
cobalt metal production, from copper mining, via electrolysis, energy allocation     CN
=================================================================================== ===========

We recommend using those rather than the original ecoinvent inventories for cobalt, provided
by the Cobalt Development Institute (CDI) since ecoinvent 3.7, which seem to lack transparency.

These inventories can be found here: LCI_cobalt_.

Lithium
-------

New inventories for lithium extraction are also added,
from the work of Schenker_ et al., 2022.
They cover lithium extraction from five different locations in Chile, Argentina and China.
They are available under the following datasets for battery production:

=================================================================================== ===========
Activity                                                                             Location
=================================================================================== ===========
market for lithium carbonate, battery grade                                          GLO
market for lithium hydroxide, battery grade                                          GLO
=================================================================================== ===========

These inventories can be found here: LCI_lithium_.

.. _Dai: https://www.mdpi.com/2313-0105/5/2/48
.. _Crenna: https://doi.org/10.1016/j.resconrec.2021.105619
.. _Schmidt: https://doi.org/10.1021/acs.est.8b05313
.. _Engels: https://doi.org/10.1016/j.jclepro.2022.130474
.. _Surovtseva: https://doi.org/10.1111/jiec.13234
.. _Elgowainy: https://greet.es.anl.gov/publication-update_cobalt
.. _Schenker: https://doi.org/10.1016/j.resconrec.2022.106611
.. _LCI_batteries: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-batteries.xlsx
.. _LCI_graphite: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-graphite.xlsx
.. _LCI_cobalt: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-cobalt.xlsx
.. _LCI_lithium: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-lithium.xlsx


Road vehicles
-------------

*premise* imports inventories for different types of on-road vehicles.

Two-wheelers
************

The following datasets for two-wheelers are imported.
Inventories are from Sacchi_ et al. 2022. The vehicles are available
for different years and emission standards. *premise* will only
import vehicles which production year is equal or inferior to
the scenario year considered. The inventories can be consulted
here: LCItwowheelers_.


.. _Sacchi: https://zenodo.org/deposit/5720779
.. _LCItwowheelers: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-two_wheelers.xlsx

 ================================================= ==================
  Two-wheeler datasets                              location
 ================================================= ==================
  transport, Kick Scooter, electric, <1kW           all IAM regions
  transport, Bicycle, conventional, urban           all IAM regions
  transport, Bicycle, electric (<25 km/h)           all IAM regions
  transport, Bicycle, electric (<45 km/h)           all IAM regions
  transport, Bicycle, electric, cargo bike          all IAM regions
  transport, Moped, gasoline, <4kW, EURO-3          all IAM regions
  transport, Moped, gasoline, <4kW, EURO-4          all IAM regions
  transport, Moped, gasoline, <4kW, EURO-5          all IAM regions
  transport, Scooter, gasoline, <4kW, EURO-3        all IAM regions
  transport, Scooter, gasoline, <4kW, EURO-4        all IAM regions
  transport, Scooter, gasoline, <4kW, EURO-5        all IAM regions
  transport, Scooter, gasoline, 4-11kW, EURO-3      all IAM regions
  transport, Scooter, gasoline, 4-11kW, EURO-4      all IAM regions
  transport, Scooter, gasoline, 4-11kW, EURO-5      all IAM regions
  transport, Scooter, electric, <4kW                all IAM regions
  transport, Scooter, electric, 4-11kW              all IAM regions
  transport, Motorbike, gasoline, 4-11kW, EURO-3    all IAM regions
  transport, Motorbike, gasoline, 4-11kW, EURO-4    all IAM regions
  transport, Motorbike, gasoline, 4-11kW, EURO-5    all IAM regions
  transport, Motorbike, gasoline, 11-35kW, EURO-3   all IAM regions
  transport, Motorbike, gasoline, 11-35kW, EURO-4   all IAM regions
  transport, Motorbike, gasoline, 11-35kW, EURO-5   all IAM regions
  transport, Motorbike, gasoline, >35kW, EURO-3     all IAM regions
  transport, Motorbike, gasoline, >35kW, EURO-4     all IAM regions
  transport, Motorbike, gasoline, >35kW, EURO-5     all IAM regions
  transport, Motorbike, electric, <4kW              all IAM regions
  transport, Motorbike, electric, 4-11kW            all IAM regions
  transport, Motorbike, electric, 11-35kW           all IAM regions
  transport, Motorbike, electric, >35kW             all IAM regions
 ================================================= ==================

These inventories do not supply inputs to other activities in the LCI database.
As such, they are optional.


Passenger cars
**************

The following datasets for passenger cars are imported.

 =============================================================================== ==================
  Passenger car datasets                                                          location
 =============================================================================== ==================
  transport, passenger car, gasoline, Large, EURO-2                               all IAM regions
  transport, passenger car, gasoline, Large, EURO-3                               all IAM regions
  transport, passenger car, gasoline, Large, EURO-4                               all IAM regions
  transport, passenger car, gasoline, Large, EURO-6ab                             all IAM regions
  transport, passenger car, gasoline, Large, EURO-6d-TEMP                         all IAM regions
  transport, passenger car, gasoline, Large, EURO-6d                              all IAM regions
  transport, passenger car, diesel, Large, EURO-2                                 all IAM regions
  transport, passenger car, diesel, Large, EURO-3                                 all IAM regions
  transport, passenger car, diesel, Large, EURO-4                                 all IAM regions
  transport, passenger car, diesel, Large, EURO-6ab                               all IAM regions
  transport, passenger car, diesel, Large, EURO-6d-TEMP                           all IAM regions
  transport, passenger car, diesel, Large, EURO-6d                                all IAM regions
  transport, passenger car, compressed gas, Large, EURO-2                         all IAM regions
  transport, passenger car, compressed gas, Large, EURO-3                         all IAM regions
  transport, passenger car, compressed gas, Large, EURO-4                         all IAM regions
  transport, passenger car, compressed gas, Large, EURO-6ab                       all IAM regions
  transport, passenger car, compressed gas, Large, EURO-6d-TEMP                   all IAM regions
  transport, passenger car, compressed gas, Large, EURO-6d                        all IAM regions
  transport, passenger car, plugin gasoline hybrid, Large, EURO-6ab               all IAM regions
  transport, passenger car, plugin gasoline hybrid, Large, EURO-6d-TEMP           all IAM regions
  transport, passenger car, plugin gasoline hybrid, Large, EURO-6d                all IAM regions
  transport, passenger car, plugin diesel hybrid, Large, EURO-6ab                 all IAM regions
  transport, passenger car, plugin diesel hybrid, Large, EURO-6d-TEMP             all IAM regions
  transport, passenger car, plugin diesel hybrid, Large, EURO-6d                  all IAM regions
  transport, passenger car, fuel cell electric, Large                             all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Large              all IAM regions
  transport, passenger car, gasoline hybrid, Large, EURO-6ab                      all IAM regions
  transport, passenger car, gasoline hybrid, Large, EURO-6d-TEMP                  all IAM regions
  transport, passenger car, gasoline hybrid, Large, EURO-6d                       all IAM regions
  transport, passenger car, diesel hybrid, Large, EURO-6ab                        all IAM regions
  transport, passenger car, diesel hybrid, Large, EURO-6d-TEMP                    all IAM regions
  transport, passenger car, diesel hybrid, Large, EURO-6d                         all IAM regions
  transport, passenger car, gasoline, Large SUV, EURO-2                           all IAM regions
  transport, passenger car, gasoline, Large SUV, EURO-3                           all IAM regions
  transport, passenger car, gasoline, Large SUV, EURO-4                           all IAM regions
  transport, passenger car, gasoline, Large SUV, EURO-6ab                         all IAM regions
  transport, passenger car, gasoline, Large SUV, EURO-6d-TEMP                     all IAM regions
  transport, passenger car, gasoline, Large SUV, EURO-6d                          all IAM regions
  transport, passenger car, diesel, Large SUV, EURO-2                             all IAM regions
  transport, passenger car, diesel, Large SUV, EURO-3                             all IAM regions
  transport, passenger car, diesel, Large SUV, EURO-4                             all IAM regions
  transport, passenger car, diesel, Large SUV, EURO-6ab                           all IAM regions
  transport, passenger car, diesel, Large SUV, EURO-6d-TEMP                       all IAM regions
  transport, passenger car, diesel, Large SUV, EURO-6d                            all IAM regions
  transport, passenger car, compressed gas, Large SUV, EURO-2                     all IAM regions
  transport, passenger car, compressed gas, Large SUV, EURO-3                     all IAM regions
  transport, passenger car, compressed gas, Large SUV, EURO-4                     all IAM regions
  transport, passenger car, compressed gas, Large SUV, EURO-6ab                   all IAM regions
  transport, passenger car, compressed gas, Large SUV, EURO-6d-TEMP               all IAM regions
  transport, passenger car, compressed gas, Large SUV, EURO-6d                    all IAM regions
  transport, passenger car, plugin gasoline hybrid, Large SUV, EURO-6ab           all IAM regions
  transport, passenger car, plugin gasoline hybrid, Large SUV, EURO-6d-TEMP       all IAM regions
  transport, passenger car, plugin gasoline hybrid, Large SUV, EURO-6d            all IAM regions
  transport, passenger car, plugin diesel hybrid, Large SUV, EURO-6ab             all IAM regions
  transport, passenger car, plugin diesel hybrid, Large SUV, EURO-6d-TEMP         all IAM regions
  transport, passenger car, plugin diesel hybrid, Large SUV, EURO-6d              all IAM regions
  transport, passenger car, fuel cell electric, Large SUV                         all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Large SUV          all IAM regions
  transport, passenger car, gasoline hybrid, Large SUV, EURO-6ab                  all IAM regions
  transport, passenger car, gasoline hybrid, Large SUV, EURO-6d-TEMP              all IAM regions
  transport, passenger car, gasoline hybrid, Large SUV, EURO-6d                   all IAM regions
  transport, passenger car, diesel hybrid, Large SUV, EURO-6ab                    all IAM regions
  transport, passenger car, diesel hybrid, Large SUV, EURO-6d-TEMP                all IAM regions
  transport, passenger car, diesel hybrid, Large SUV, EURO-6d                     all IAM regions
  transport, passenger car, gasoline, Lower medium, EURO-2                        all IAM regions
  transport, passenger car, gasoline, Lower medium, EURO-3                        all IAM regions
  transport, passenger car, gasoline, Lower medium, EURO-4                        all IAM regions
  transport, passenger car, gasoline, Lower medium, EURO-6ab                      all IAM regions
  transport, passenger car, gasoline, Lower medium, EURO-6d-TEMP                  all IAM regions
  transport, passenger car, gasoline, Lower medium, EURO-6d                       all IAM regions
  transport, passenger car, diesel, Lower medium, EURO-2                          all IAM regions
  transport, passenger car, diesel, Lower medium, EURO-3                          all IAM regions
  transport, passenger car, diesel, Lower medium, EURO-4                          all IAM regions
  transport, passenger car, diesel, Lower medium, EURO-6ab                        all IAM regions
  transport, passenger car, diesel, Lower medium, EURO-6d-TEMP                    all IAM regions
  transport, passenger car, diesel, Lower medium, EURO-6d                         all IAM regions
  transport, passenger car, compressed gas, Lower medium, EURO-2                  all IAM regions
  transport, passenger car, compressed gas, Lower medium, EURO-3                  all IAM regions
  transport, passenger car, compressed gas, Lower medium, EURO-4                  all IAM regions
  transport, passenger car, compressed gas, Lower medium, EURO-6ab                all IAM regions
  transport, passenger car, compressed gas, Lower medium, EURO-6d-TEMP            all IAM regions
  transport, passenger car, compressed gas, Lower medium, EURO-6d                 all IAM regions
  transport, passenger car, plugin gasoline hybrid, Lower medium, EURO-6ab        all IAM regions
  transport, passenger car, plugin gasoline hybrid, Lower medium, EURO-6d-TEMP    all IAM regions
  transport, passenger car, plugin gasoline hybrid, Lower medium, EURO-6d         all IAM regions
  transport, passenger car, plugin diesel hybrid, Lower medium, EURO-6ab          all IAM regions
  transport, passenger car, plugin diesel hybrid, Lower medium, EURO-6d-TEMP      all IAM regions
  transport, passenger car, plugin diesel hybrid, Lower medium, EURO-6d           all IAM regions
  transport, passenger car, fuel cell electric, Lower medium                      all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Lower medium       all IAM regions
  transport, passenger car, gasoline hybrid, Lower medium, EURO-6ab               all IAM regions
  transport, passenger car, gasoline hybrid, Lower medium, EURO-6d-TEMP           all IAM regions
  transport, passenger car, gasoline hybrid, Lower medium, EURO-6d                all IAM regions
  transport, passenger car, diesel hybrid, Lower medium, EURO-6ab                 all IAM regions
  transport, passenger car, diesel hybrid, Lower medium, EURO-6d-TEMP             all IAM regions
  transport, passenger car, diesel hybrid, Lower medium, EURO-6d                  all IAM regions
  transport, passenger car, gasoline, Medium, EURO-2                              all IAM regions
  transport, passenger car, gasoline, Medium, EURO-3                              all IAM regions
  transport, passenger car, gasoline, Medium, EURO-4                              all IAM regions
  transport, passenger car, gasoline, Medium, EURO-6ab                            all IAM regions
  transport, passenger car, gasoline, Medium, EURO-6d-TEMP                        all IAM regions
  transport, passenger car, gasoline, Medium, EURO-6d                             all IAM regions
  transport, passenger car, diesel, Medium, EURO-2                                all IAM regions
  transport, passenger car, diesel, Medium, EURO-3                                all IAM regions
  transport, passenger car, diesel, Medium, EURO-4                                all IAM regions
  transport, passenger car, diesel, Medium, EURO-6ab                              all IAM regions
  transport, passenger car, diesel, Medium, EURO-6d-TEMP                          all IAM regions
  transport, passenger car, diesel, Medium, EURO-6d                               all IAM regions
  transport, passenger car, compressed gas, Medium, EURO-2                        all IAM regions
  transport, passenger car, compressed gas, Medium, EURO-3                        all IAM regions
  transport, passenger car, compressed gas, Medium, EURO-4                        all IAM regions
  transport, passenger car, compressed gas, Medium, EURO-6ab                      all IAM regions
  transport, passenger car, compressed gas, Medium, EURO-6d-TEMP                  all IAM regions
  transport, passenger car, compressed gas, Medium, EURO-6d                       all IAM regions
  transport, passenger car, plugin gasoline hybrid, Medium, EURO-6ab              all IAM regions
  transport, passenger car, plugin gasoline hybrid, Medium, EURO-6d-TEMP          all IAM regions
  transport, passenger car, plugin gasoline hybrid, Medium, EURO-6d               all IAM regions
  transport, passenger car, plugin diesel hybrid, Medium, EURO-6ab                all IAM regions
  transport, passenger car, plugin diesel hybrid, Medium, EURO-6d-TEMP            all IAM regions
  transport, passenger car, plugin diesel hybrid, Medium, EURO-6d                 all IAM regions
  transport, passenger car, fuel cell electric, Medium                            all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Medium             all IAM regions
  transport, passenger car, gasoline hybrid, Medium, EURO-6ab                     all IAM regions
  transport, passenger car, gasoline hybrid, Medium, EURO-6d-TEMP                 all IAM regions
  transport, passenger car, gasoline hybrid, Medium, EURO-6d                      all IAM regions
  transport, passenger car, diesel hybrid, Medium, EURO-6ab                       all IAM regions
  transport, passenger car, diesel hybrid, Medium, EURO-6d-TEMP                   all IAM regions
  transport, passenger car, diesel hybrid, Medium, EURO-6d                        all IAM regions
  transport, passenger car, gasoline, Medium SUV, EURO-2                          all IAM regions
  transport, passenger car, gasoline, Medium SUV, EURO-3                          all IAM regions
  transport, passenger car, gasoline, Medium SUV, EURO-4                          all IAM regions
  transport, passenger car, gasoline, Medium SUV, EURO-6ab                        all IAM regions
  transport, passenger car, gasoline, Medium SUV, EURO-6d-TEMP                    all IAM regions
  transport, passenger car, gasoline, Medium SUV, EURO-6d                         all IAM regions
  transport, passenger car, diesel, Medium SUV, EURO-2                            all IAM regions
  transport, passenger car, diesel, Medium SUV, EURO-3                            all IAM regions
  transport, passenger car, diesel, Medium SUV, EURO-4                            all IAM regions
  transport, passenger car, diesel, Medium SUV, EURO-6ab                          all IAM regions
  transport, passenger car, diesel, Medium SUV, EURO-6d-TEMP                      all IAM regions
  transport, passenger car, diesel, Medium SUV, EURO-6d                           all IAM regions
  transport, passenger car, compressed gas, Medium SUV, EURO-2                    all IAM regions
  transport, passenger car, compressed gas, Medium SUV, EURO-3                    all IAM regions
  transport, passenger car, compressed gas, Medium SUV, EURO-4                    all IAM regions
  transport, passenger car, compressed gas, Medium SUV, EURO-6ab                  all IAM regions
  transport, passenger car, compressed gas, Medium SUV, EURO-6d-TEMP              all IAM regions
  transport, passenger car, compressed gas, Medium SUV, EURO-6d                   all IAM regions
  transport, passenger car, plugin gasoline hybrid, Medium SUV, EURO-6ab          all IAM regions
  transport, passenger car, plugin gasoline hybrid, Medium SUV, EURO-6d-TEMP      all IAM regions
  transport, passenger car, plugin gasoline hybrid, Medium SUV, EURO-6d           all IAM regions
  transport, passenger car, plugin diesel hybrid, Medium SUV, EURO-6ab            all IAM regions
  transport, passenger car, plugin diesel hybrid, Medium SUV, EURO-6d-TEMP        all IAM regions
  transport, passenger car, plugin diesel hybrid, Medium SUV, EURO-6d             all IAM regions
  transport, passenger car, fuel cell electric, Medium SUV                        all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Medium SUV         all IAM regions
  transport, passenger car, gasoline hybrid, Medium SUV, EURO-6ab                 all IAM regions
  transport, passenger car, gasoline hybrid, Medium SUV, EURO-6d-TEMP             all IAM regions
  transport, passenger car, gasoline hybrid, Medium SUV, EURO-6d                  all IAM regions
  transport, passenger car, diesel hybrid, Medium SUV, EURO-6ab                   all IAM regions
  transport, passenger car, diesel hybrid, Medium SUV, EURO-6d-TEMP               all IAM regions
  transport, passenger car, diesel hybrid, Medium SUV, EURO-6d                    all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Micro              all IAM regions
  transport, passenger car, gasoline, Mini, EURO-2                                all IAM regions
  transport, passenger car, gasoline, Mini, EURO-3                                all IAM regions
  transport, passenger car, gasoline, Mini, EURO-4                                all IAM regions
  transport, passenger car, gasoline, Mini, EURO-6ab                              all IAM regions
  transport, passenger car, gasoline, Mini, EURO-6d-TEMP                          all IAM regions
  transport, passenger car, gasoline, Mini, EURO-6d                               all IAM regions
  transport, passenger car, diesel, Mini, EURO-2                                  all IAM regions
  transport, passenger car, diesel, Mini, EURO-3                                  all IAM regions
  transport, passenger car, diesel, Mini, EURO-4                                  all IAM regions
  transport, passenger car, diesel, Mini, EURO-6ab                                all IAM regions
  transport, passenger car, diesel, Mini, EURO-6d-TEMP                            all IAM regions
  transport, passenger car, diesel, Mini, EURO-6d                                 all IAM regions
  transport, passenger car, compressed gas, Mini, EURO-2                          all IAM regions
  transport, passenger car, compressed gas, Mini, EURO-3                          all IAM regions
  transport, passenger car, compressed gas, Mini, EURO-4                          all IAM regions
  transport, passenger car, compressed gas, Mini, EURO-6ab                        all IAM regions
  transport, passenger car, compressed gas, Mini, EURO-6d-TEMP                    all IAM regions
  transport, passenger car, compressed gas, Mini, EURO-6d                         all IAM regions
  transport, passenger car, plugin gasoline hybrid, Mini, EURO-6ab                all IAM regions
  transport, passenger car, plugin gasoline hybrid, Mini, EURO-6d-TEMP            all IAM regions
  transport, passenger car, plugin gasoline hybrid, Mini, EURO-6d                 all IAM regions
  transport, passenger car, plugin diesel hybrid, Mini, EURO-6ab                  all IAM regions
  transport, passenger car, plugin diesel hybrid, Mini, EURO-6d-TEMP              all IAM regions
  transport, passenger car, plugin diesel hybrid, Mini, EURO-6d                   all IAM regions
  transport, passenger car, fuel cell electric, Mini                              all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Mini               all IAM regions
  transport, passenger car, gasoline hybrid, Mini, EURO-6ab                       all IAM regions
  transport, passenger car, gasoline hybrid, Mini, EURO-6d-TEMP                   all IAM regions
  transport, passenger car, gasoline hybrid, Mini, EURO-6d                        all IAM regions
  transport, passenger car, diesel hybrid, Mini, EURO-6ab                         all IAM regions
  transport, passenger car, diesel hybrid, Mini, EURO-6d-TEMP                     all IAM regions
  transport, passenger car, diesel hybrid, Mini, EURO-6d                          all IAM regions
  transport, passenger car, gasoline, Small, EURO-2                               all IAM regions
  transport, passenger car, gasoline, Small, EURO-3                               all IAM regions
  transport, passenger car, gasoline, Small, EURO-4                               all IAM regions
  transport, passenger car, gasoline, Small, EURO-6ab                             all IAM regions
  transport, passenger car, gasoline, Small, EURO-6d-TEMP                         all IAM regions
  transport, passenger car, gasoline, Small, EURO-6d                              all IAM regions
  transport, passenger car, diesel, Small, EURO-2                                 all IAM regions
  transport, passenger car, diesel, Small, EURO-3                                 all IAM regions
  transport, passenger car, diesel, Small, EURO-4                                 all IAM regions
  transport, passenger car, diesel, Small, EURO-6ab                               all IAM regions
  transport, passenger car, diesel, Small, EURO-6d-TEMP                           all IAM regions
  transport, passenger car, diesel, Small, EURO-6d                                all IAM regions
  transport, passenger car, compressed gas, Small, EURO-2                         all IAM regions
  transport, passenger car, compressed gas, Small, EURO-3                         all IAM regions
  transport, passenger car, compressed gas, Small, EURO-4                         all IAM regions
  transport, passenger car, compressed gas, Small, EURO-6ab                       all IAM regions
  transport, passenger car, compressed gas, Small, EURO-6d-TEMP                   all IAM regions
  transport, passenger car, compressed gas, Small, EURO-6d                        all IAM regions
  transport, passenger car, plugin gasoline hybrid, Small, EURO-6ab               all IAM regions
  transport, passenger car, plugin gasoline hybrid, Small, EURO-6d-TEMP           all IAM regions
  transport, passenger car, plugin gasoline hybrid, Small, EURO-6d                all IAM regions
  transport, passenger car, plugin diesel hybrid, Small, EURO-6ab                 all IAM regions
  transport, passenger car, plugin diesel hybrid, Small, EURO-6d-TEMP             all IAM regions
  transport, passenger car, plugin diesel hybrid, Small, EURO-6d                  all IAM regions
  transport, passenger car, fuel cell electric, Small                             all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Small              all IAM regions
  transport, passenger car, gasoline hybrid, Small, EURO-6ab                      all IAM regions
  transport, passenger car, gasoline hybrid, Small, EURO-6d-TEMP                  all IAM regions
  transport, passenger car, gasoline hybrid, Small, EURO-6d                       all IAM regions
  transport, passenger car, diesel hybrid, Small, EURO-6ab                        all IAM regions
  transport, passenger car, diesel hybrid, Small, EURO-6d-TEMP                    all IAM regions
  transport, passenger car, diesel hybrid, Small, EURO-6d                         all IAM regions
  transport, passenger car, gasoline, Van, EURO-2                                 all IAM regions
  transport, passenger car, gasoline, Van, EURO-3                                 all IAM regions
  transport, passenger car, gasoline, Van, EURO-4                                 all IAM regions
  transport, passenger car, gasoline, Van, EURO-6ab                               all IAM regions
  transport, passenger car, gasoline, Van, EURO-6d-TEMP                           all IAM regions
  transport, passenger car, gasoline, Van, EURO-6d                                all IAM regions
  transport, passenger car, diesel, Van, EURO-2                                   all IAM regions
  transport, passenger car, diesel, Van, EURO-3                                   all IAM regions
  transport, passenger car, diesel, Van, EURO-4                                   all IAM regions
  transport, passenger car, diesel, Van, EURO-6ab                                 all IAM regions
  transport, passenger car, diesel, Van, EURO-6d-TEMP                             all IAM regions
  transport, passenger car, diesel, Van, EURO-6d                                  all IAM regions
  transport, passenger car, compressed gas, Van, EURO-2                           all IAM regions
  transport, passenger car, compressed gas, Van, EURO-3                           all IAM regions
  transport, passenger car, compressed gas, Van, EURO-4                           all IAM regions
  transport, passenger car, compressed gas, Van, EURO-6ab                         all IAM regions
  transport, passenger car, compressed gas, Van, EURO-6d-TEMP                     all IAM regions
  transport, passenger car, compressed gas, Van, EURO-6d                          all IAM regions
  transport, passenger car, plugin gasoline hybrid, Van, EURO-6ab                 all IAM regions
  transport, passenger car, plugin gasoline hybrid, Van, EURO-6d-TEMP             all IAM regions
  transport, passenger car, plugin gasoline hybrid, Van, EURO-6d                  all IAM regions
  transport, passenger car, plugin diesel hybrid, Van, EURO-6ab                   all IAM regions
  transport, passenger car, plugin diesel hybrid, Van, EURO-6d-TEMP               all IAM regions
  transport, passenger car, plugin diesel hybrid, Van, EURO-6d                    all IAM regions
  transport, passenger car, fuel cell electric, Van                               all IAM regions
  transport, passenger car, battery electric, NMC-622 battery, Van                all IAM regions
  transport, passenger car, gasoline hybrid, Van, EURO-6ab                        all IAM regions
  transport, passenger car, gasoline hybrid, Van, EURO-6d-TEMP                    all IAM regions
  transport, passenger car, gasoline hybrid, Van, EURO-6d                         all IAM regions
  transport, passenger car, diesel hybrid, Van, EURO-6ab                          all IAM regions
  transport, passenger car, diesel hybrid, Van, EURO-6d-TEMP                      all IAM regions
  transport, passenger car, diesel hybrid, Van, EURO-6d                           all IAM regions
 =============================================================================== ==================

Inventories are from Sacchi2_ et al. 2022. The vehicles are available
for different years and emission standards and for each IAM region. *premise* will only
import vehicles which production year is equal or inferior to
the scenario year considered. *premise* will create fleet average vehicles
during the *Transport* transformation for each IAM region. The inventories can be consulted
here: LCIpasscars_.

.. _Sacchi2: https://www.psi.ch/en/media/72391/download
.. _LCIpasscars: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-pass_cars.xlsx

At the moment. these inventories do not supply inputs to other activities in the LCI database.
As such, they are optional.


Medium and heavy duty trucks
****************************

The following datasets for medium and heavy-duty trucks are imported.

 ================================================================================== ==================
  Truck datasets                                                                     location
 ================================================================================== ==================
  transport, freight, lorry, battery electric, NMC-622 battery, 3.5t gross weight    all IAM regions
  transport, freight, lorry, fuel cell electric, 3.5t gross weight                   all IAM regions
  transport, freight, lorry, diesel hybrid, 3.5t gross weight, EURO-VI               all IAM regions
  transport, freight, lorry, diesel, 3.5t gross weight, EURO-III                     all IAM regions
  transport, freight, lorry, diesel, 3.5t gross weight, EURO-IV                      all IAM regions
  transport, freight, lorry, diesel, 3.5t gross weight, EURO-V                       all IAM regions
  transport, freight, lorry, diesel, 3.5t gross weight, EURO-VI                      all IAM regions
  transport, freight, lorry, compressed gas, 3.5t gross weight, EURO-III             all IAM regions
  transport, freight, lorry, compressed gas, 3.5t gross weight, EURO-IV              all IAM regions
  transport, freight, lorry, compressed gas, 3.5t gross weight, EURO-V               all IAM regions
  transport, freight, lorry, compressed gas, 3.5t gross weight, EURO-VI              all IAM regions
  transport, freight, lorry, plugin diesel hybrid, 3.5t gross weight, EURO-VI        all IAM regions
  transport, freight, lorry, battery electric, NMC-622 battery, 7.5t gross weight    all IAM regions
  transport, freight, lorry, fuel cell electric, 7.5t gross weight                   all IAM regions
  transport, freight, lorry, diesel hybrid, 7.5t gross weight, EURO-VI               all IAM regions
  transport, freight, lorry, diesel, 7.5t gross weight, EURO-III                     all IAM regions
  transport, freight, lorry, diesel, 7.5t gross weight, EURO-IV                      all IAM regions
  transport, freight, lorry, diesel, 7.5t gross weight, EURO-V                       all IAM regions
  transport, freight, lorry, diesel, 7.5t gross weight, EURO-VI                      all IAM regions
  transport, freight, lorry, compressed gas, 7.5t gross weight, EURO-III             all IAM regions
  transport, freight, lorry, compressed gas, 7.5t gross weight, EURO-IV              all IAM regions
  transport, freight, lorry, compressed gas, 7.5t gross weight, EURO-V               all IAM regions
  transport, freight, lorry, compressed gas, 7.5t gross weight, EURO-VI              all IAM regions
  transport, freight, lorry, plugin diesel hybrid, 7.5t gross weight, EURO-VI        all IAM regions
  transport, freight, lorry, battery electric, NMC-622 battery, 18t gross weight     all IAM regions
  transport, freight, lorry, fuel cell electric, 18t gross weight                    all IAM regions
  transport, freight, lorry, diesel hybrid, 18t gross weight, EURO-VI                all IAM regions
  transport, freight, lorry, diesel, 18t gross weight, EURO-III                      all IAM regions
  transport, freight, lorry, diesel, 18t gross weight, EURO-IV                       all IAM regions
  transport, freight, lorry, diesel, 18t gross weight, EURO-V                        all IAM regions
  transport, freight, lorry, diesel, 18t gross weight, EURO-VI                       all IAM regions
  transport, freight, lorry, compressed gas, 18t gross weight, EURO-III              all IAM regions
  transport, freight, lorry, compressed gas, 18t gross weight, EURO-IV               all IAM regions
  transport, freight, lorry, compressed gas, 18t gross weight, EURO-V                all IAM regions
  transport, freight, lorry, compressed gas, 18t gross weight, EURO-VI               all IAM regions
  transport, freight, lorry, plugin diesel hybrid, 18t gross weight, EURO-VI         all IAM regions
  transport, freight, lorry, battery electric, NMC-622 battery, 26t gross weight     all IAM regions
  transport, freight, lorry, fuel cell electric, 26t gross weight                    all IAM regions
  transport, freight, lorry, diesel hybrid, 26t gross weight, EURO-VI                all IAM regions
  transport, freight, lorry, diesel, 26t gross weight, EURO-III                      all IAM regions
  transport, freight, lorry, diesel, 26t gross weight, EURO-IV                       all IAM regions
  transport, freight, lorry, diesel, 26t gross weight, EURO-V                        all IAM regions
  transport, freight, lorry, diesel, 26t gross weight, EURO-VI                       all IAM regions
  transport, freight, lorry, compressed gas, 26t gross weight, EURO-III              all IAM regions
  transport, freight, lorry, compressed gas, 26t gross weight, EURO-IV               all IAM regions
  transport, freight, lorry, compressed gas, 26t gross weight, EURO-V                all IAM regions
  transport, freight, lorry, compressed gas, 26t gross weight, EURO-VI               all IAM regions
  transport, freight, lorry, plugin diesel hybrid, 26t gross weight, EURO-VI         all IAM regions
  transport, freight, lorry, battery electric, NMC-622 battery, 32t gross weight     all IAM regions
  transport, freight, lorry, fuel cell electric, 32t gross weight                    all IAM regions
  transport, freight, lorry, diesel hybrid, 32t gross weight, EURO-VI                all IAM regions
  transport, freight, lorry, diesel, 32t gross weight, EURO-III                      all IAM regions
  transport, freight, lorry, diesel, 32t gross weight, EURO-IV                       all IAM regions
  transport, freight, lorry, diesel, 32t gross weight, EURO-V                        all IAM regions
  transport, freight, lorry, diesel, 32t gross weight, EURO-VI                       all IAM regions
  transport, freight, lorry, compressed gas, 32t gross weight, EURO-III              all IAM regions
  transport, freight, lorry, compressed gas, 32t gross weight, EURO-IV               all IAM regions
  transport, freight, lorry, compressed gas, 32t gross weight, EURO-V                all IAM regions
  transport, freight, lorry, compressed gas, 32t gross weight, EURO-VI               all IAM regions
  transport, freight, lorry, plugin diesel hybrid, 32t gross weight, EURO-VI         all IAM regions
  transport, freight, lorry, battery electric, NMC-622 battery, 40t gross weight     all IAM regions
  transport, freight, lorry, fuel cell electric, 40t gross weight                    all IAM regions
  transport, freight, lorry, diesel hybrid, 40t gross weight, EURO-VI                all IAM regions
  transport, freight, lorry, diesel, 40t gross weight, EURO-III                      all IAM regions
  transport, freight, lorry, diesel, 40t gross weight, EURO-IV                       all IAM regions
  transport, freight, lorry, diesel, 40t gross weight, EURO-V                        all IAM regions
  transport, freight, lorry, diesel, 40t gross weight, EURO-VI                       all IAM regions
  transport, freight, lorry, compressed gas, 40t gross weight, EURO-III              all IAM regions
  transport, freight, lorry, compressed gas, 40t gross weight, EURO-IV               all IAM regions
  transport, freight, lorry, compressed gas, 40t gross weight, EURO-V                all IAM regions
  transport, freight, lorry, compressed gas, 40t gross weight, EURO-VI               all IAM regions
  transport, freight, lorry, plugin diesel hybrid, 40t gross weight, EURO-VI         all IAM regions
 ================================================================================== ==================


Inventories are from Sacchi3_ et al. 2021. The vehicles are available
for different years and emission standards and for each IAM region. *premise* will only
import vehicles which production year is equal or inferior to
the scenario year considered. *premise* will create fleet average vehicles
during the *Transport* transformation for each IAM region. The inventories can be consulted
here: LCItrucks_.

.. _LCItrucks: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-trucks.xlsx
.. _Sacchi3: https://pubs.acs.org/doi/abs/10.1021/acs.est.0c07773



Buses
*****

The following datasets for city and coach buses are imported.

  =================================================================================================================== ==================
  Bus datasets                                                                                                        location
 =================================================================================================================== ==================
  transport, passenger bus, battery electric - overnight charging, NMC-622 battery, 9m midibus                        all IAM regions
  transport, passenger bus, battery electric - opportunity charging, LTO battery, 9m midibus                          all IAM regions
  transport, passenger bus, fuel cell electric, 9m midibus                                                            all IAM regions
  transport, passenger bus, diesel hybrid, 9m midibus, EURO-VI                                                        all IAM regions
  transport, passenger bus, diesel, 9m midibus, EURO-III                                                              all IAM regions
  transport, passenger bus, diesel, 9m midibus, EURO-IV                                                               all IAM regions
  transport, passenger bus, diesel, 9m midibus, EURO-V                                                                all IAM regions
  transport, passenger bus, diesel, 9m midibus, EURO-VI                                                               all IAM regions
  transport, passenger bus, compressed gas, 9m midibus, EURO-III                                                      all IAM regions
  transport, passenger bus, compressed gas, 9m midibus, EURO-IV                                                       all IAM regions
  transport, passenger bus, compressed gas, 9m midibus, EURO-V                                                        all IAM regions
  transport, passenger bus, compressed gas, 9m midibus, EURO-VI                                                       all IAM regions
  transport, passenger bus, battery electric - overnight charging, NMC-622 battery, 13m single deck urban bus         all IAM regions
  transport, passenger bus, battery electric - battery-equipped trolleybus, LTO battery, 13m single deck urban bus    all IAM regions
  transport, passenger bus, battery electric - opportunity charging, LTO battery, 13m single deck urban bus           all IAM regions
  transport, passenger bus, fuel cell electric, 13m single deck urban bus                                             all IAM regions
  transport, passenger bus, diesel hybrid, 13m single deck urban bus, EURO-VI                                         all IAM regions
  transport, passenger bus, diesel, 13m single deck urban bus, EURO-III                                               all IAM regions
  transport, passenger bus, diesel, 13m single deck urban bus, EURO-IV                                                all IAM regions
  transport, passenger bus, diesel, 13m single deck urban bus, EURO-V                                                 all IAM regions
  transport, passenger bus, diesel, 13m single deck urban bus, EURO-VI                                                all IAM regions
  transport, passenger bus, compressed gas, 13m single deck urban bus, EURO-III                                       all IAM regions
  transport, passenger bus, compressed gas, 13m single deck urban bus, EURO-IV                                        all IAM regions
  transport, passenger bus, compressed gas, 13m single deck urban bus, EURO-V                                         all IAM regions
  transport, passenger bus, compressed gas, 13m single deck urban bus, EURO-VI                                        all IAM regions
  transport, passenger bus, fuel cell electric, 13m single deck coach bus                                             all IAM regions
  transport, passenger bus, diesel hybrid, 13m single deck coach bus, EURO-VI                                         all IAM regions
  transport, passenger bus, diesel, 13m single deck coach bus, EURO-III                                               all IAM regions
  transport, passenger bus, diesel, 13m single deck coach bus, EURO-IV                                                all IAM regions
  transport, passenger bus, diesel, 13m single deck coach bus, EURO-V                                                 all IAM regions
  transport, passenger bus, diesel, 13m single deck coach bus, EURO-VI                                                all IAM regions
  transport, passenger bus, compressed gas, 13m single deck coach bus, EURO-III                                       all IAM regions
  transport, passenger bus, compressed gas, 13m single deck coach bus, EURO-IV                                        all IAM regions
  transport, passenger bus, compressed gas, 13m single deck coach bus, EURO-V                                         all IAM regions
  transport, passenger bus, compressed gas, 13m single deck coach bus, EURO-VI                                        all IAM regions
  transport, passenger bus, battery electric - overnight charging, NMC-622 battery, 13m double deck urban bus         all IAM regions
  transport, passenger bus, battery electric - opportunity charging, LTO battery, 13m double deck urban bus           all IAM regions
  transport, passenger bus, fuel cell electric, 13m double deck urban bus                                             all IAM regions
  transport, passenger bus, diesel hybrid, 13m double deck urban bus, EURO-VI                                         all IAM regions
  transport, passenger bus, diesel, 13m double deck urban bus, EURO-III                                               all IAM regions
  transport, passenger bus, diesel, 13m double deck urban bus, EURO-IV                                                all IAM regions
  transport, passenger bus, diesel, 13m double deck urban bus, EURO-V                                                 all IAM regions
  transport, passenger bus, diesel, 13m double deck urban bus, EURO-VI                                                all IAM regions
  transport, passenger bus, compressed gas, 13m double deck urban bus, EURO-III                                       all IAM regions
  transport, passenger bus, compressed gas, 13m double deck urban bus, EURO-IV                                        all IAM regions
  transport, passenger bus, compressed gas, 13m double deck urban bus, EURO-V                                         all IAM regions
  transport, passenger bus, compressed gas, 13m double deck urban bus, EURO-VI                                        all IAM regions
  transport, passenger bus, fuel cell electric, 13m double deck coach bus                                             all IAM regions
  transport, passenger bus, diesel hybrid, 13m double deck coach bus, EURO-VI                                         all IAM regions
  transport, passenger bus, diesel, 13m double deck coach bus, EURO-III                                               all IAM regions
  transport, passenger bus, diesel, 13m double deck coach bus, EURO-IV                                                all IAM regions
  transport, passenger bus, diesel, 13m double deck coach bus, EURO-V                                                 all IAM regions
  transport, passenger bus, diesel, 13m double deck coach bus, EURO-VI                                                all IAM regions
  transport, passenger bus, compressed gas, 13m double deck coach bus, EURO-III                                       all IAM regions
  transport, passenger bus, compressed gas, 13m double deck coach bus, EURO-IV                                        all IAM regions
  transport, passenger bus, compressed gas, 13m double deck coach bus, EURO-V                                         all IAM regions
  transport, passenger bus, compressed gas, 13m double deck coach bus, EURO-VI                                        all IAM regions
  transport, passenger bus, battery electric - overnight charging, NMC-622 battery, 18m articulated urban bus         all IAM regions
  transport, passenger bus, battery electric - battery-equipped trolleybus, LTO battery, 18m articulated urban bus    all IAM regions
  transport, passenger bus, battery electric - opportunity charging, LTO battery, 18m articulated urban bus           all IAM regions
  transport, passenger bus, fuel cell electric, 18m articulated urban bus                                             all IAM regions
  transport, passenger bus, diesel hybrid, 18m articulated urban bus, EURO-VI                                         all IAM regions
  transport, passenger bus, diesel, 18m articulated urban bus, EURO-III                                               all IAM regions
  transport, passenger bus, diesel, 18m articulated urban bus, EURO-IV                                                all IAM regions
  transport, passenger bus, diesel, 18m articulated urban bus, EURO-V                                                 all IAM regions
  transport, passenger bus, diesel, 18m articulated urban bus, EURO-VI                                                all IAM regions
  transport, passenger bus, compressed gas, 18m articulated urban bus, EURO-III                                       all IAM regions
  transport, passenger bus, compressed gas, 18m articulated urban bus, EURO-IV                                        all IAM regions
  transport, passenger bus, compressed gas, 18m articulated urban bus, EURO-V                                         all IAM regions
  transport, passenger bus, compressed gas, 18m articulated urban bus, EURO-VI                                        all IAM regions
 =================================================================================================================== ==================

Inventories are from Sacchi_ et al. 2021. The vehicles are available
for different years and emission standards and for each IAM region. *premise* will only
import vehicles which production year is equal or inferior to
the scenario year considered. *premise* will create fleet average vehicles
during the *Transport* transformation for each IAM region. The inventories can be consulted
here: LCIbuses_.

.. _LCIbuses: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-buses.xlsx

At the moment. these inventories do not supply inputs to other activities in the LCI database.
As such, they are optional.


Migration between ecoinvent versions
------------------------------------

Because the additional inventories that are imported may be composed
of exchanges meant to link with an ecoinvent version different
than what the user specifies to *premise* upon the database creation,
it is necessary to be able to "translate" the imported inventories
so that they correctly link to any ecoinvent version *premise* is
compatible with.

Therefore, *premise* has a migration map that is used to convert
certain exchanges to be compatible with a given ecoinvent version.

This migration map is provided here: migrationmap_.

.. _migrationmap: https://github.com/polca/premise/blob/master/premise/data/additional_inventories/migration_map.csv

IAM data collection
"""""""""""""""""""

After extracting the ecoinvent database and additional inventories,
*premise* instantiates the class *IAMDataCollection*, which collects
all sorts of data from the IAM output file and store it into
multi-dimensional arrays.


Production volumes
------------------

Production volumes for different commodities are collected, for the
year and scenario specified by the user. Production volumes are used to
build regional markets. For example, for the global market, the volume-based
shares of each region are used to reflect their respective supply importance.
Another example is for building electricity markets: the respective
production volumes of each electricity-producing technology is used to
determine the gross supply mix of the market.


The table below shows the correspondence between *premise*, REMIND, IMAGE
and LCI terminology, regarding electricity producing technologies. *premise*
production volumes given for secondary energy carriers for electricity.
The mapping file is available in the library root folder: mappingElec_.

.. _mappingElec: https://github.com/polca/premise/blob/master/premise/data/electricity/electricity_tech_vars.yml


 ========================== ===================================== ================================================= ===================================================================================================
  name in premise            name in REMIND                         name in IMAGE                                    name in LCI database (only first of several shown)
 ========================== ===================================== ================================================= ===================================================================================================
  Biomass CHP                SE|Electricity|Biomass|CHP|w/o CCS    Secondary Energy|Electricity|Biomass|w/o CCS|3    heat and power co-generation, wood chips
  Biomass CHP CCS                                                  Secondary Energy|Electricity|Biomass|w/ CCS|2     electricity production, at co-generation power plant/wood, post, pipeline 200km, storage 1000m
  Biomass ST                                                       Secondary Energy|Electricity|Biomass|w/o CCS|1    electricity production, at wood burning power plant 20 MW, truck 25km, no CCS
  Biomass IGCC CCS           SE|Electricity|Biomass|IGCCC|w/ CCS   Secondary Energy|Electricity|Biomass|w/ CCS|1     electricity production, from CC plant, 100% SNG, truck 25km, post, pipeline 200km, storage 1000m
  Biomass IGCC               SE|Electricity|Biomass|IGCC|w/o CCS   Secondary Energy|Electricity|Biomass|w/o CCS|2    electricity production, at BIGCC power plant 450MW, no CCS
  Coal PC                    SE|Electricity|Coal|PC|w/o CCS        Secondary Energy|Electricity|Coal|w/o CCS|1       electricity production, hard coal
  Coal IGCC                  SE|Electricity|Coal|IGCC|w/o CCS      Secondary Energy|Electricity|Coal|w/o CCS|2       electricity production, at power plant/hard coal, IGCC, no CCS
  Coal PC CCS                SE|Electricity|Coal|PCC|w/ CCS                                                          electricity production, at power plant/hard coal, post, pipeline 200km, storage 1000m
  Coal IGCC CCS              SE|Electricity|Coal|IGCCC|w/ CCS      Secondary Energy|Electricity|Coal|w/ CCS|1        electricity production, at power plant/hard coal, pre, pipeline 200km, storage 1000m
  Coal CHP                   SE|Electricity|Coal|CHP|w/o CCS       Secondary Energy|Electricity|Coal|w/o CCS|3       heat and power co-generation, hard coal
  Coal CHP CCS                                                     Secondary Energy|Electricity|Coal|w/ CCS|2        electricity production, at co-generation power plant/hard coal, oxy, pipeline
  Gas OC                     SE|Electricity|Gas|GT                 Secondary Energy|Electricity|Gas|w/o CCS|1        electricity production, natural gas, conventional power plant
  Gas CC                     SE|Electricity|Gas|CC|w/o CCS         Secondary Energy|Electricity|Gas|w/o CCS|2        electricity production, natural gas, combined cycle power plant
  Gas CHP                    SE|Electricity|Gas|CHP|w/o CCS        Secondary Energy|Electricity|Gas|w/o CCS|3        heat and power co-generation, natural gas, combined cycle power plant, 400MW electrical
  Gas CHP CCS                                                      Secondary Energy|Electricity|Gas|w/ CCS|2         electricity production, at co-generation power plant/natural gas, post, pipeline
  Gas CC CCS                 SE|Electricity|Gas|w/ CCS             Secondary Energy|Electricity|Gas|w/ CCS|1         electricity production, at power plant/natural gas, pre, pipeline
  Geothermal                 SE|Electricity|Geothermal             Secondary Energy|Electricity|Other                electricity production, deep geothermal
  Hydro                      SE|Electricity|Hydro                  Secondary Energy|Electricity|Hydro                electricity production, hydro, reservoir
  Nuclear                    SE|Electricity|Nuclear                Secondary Energy|Electricity|Nuclear              electricity production, nuclear
  Oil ST                     SE|Electricity|Oil|w/o CCS            Secondary Energy|Electricity|Oil|w/o CCS|1        electricity production, oil
  Oil CC                                                           Secondary Energy|Electricity|Oil|w/o CCS|2        electricity production, oil
  Oil CC CCS                                                       Secondary Energy|Electricity|Oil|w/ CCS|1         electricity production, at co-generation power plant/oil, post, pipeline 200km, storage 1000m
  Oil CHP                                                          Secondary Energy|Electricity|Oil|w/o CCS|3        heat and power co-generation, oil
  Oil CHP CCS                                                      Secondary Energy|Electricity|Oil|w/ CCS|2         electricity production, at co-generation power plant/oil, post, pipeline 200km, storage 1000m
  Solar CSP                  SE|Electricity|Solar|CSP              Secondary Energy|Electricity|Solar|CSP            electricity production, solar thermal parabolic trough, 50 MW
  Solar PV Centralized       SE|Electricity|Solar|PV               Secondary Energy|Electricity|Solar|PV|1           electricity production, photovoltaic, commercial
  Solar PV Residential                                             Secondary Energy|Electricity|Solar|PV|2           electricity production, photovoltaic, residential
  Wind Onshore               SE|Electricity|Wind|Onshore           Secondary Energy|Electricity|Wind|1               electricity production, wind, <1MW turbine, onshore
  Wind Offshore              SE|Electricity|Wind|Offshore          Secondary Energy|Electricity|Wind|2               electricity production, wind, 1-3MW turbine, offshore
  biomass - purpose grown    SE|Electricity|Biomass|Energy Crops   Primary Energy|Biomass|Energy Crops               market for wood chips
  biomass - residual         SE|Electricity|Biomass|Residues       Primary Energy|Biomass|Residues                   Supply of forest residue
 ========================== ===================================== ================================================= ===================================================================================================

.. note::

    IAMs do not necessarily display the same variety of technologies.
    For example, REMIND does not provide a variable for residential PV production.


.. note::

    Because of a lack of more diverse inventories, wind power is only represented
    with relatively small installations (< 1MW, 1-3 MW and >3 MW), in respect to today's
    standard. This can lead to overestimate the associated environmental burden.


The table below shows the correspondence between *premise*, REMIND, IMAGE
and LCI terminology, regarding steel and cement producing technologies. The mapping files are
available in the library root folder: mappingCement_ and mappingSteel_.


 ==================== ====================================== ============================= ==============================
  name in premise      name in REMIND                          name in IMAGE                name in LCI database
 ==================== ====================================== ============================= ==============================
  cement               Production|Industry|Cement             Production|Cement             cement production, Portland
  steel - primary      Production|Industry|Steel|Primary      Production|Steel|Primary      steel production, converter
  steel - secondary    Production|Industry|Steel|Secondary    Production|Steel|Secondary    steel production, electric
 ==================== ====================================== ============================= ==============================

The table below shows the correspondence between *premise*, REMIND, IMAGE
and LCI terminology, regarding fuel producing technologies. The mapping file is
available in the library root folder: mappingFuels_.


 ==================================== =============================================== ========================================================================= ================================================================================================================================================
  name in premise                      name in REMIND                                   name in IMAGE                                                            name in LCI database (only first of several shown)
 ==================================== =============================================== ========================================================================= ================================================================================================================================================
  natural gas                          SE|Gases|Non-Biomass                                                                                                      natural gas, high pressure
  biomethane                           SE|Gases|Biomass                                                                                                          biomethane, gaseous
  diesel                               SE|Liquids|Oil                                  Secondary Energy|Consumption|Liquids|Fossil                               diesel production, low-sulfur
  gasoline                             SE|Liquids|Oil                                  Secondary Energy|Consumption|Liquids|Fossil                               petrol production, low-sulfur
  petrol, synthetic, hydrogen          SE|Liquids|Hydrogen                                                                                                       gasoline production, synthetic, from methanol, hydrogen from electrolysis, CO2 from DAC, energy allocation, at fuelling station
  petrol, synthetic, coal              SE|Liquids|Coal|w/o CCS                                                                                                   gasoline production, synthetic, from methanol, hydrogen from coal gasification, CO2 from DAC, energy allocation, at fuelling station
  diesel, synthetic, hydrogen          SE|Liquids|Hydrogen                                                                                                       diesel production, synthetic, from Fischer Tropsch process, hydrogen from electrolysis, energy allocation, at fuelling station
  diesel, synthetic, coal              SE|Liquids|Coal|w/o CCS                                                                                                   diesel production, synthetic, from Fischer Tropsch process, hydrogen from coal gasification, energy allocation, at fuelling station
  diesel, synthetic, wood              SE|Liquids|Biomass|Biofuel|BioFTR|w/o CCS       Secondary Energy|Consumption|Liquids|Biomass|FT Diesel|Woody|w/oCCS       diesel production, synthetic, from Fischer Tropsch process, hydrogen from wood gasification, energy allocation, at fuelling station
  diesel, synthetic, wood, with CCS    SE|Liquids|Biomass|Biofuel|BioFTRC|w/ CCS       Secondary Energy|Consumption|Liquids|Biomass|FT Diesel|Woody|w/CCS        diesel production, synthetic, from Fischer Tropsch process, hydrogen from wood gasification, with CCS, energy allocation, at fuelling station
  diesel, synthetic, grass                                                             Secondary Energy|Consumption|Liquids|Biomass|FT Diesel|Grassy|w/oCCS      diesel production, synthetic, from Fischer Tropsch process, hydrogen from wood gasification, energy allocation, at fuelling station
  diesel, synthetic, grass, with CCS                                                   Secondary Energy|Consumption|Liquids|Biomass|FT Diesel|Grassy|w/CCS       diesel production, synthetic, from Fischer Tropsch process, hydrogen from wood gasification, with CCS, energy allocation, at fuelling station
  hydrogen, electrolysis               SE|Hydrogen|Electricity                                                                                                   hydrogen supply, from electrolysis
  hydrogen, biomass                    SE|Hydrogen|Biomass|w/o CCS                                                                                               hydrogen supply, from gasification of biomass, by
  hydrogen, biomass, with CCS          SE|Hydrogen|Biomass|w/ CCS                                                                                                hydrogen supply, from gasification of biomass by heatpipe reformer, with CCS
  hydrogen, coal                       SE|Hydrogen|Coal|w/o CCS                                                                                                  hydrogen supply, from coal gasification, by truck, as gaseous, over 500 km
  hydrogen, nat. gas                   SE|Hydrogen|Gas|w/o CCS                                                                                                   hydrogen supply, from SMR of nat. gas, by truck, as gaseous, over 500 km
  hydrogen, nat. gas, with CCS         SE|Hydrogen|Gas|w/ CCS                                                                                                    hydrogen supply, from SMR of nat. gas, with CCS, by truck, as gaseous, over 500 km
  biodiesel, oil                       SE|Liquids|Biomass|Biofuel|Biodiesel|w/o CCS    Secondary Energy|Consumption|Liquids|Biomass|Biodiesel|Oilcrops|w/oCCS    biodiesel production, via transesterification
  biodiesel, oil, with CCS                                                             Secondary Energy|Consumption|Liquids|Biomass|Biodiesel|Oilcrops|w/CCS     biodiesel production, via transesterification
  bioethanol, wood                     SE|Liquids|Biomass|Cellulosic|w/o CCS           Secondary Energy|Consumption|Liquids|Biomass|Ethanol|Woody|w/oCCS         ethanol production, via fermentation, from forest
  bioethanol, wood, with CCS           SE|Liquids|Biomass|Cellulosic|w/ CCS            Secondary Energy|Consumption|Liquids|Biomass|Ethanol|Woody|w/CCS          ethanol production, via fermentation, from forest, with carbon capture and storage
  bioethanol, grass                    SE|Liquids|Biomass|Non-Cellulosic               Secondary Energy|Consumption|Liquids|Biomass|Ethanol|Grassy|w/oCCS        ethanol production, via fermentation, from switchgrass
  bioethanol, grass, with CCS                                                          Secondary Energy|Consumption|Liquids|Biomass|Ethanol|Grassy|w/CCS         ethanol production, via fermentation, from switchgrass, with carbon capture and storage
  bioethanol, grain                    SE|Liquids|Biomass|Conventional Ethanol         Secondary Energy|Consumption|Liquids|Biomass|Ethanol|Maize|w/oCCS         ethanol production, via fermentation, from wheat grains
  bioethanol, grain, with CCS                                                          Secondary Energy|Consumption|Liquids|Biomass|Ethanol|Maize|w/CCS          ethanol production, via fermentation, from corn, with carbon capture and storage
  bioethanol, sugar                    SE|Liquids|Biomass|Conventional Ethanol         Secondary Energy|Consumption|Liquids|Biomass|Ethanol|Sugar|w/oCCS         ethanol production, via fermentation, from sugarbeet
  bioethanol, sugar, with CCS                                                          Secondary Energy|Consumption|Liquids|Biomass|Ethanol|Sugar|w/CCS          ethanol production, via fermentation, from sugarbeet, with carbon capture and storage
  methanol, wood                                                                       Secondary Energy|Consumption|Liquids|Biomass|Methanol|Woody|w/oCCS        market for methanol, from biomass
  methanol, grass                                                                      Secondary Energy|Consumption|Liquids|Biomass|Methanol|Grassy|w/oCCS       market for methanol, from biomass
  methanol, wood, with CCS                                                             Secondary Energy|Consumption|Liquids|Biomass|Methanol|Woody|w/CCS         market for methanol, from biomass
  methanol, grass, with CCS                                                            Secondary Energy|Consumption|Liquids|Biomass|Methanol|Grassy|w/CCS        market for methanol, from biomass
 ==================================== =============================================== ========================================================================= ================================================================================================================================================

.. warning::

    Some fuel types are not properly represented in the LCI database.
    Available inventories for biomass-based methanol production do not differentiate
    between wood and grass as the feedstock.

.. note::

    **Modelling choice**: *premise* builds several potential supply chains for hydrogen.
    Because the logistics to supply hydrogen in the future is not known or indicated by the IAM,
    the choice is made to supply it by truck over 500 km, in a gaseous state.


The production volumes considered for a given scenario can be consulted, like so:

.. code-block:: python

    ndb.scenarios[0]["iam data"].production_volumes


Efficiencies
------------

The efficiency of the different technologies producing
commodities (e.g., electricity, steel, cement, fuel) is modelled to change over time
by the IAM. *premise* stores the relative change in efficiency of such technologies.

The table below shows the correspondence between *premise*, REMIND, IMAGE,
regarding efficiency variables for electricity producing technologies. The mapping file is
available in the library root folder: mappingElec_.

.. _mappingElec: https://github.com/polca/premise/blob/master/premise/data/electricity/electricity_tech_vars.yml

 ================== ================================================== ===========================================
  name in premise    name in REMIND                                      name in IMAGE
 ================== ================================================== ===========================================
  Biomass CHP        Tech|Electricity|Biomass|CHP|w/o CCS|Efficiency    Efficiency|Electricity|Biomass|w/o CCS|3
  Biomass CHP CCS                                                       Efficiency|Electricity|Biomass|w/ CCS|2
  Biomass ST                                                            Efficiency|Electricity|Biomass|w/o CCS|1
  Biomass IGCC CCS   Tech|Electricity|Biomass|IGCCC|w/ CCS|Efficiency   Efficiency|Electricity|Biomass|w/ CCS|1
  Biomass IGCC       Tech|Electricity|Biomass|IGCC|w/o CCS|Efficiency   Efficiency|Electricity|Biomass|w/o CCS|2
  Coal PC            Tech|Electricity|Coal|PC|w/o CCS|Efficiency        Efficiency|Electricity|Coal|w/o CCS|1
  Coal IGCC          Tech|Electricity|Coal|IGCC|w/o CCS|Efficiency      Efficiency|Electricity|Coal|w/o CCS|2
  Coal PC CCS        Tech|Electricity|Coal|PCC|w/ CCS|Efficiency
  Coal IGCC CCS      Tech|Electricity|Coal|IGCCC|w/ CCS|Efficiency      Efficiency|Electricity|Coal|w/ CCS|1
  Coal CHP           Tech|Electricity|Coal|CHP|w/o CCS|Efficiency       Efficiency|Electricity|Coal|w/o CCS|3
  Coal CHP CCS                                                          Efficiency|Electricity|Coal|w/ CCS|2
  Gas OC             Tech|Electricity|Gas|GT|Efficiency                 Efficiency|Electricity|Gas|w/o CCS|1
  Gas CC             Tech|Electricity|Gas|CC|w/o CCS|Efficiency         Efficiency|Electricity|Gas|w/o CCS|2
  Gas CHP            Tech|Electricity|Gas|CHP|w/o CCS|Efficiency        Efficiency|Electricity|Gas|w/o CCS|3
  Gas CHP CCS                                                           Efficiency|Electricity|Gas|w/ CCS|2
  Gas CC CCS         Tech|Electricity|Gas|CCC|w/ CCS|Efficiency         Efficiency|Electricity|Gas|w/ CCS|1
  Nuclear                                                               Efficiency|Electricity|Nuclear
  Oil ST             Tech|Electricity|Oil|DOT|Efficiency                Efficiency|Electricity|Oil|w/o CCS|1
  Oil CC                                                                Efficiency|Electricity|Oil|w/o CCS|2
  Oil CC CCS                                                            Efficiency|Electricity|Oil|w/ CCS|1
  Oil CHP                                                               Efficiency|Electricity|Oil|w/o CCS|3
  Oil CHP CCS                                                           Efficiency|Electricity|Oil|w/ CCS|2
 ================== ================================================== ===========================================

The table below shows the correspondence between *premise*, REMIND, IMAGE,
regarding efficiency variables for cement and steel
producing technologies. For cement and steel, it is different, as *premise*
derives efficiencies by dividing the the final energy demand by the production volume
(to obtain GJ/t steel or cement). This is because efficiency variables for cement
and steel is not always given as such. The mapping files are
available in the library root folder: mappingCement_ and mappingSteel_.

.. _mappingCement: https://github.com/polca/premise/blob/master/premise/data/cement/cement_tech_vars.yml
.. _mappingSteel: https://github.com/polca/premise/blob/master/premise/data/steel/steel_tech_vars.yml

 ==================== ========================================== ==============================
  name in premise      name in REMIND                              name in IMAGE
 ==================== ========================================== ==============================
  cement               Final Energy|Industry|Cement               FE|Industry|Cement
  steel - primary      Final Energy|Industry|Steel                FE|Industry|Steel|Primary
  steel - secondary    Final Energy|Industry|Steel|Electricity    FE|Industry|Steel|Secondary
 ==================== ========================================== ==============================

The table below shows the correspondence between *premise*, REMIND, IMAGE,
regarding efficiency variables for fuels producing technologies. The mapping file is
available in the library root folder: mappingFuels_.

.. _mappingFuels: https://github.com/polca/premise/blob/master/premise/data/fuels/fuel_tech_vars.yml

 ==================================== ======================================================================= ========================================================
  name in premise                      name in REMIND                                                           name in IMAGE
 ==================================== ======================================================================= ========================================================
  biomethane                           Tech|Gases|Biomass|w/o CCS|Efficiency
  diesel                               Tech|Liquids|Oil|Efficiency
  gasoline                             Tech|Liquids|Oil|Efficiency
  diesel, synthetic, wood                                                                                      Efficiency|Liquids|Biomass|FT Diesel|Woody|w/o CCS
  diesel, synthetic, wood, with CCS                                                                            Efficiency|Liquids|Biomass|FT Diesel|Woody|w/ CCS
  diesel, synthetic, grass                                                                                     Efficiency|Liquids|Biomass|FT Diesel|Woody|w/o CCS
  diesel, synthetic, grass, with CCS                                                                           Efficiency|Liquids|Biomass|FT Diesel|Woody|w/ CCS
  biodiesel, oil                       Tech|Liquids|Biomass|Biofuel|Biodiesel|w/o CCS|Efficiency               Efficiency|Liquids|Biomass|Biodiesel|Oilcrops|w/o CCS
  biodiesel, oil, with CCS                                                                                     Efficiency|Liquids|Biomass|Biodiesel|Oilcrops|w/ CCS
  bioethanol, wood                     Tech|Liquids|Biomass|Biofuel|Ethanol|Cellulosic|w/o CCS|Efficiency      Efficiency|Liquids|Biomass|Ethanol|Woody|w/o CCS
  bioethanol, wood, with CCS                                                                                   Efficiency|Liquids|Biomass|Ethanol|Woody|w/ CCS
  bioethanol, grass                    Tech|Liquids|Biomass|Biofuel|Ethanol|Cellulosic|w/o CCS|Efficiency      Efficiency|Liquids|Biomass|Ethanol|Grassy|w/o CCS
  bioethanol, grass, with CCS                                                                                  Efficiency|Liquids|Biomass|Ethanol|Grassy|w/ CCS
  bioethanol, grain                    Tech|Liquids|Biomass|Biofuel|Ethanol|Conventional|w/o CCS|Efficiency    Efficiency|Liquids|Biomass|Ethanol|Maize|w/o CCS
  bioethanol, grain, with CCS                                                                                  Efficiency|Liquids|Biomass|Ethanol|Maize|w/ CCS
  bioethanol, sugar                    Tech|Liquids|Biomass|Biofuel|Ethanol|Conventional|w/o CCS|Efficiency    Efficiency|Liquids|Biomass|Ethanol|Sugar|w/o CCS
  bioethanol, sugar, with CCS                                                                                  Efficiency|Liquids|Biomass|Ethanol|Sugar|w/ CCS
  methanol, wood                                                                                               Efficiency|Liquids|Biomass|Methanol|Woody|w/o CCS
  methanol, grass                                                                                              Efficiency|Liquids|Biomass|Methanol|Grassy|w/o CCS
  methanol, wood, with CCS                                                                                     Efficiency|Liquids|Biomass|Methanol|Woody|w/ CCS
  methanol, grass, with CCS                                                                                    Efficiency|Liquids|Biomass|Methanol|Grassy|w/ CCS
 ==================================== ======================================================================= ========================================================


*premise* stores the change in efficiency (called *scaling factor*) of a given technology
relative to 2020. This is based on the fact that the efficiency of ecoinvent datasets
are believed to reflect current (2020) efficiency.

.. note::

    If a technology, in a given region, is given a *scaling factor* of 1.2 in 2030,
    this means that the corresponding ecoinvent dataset is adjusted so that its
    efficiency is improved by 20% (by multiplying the dataset inputs by 1/1.2).
    In other words, *premise* does not use the efficiency given by the IAM,
    but rather its change over time relative to 2020.

The *scaling factors* considered for a given scenario can be consulted, like so:

.. code-block:: python

    ndb.scenarios[0]["iam data"].efficiency

Land use and land use change
----------------------------

When building prospective databases using the IAM IMAGE model, the latter provides
additional variables relating to average *land use* and *land use change* emissions, for each type of
crop grown to be used in biofuel production.
Upon the creation of biofuel supply chains in the *Fuels* transformation function, such information
is used to adjust the inventories of crop farming datasets. The table below shows the IMAGE variables
used to that effect. The mapping file is
available in the library root folder: mappingCrops_.

.. _mappingCrops: https://github.com/polca/premise/blob/master/premise/data/fuels/crops_properties.yml

 ========================= ========================== ========================================== =============================================================
  Crop family in premise    Crop type in premise       Land use variable in IMAGE [Ha/GJ-Prim]    Land use change variable in IMAGE [kg CO2/GJ-Prim]
 ========================= ========================== ========================================== =============================================================
  sugar                     sugarbeet, sugarcane       Land Use|Average|Biomass|Sugar             Emission Factor|CO2|Energy|Supply|Biomass|Average|Sugar
  oil                       rapeseed, palm oil         Land Use|Average|Biomass|OilCrop           Emission Factor|CO2|Energy|Supply|Biomass|Average|Oilcrops
  wood                      poplar, eucalyptus         Land Use|Average|Biomass|Woody             Emission Factor|CO2|Energy|Supply|Biomass|Average|Woody
  grass                     switchgrass, miscanthus    Land Use|Average|Biomass|Grassy            Emission Factor|CO2|Energy|Supply|Biomass|Average|Grassy
  grain                     corn                       Land Use|Average|Biomass|Maize             Emission Factor|CO2|Energy|Supply|Biomass|Average|Maize
 ========================= ========================== ========================================== =============================================================

The *land use* and *land use change* emissions considered for a given scenario
can be consulted, like so:

.. code-block:: python

    ndb.scenarios[0]["iam data"].land_use
    ndb.scenarios[0]["iam data"].land_use_change

Carbon Capture and Storage
--------------------------

Some scenarios involve the capture and storage of CO2 emissions
of certain sectors (e.g., cement and steel).
The capture rate of a given sector is calculated
from the IAM data file, as::

    rate = amount of CO2 captured / (amount of CO2 captured + amount of CO2 not captured)

The table below lists the variables needed to calculate those rates.

 ============================== =============================== ============================================
  name in premise                name in REMIND                  name in IMAGE
 ============================== =============================== ============================================
  cement - CO2 (not captured)    Emi|CO2|FFaI|Industry|Cement    Emissions|CO2|Industry|Cement|Gross
  cement - CCO2 (captured)       Emi|CCO2|FFaI|Industry|Cement   Emissions|CO2|Industry|Cement|Sequestered
  steel - CO2 (not captured)     Emi|CO2|FFaI|Industry|Steel     Emissions|CO2|Industry|Steel|Gross
  steel - CCO2 (captured)        Emi|CCO2|FFaI|Industry|Steel    Emissions|CO2|Industry|Steel|Sequestered
 ============================== =============================== ============================================


The *carbon capture rates* which are floating values
comprised between 0 and 1, can be consulted like so:

.. code-block:: python

    ndb.scenarios[0]["iam data"].carbon_capture_rate


Data sources external to the IAM
--------------------------------

*premise* tries to adhere to the IAM scenario data as much as possible. There are
however a number of cases where external data sources are used. This is notably the case
for non-CO2 pollutants emissions for different sectors (electricity, steel and cement),
fuel mixes and power generation for the cement industry, as well as expected efficiency gains
for photovoltaic panels.

Air emissions
*************

*premise* relies on projections from the air emissions models GAINS-EU_ and GAINS-IAM_
to adjust the emissions of pollutants for different sectors.
As with efficiencies, *premise* stores the change in emissions (called *scaling factor*)
of a given technology relative to 2020. This is based on the fact that the emissions of
ecoinvent datasets are believed to reflect the current (2020) situation.
Hence, if a technology, in a given region, has a *scaling factor* of 1.2 in 2030,
this means that the corresponding ecoinvent dataset is adjusted so that its emissions
of a given substance is improved by 20%. In other words, *premise* does not use
the emissions level given by GAINS, but rather its change over time relative to 2020.

For more information about this step, refer to sub-section "GAINS emission factors" in the
EXTRACT section.

.. _GAINS-EU: https://gains.iiasa.ac.at/gains/EUN/index.login
.. _GAINS-IAM: https://gains.iiasa.ac.at/gains/IAM/index.login


Cement production
*****************

A number of parameters to model future clinker/cement production is sourced from the
IAM file, such as:

* The expected change in fuel efficiency for clinker production.


Photovoltaic panels
*******************

Module efficiencies in 2010 for micro-Si and single-Si are from IEA_ Task 12
report. For multi-Si, CIGS, CIS and CdTe, they are from IEA2_ road map report
on PV panels.

.. _IEA2: https://iea.blob.core.windows.net/assets/3a99654f-ffff-469f-b83c-bf0386ed8537/pv_roadmap.pdf

Current (2020) module efficiencies for all PV types are given by a 2021 report
from the Fraunhofer_ Institute.

The efficiencies indicated for 2050 are what has been obtained in laboratory
conditions by the Fraunhofer_ Institute. In other words, it is assumed that by 2050,
solar PVs will reach production level efficiencies equal to those observed today
in laboratories.

.. _Fraunhofer: https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/Photovoltaics-Report.pdff

 ====================== =========== ============ =========== ======= ====== =======
  % module efficiency    micro-Si    single-Si    multi-Si    CIGS    CIS    CdTe
 ====================== =========== ============ =========== ======= ====== =======
  2010                   10          15.1         14          11      11     10
  2020                   11.9        17.9         16.8        14      14     16.8
  2050                   12.5        26.7         24.4        23.4    23.4   21
 ====================== =========== ============ =========== ======= ====== =======




