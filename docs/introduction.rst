In a nutshell
"""""""""""""

Purpose
-------

*premise* enables the alignment of life cycle inventories within the ecoinvent_
3.6-3.10 database, using either a "cut-off" or "consequential"
system model, to match the output results of Integrated
Assessment Models (IAMs) such as REMIND_, IMAGE_ or TIAM-UCL_.
This allows for the creation of life cycle inventory databases
under future policy scenarios for any year between 2005 and 2100.

.. _ecoinvent: https://ecoinvent.org/
.. _REMIND: https://www.pik-potsdam.de/en/institute/departments/transformation-pathways/models/remind
.. _IMAGE: https://models.pbl.nl/image/index.php/Welcome_to_IMAGE_3.2_Documentation
.. _TIAM-UCL: https://www.ucl.ac.uk/energy-models/models/tiam-ucl


.. note::

    The ecoinvent database is not included in this package. You need to have a valid license for ecoinvent 3.6-3.9.1 to use *premise*.
    Also, please read carefully ecoinvent's EULA_ before using *premise*.

.. _EULA: https://ecoinvent.org/app/uploads/2024/01/EULA_new_branding_08_11_2023.pdf

Publication
-----------

The methodology behind *premise* is described in the following publication:

R. Sacchi, T. Terlouw, K. Siala, A. Dirnaichner, C. Bauer, B. Cox, C. Mutel, V. Daioglou, G. Luderer,
PRospective EnvironMental Impact asSEment (premise): A streamlined approach to producing databases for prospective life cycle assessment using integrated assessment models,
Renewable and Sustainable Energy Reviews, 2022, https://doi.org/10.1016/j.rser.2022.112311.

.. note::

    If you use *premise* in your research, please cite the above publication.

Additionally, you may want to cite the ecoinvent database:

Wernet, G. et al. The ecoinvent database version 3 (part I): overview and methodology. Int. J. Life Cycle Assess. 21, 1218–1230 (2016) . http://link.springer.com/10.1007/s11367-016-1087-8.

Finally, you may want to cite the IAM model used with *premise*:

* REMIND: Baumstark et al. REMIND2.1: transformation and innovation dynamics of the energy-economic system within climate and sustainability limits, Geoscientific Model Development, 2021.
* IMAGE: Stehfest, Elke, et al. Integrated assessment of global environmental change with IMAGE 3.0: Model description and policy applications. Netherlands Environmental Assessment Agency (PBL), 2014.
* TIAM-UCL: L. Clarke, et al. International climate policy architectures: Overview of the EMF 22 International Scenarios, Energy Economics, 2009.


Models
------

+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Model       | Description                                                                                                                                                                                                                                               |
+=============+===========================================================================================================================================================================================================================================================+
| REMIND      | REMIND (Regionalized Model of Investment and Development) is an integrated assessment model that combines macroeconomic growth, energy system, and climate policy analysis. It is designed to analyze long-term energy transition pathways, accounting for  |
|             | technological, economic, and environmental factors. REMIND simulates how regions invest in different technologies and energy resources to balance economic growth and climate targets, while considering factors like energy efficiency, emissions, and      |
|             | resource availability. The model is particularly strong in its detailed representation of energy markets and macroeconomic interactions across regions, making it valuable for global climate policy assessments.                                            |
+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| IMAGE       | IMAGE (Integrated Model to Assess the Global Environment) is a comprehensive IAM developed to explore the interactions between human development, energy consumption, and environmental systems over the long term. It focuses on assessing how land use,   |
|             | food systems, energy systems, and climate change interact under different policy scenarios. The model integrates biophysical processes, such as land-use change and greenhouse gas emissions, with socio-economic drivers like population growth and economic |
|             | development. IMAGE is commonly used for analyzing sustainable development strategies, climate impacts, biodiversity loss, and exploring mitigation and adaptation options.                                                                                   |
+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TIAM-UCL    | TIAM-UCL (TIMES Integrated Assessment Model by University College London) is a global energy system model based on the TIMES (The Integrated MARKAL-EFOM System) framework, developed to evaluate long-term decarbonization pathways for global energy       |
|             | systems. It provides detailed insights into energy technology options, resource availability, and emission reduction strategies under various climate policy scenarios. The model focuses on the trade-offs and synergies between energy security, economic    |
|             | costs, and environmental outcomes. TIAM-UCL is frequently used to analyze scenarios consistent with the Paris Agreement and examine technological innovation's role in mitigating climate change globally.                                                    |
+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Workflow
--------

.. image:: main_workflow.png

As illustrated in the workflow diagram above, *premise* follows an Extract, Transform, Load (ETL_) process:

Extract the ecoinvent database from a Brightway_ project or from ecospold2_ files.
Expand the database by adding additional inventories for future production pathways for certain commodities, such as electricity, steel, cement, etc.
Modify the ecoinvent database, focusing primarily on process efficiency improvements and market adjustments.
Load the updated database back into a Brightway project or export it as a set of CSV files, such as Simapro CSV files.

.. _brightway: https://brightway.dev/
.. _ecospold2: https://ecoinvent.org/the-ecoinvent-database/data-formats/ecospold2/
.. _ETL: https://www.guru99.com/etl-extract-load-process.html#:~:text=ETL%20is%20a%20process%20that,is%20Extract%2C%20Transform%20and%20Load.

Default IAM scenarios
---------------------

Provided a decryption key (ask the maintainers_), the following IAM scenarios are available when installing *premise*:

+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP/RCP scenario | GMST increase by 2100 | Society/economy trend                                                              | Climate policy                              | REMIND          | IMAGE      | TIAM-UCL   |
+==================+=======================+====================================================================================+=============================================+=================+============+-===========+
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP1-None        | 2.3-2.8 °C            | Optimistic trends for human develop. and economy, driven by sustainable practices. | None                                        | SSP1-Base       | SSP1-Base  |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP1-None        | ~2.2 °C               | Optimistic trends for human develop. and economy, driven by sustainable practices. | National Policies Implemented (NPI).        | SSP1-NPi        |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP1-None        | ~1.9 °C               | Optimistic trends for human develop. and economy, driven by sustainable practices. | Nationally Determined Contributions (NDCs). | SSP1-NDC        |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP1-RCP2.6      | ~1.7 °C               | Optimistic trends for human develop. and economy, driven by sustainable practices. | Paris Agreement objective.                  | SSP1-PkBudg1150 |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP1-RCP1.9      | ~1.3 °C               | Optimistic trends for human develop. and economy, driven by sustainable practices. | Paris Agreement objective.                  | SSP1-PkBudg500  |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP2-None        | ~3.5 °C               | Extrapolation from historical developments.                                        | None (eq. to RCP6)                          | SSP2-Base       | SSP2-Base  | SSP2-Base  |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP2-None        | ~3.3 °C               | Extrapolation from historical developments.                                        | National Policies Implemented (NPI).        | SSP2-NPi        |            | SSP2-RCP45 |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP2-None        | ~2.5 °C               | Extrapolation from historical developments.                                        | Nationally Determined Contributions (NDCs). | SSP2-NDC        |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP2-RCP2.6      | 1.6-1.8 °C            | Extrapolation from historical developments.                                        | Paris Agreement objective.                  | SSP2-PkBudg1150 | SSP2-RCP26 | SSP2-RCP26 |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP2-RCP1.9      | 1.2-1.4 °C            | Extrapolation from historical developments.                                        | Paris Agreement objective.                  | SSP2-PkBudg500  | SSP2-RCP19 | SSP2-RCP19 |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP5-None        | ~4.5 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | None                                        | SSP5-Base       |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP5-None        | ~4.0 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | National Policies Implemented (NPI).        | SSP5-NPi        |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP5-None        | ~3.0 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | Nationally Determined Contributions (NDCs). | SSP5-NDC        |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP5-RCP2.6      | ~1.7 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | Paris Agreement objective.                  | SSP5-PkBudg1150 |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+
| SSP5-RCP1.9      | ~1.0 °C               | Optimistic trends for human develop. and economy, driven by fossil fuels.          | Paris Agreement objective.                  | SSP5-PkBudg500  |            |            |
+------------------+-----------------------+------------------------------------------------------------------------------------+---------------------------------------------+-----------------+------------+------------+

CarbonBrief_ wrote a good article explaining the meaning of the SSP/RCP system.

Additionally, we provided a summary of the main characteristics of each scenario `here <https://premisedash-6f5a0259c487.herokuapp.com/>`_.


.. _CarbonBrief: https://www.carbonbrief.org/explainer-how-shared-socioeconomic-pathways-explore-future-climate-change

You can however use any other scenario files generated by REMIND or IMAGE.
If you wish to use an IAM file which has not been generated by either of these
two models, you should refer to the **Mapping** section.

.. _maintainers: mailto:romain.sacchi@psi.ch


Requirements
------------
* Python language interpreter **>=3.9**
* License for ecoinvent 3
* Brightway2 (optional)

.. note::

    If you wish to export Brightway 2.5-compatible databases, you will need ot upgrade `bw2data` to >= 4.0.0.

How to install this package?
----------------------------

Two options:

A development version with the latest advancements (but with the risks of unseen bugs),
is available on Anaconda Cloud:

.. code-block:: python

    conda install -c conda-forge premise

For a more stable and proven version, from Pypi:

.. code-block:: python

    pip install premise

This will install the package and the required dependencies.

How to use it?
--------------

Examples notebook
*****************

`This notebook <https://github.com/polca/premise/blob/master/examples/examples.ipynb>`_ will show
you everything you need to know to use *premise*.

ScenarioLink plugin
*******************
There now exists a plugin for Activity Browser, called ScenarioLink, which allows you to
directly download IAM scenario-based premise databases from the browser, without the use of premise.
You can find it `here <https://github.com/polca/ScenarioLink>`_.

Active contributors
-------------------

* `Romain Sacchi <https://github.com/romainsacchi>`_
* `Alvaro Hahn Menacho <https://github.com/alvarojhahn>`_

Historical contributors
-----------------------

* `Alois Dirnaichner <https://github.com/Loisel>`_
* `Chris Mutel <https://github.com/cmutel>`_
* `Brian Cox <https://github.com/brianlcox>`_
