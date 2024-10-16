# ``premise``

<div style="text-align:center">
<img src="https://github.com/polca/premise/raw/master/docs/large.png" height="300"/>
</div>

# **PR**ospective **E**nviron**M**ental **I**mpact As**SE**ssment
## Coupling the ecoinvent database with projections from Integrated Assessment Models (IAM)


<p align="center">
  <a href="https://badge.fury.io/py/premise" target="_blank"><img src="https://badge.fury.io/py/premise.svg"></a>
  <a href="https://anaconda.org/conda-forge/premise" target="_blank"><img src="https://img.shields.io/conda/vn/conda-forge/premise.svg"></a>
  <a href="https://github.com/polca/premise" target="_blank"><img src="https://github.com/polca/premise/actions/workflows/main.yml/badge.svg?branch=master"></a>
  <a href="https://coveralls.io/github/polca/premise" target="_blank"><img src="https://coveralls.io/repos/github/polca/premise/badge.svg"></a>
  <a href="https://premise.readthedocs.io/en/latest/" target="_blank"><img src="https://readthedocs.org/projects/premise/badge/?version=latest"></a>
</p>

``premise`` is a Python tool for prospective life cycle assessment. 
It allows users to project the ecoinvent 3 database into the future, 
using scenarios from Integrated Assessment Models (IAMs). It does so by 
modifying the ecoinvent database to reflect projected energy policy trajectories, include emerging
technologies, modify market shares as well as technologies' efficiency.

Among others, it can be used to assess the environmental impacts of future energy systems,
and to compare different energy policies. It includes a set of IAM scenarios
and a set of tools to create custom scenarios.

The tool was designed to be user-friendly and to allow for reproducible results. 
While it is built on the [brightway2 framework](https://docs.brightway.dev/en/latest/), 
its outputs can naturally be used in [Activity Browser](https://github.com/LCA-ActivityBrowser/activity-browser), 
but also in other LCA software, such as SimaPro, OpenLCA, or directly in Python.

The tool is described in the following scientific publication: [Sacchi et al, 2022](https://doi.org/10.1016/j.rser.2022.112311).
If this tool helps you in your research, please consider citing this publication.

Also, use the following references to cite the scenarios used with the tool:

- REMIND scenarios: Baumstark et al. REMIND2.1: transformation and innovation dynamics of the energy-economic system within climate and sustainability limits, Geoscientific Model Development, 2021.
- IMAGE scenarios: Stehfest, Elke, et al. Integrated assessment of global environmental change with IMAGE 3.0: Model description and policy applications. Netherlands Environmental Assessment Agency (PBL), 2014.
- TIAM-UCL scenarios: []()

Models
------

The tool currently supports the following IAMs:

| Model    | Description |
|----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| REMIND   | REMIND (Regionalized Model of Investment and Development) is an integrated assessment model that combines macroeconomic growth, energy system, and climate policy analysis. It is designed to analyze long-term energy transition pathways, accounting for technological, economic, and environmental factors. REMIND simulates how regions invest in different technologies and energy resources to balance economic growth and climate targets, while considering factors like energy efficiency, emissions, and resource availability. The model is particularly strong in its detailed representation of energy markets and macroeconomic interactions across regions, making it valuable for global climate policy assessments. |
| IMAGE    | IMAGE (Integrated Model to Assess the Global Environment) is a comprehensive IAM developed to explore the interactions between human development, energy consumption, and environmental systems over the long term. It focuses on assessing how land use, food systems, energy systems, and climate change interact under different policy scenarios. The model integrates biophysical processes, such as land-use change and greenhouse gas emissions, with socio-economic drivers like population growth and economic development. IMAGE is commonly used for analyzing sustainable development strategies, climate impacts, biodiversity loss, and exploring mitigation and adaptation options. |
| TIAM-UCL | TIAM-UCL (TIMES Integrated Assessment Model by University College London) is a global energy system model based on the TIMES (The Integrated MARKAL-EFOM System) framework, developed to evaluate long-term decarbonization pathways for global energy systems. It provides detailed insights into energy technology options, resource availability, and emission reduction strategies under various climate policy scenarios. The model focuses on the trade-offs and synergies between energy security, economic costs, and environmental outcomes. TIAM-UCL is frequently used to analyze scenarios consistent with the Paris Agreement and examine technological innovation's role in mitigating climate change globally. |



What's new in 2.1.4?
====================

- Integration of TIAM-UCL scenarios


What's new in 2.1.1?
====================

- Ecoinvent 3.10 support
- Adds inventories on ammonia production
- Fixes issues with scaling applied to PV systems

What's new in 2.1.0?
====================

- More efficient use of memory.
- Easier syntax for using custom scenarios.

What's new in 2.0.0?
====================

- Adds .write_db_to_olca(), which produces a slightly modified version of a Simapro CSV database file which can then be imported in OpenLCA. The use of the SimaPro_Import.csv mapping file must be selected.
- Marginal mixes were wrong because the average lifetime of the mix was calculated using the technology-specific lead time values instead of lifetime values.
- Fix issue with CCS implementation in IMAGE scenarios
- Fix several issues with external/custom scenarios linking algorithm.
- Drops DAC efficiency improvement based on cumulated deployment. Uses directly efficiency variables now (if provided).
- Improves documentation for consequential modelling.
- Code-breaking change: update functions are now called like so: .update(xxx). For example, to update the electricity sector: .update(["electricity",]). To update all sectors: .update().
- Changes minimum Python requirement to 3.10

What's new in 1.8.0?
====================

-  Added support for brightway 2.5 (requires `bw2data` >= 4.0.0)
-  Added support for Python 3.11
-  Uses bw2io 0.8.10
-  Adds electricity storage in electricity markets -- see [docs](https://premise.readthedocs.io/en/latest/transform.html#storage)
-  Adds [scenario explorer dashboard](https://premisedash-6f5a0259c487.herokuapp.com/)

What's new in 1.5.0?
====================

-  Added support for ecoinvent 3.9 and 3.9.1
-  Added support for ecoinvent 3.8 and 3.9/3.9.1 consequential -- see [docs](https://premise.readthedocs.io/en/latest/consequential.html)
-  Added REMIND SSP1 and SSP5 scenarios -- see [docs](https://premise.readthedocs.io/en/latest/introduction.html#default-iam-scenarios)
-  Updated GAINS emission factors, using GAINS-EU and GAINS-IAM -- see [docs](https://premise.readthedocs.io/en/latest/transform.html#gains-emission-factors)
-  Added new inventories for DAC and DACCS -- see [docs](https://premise.readthedocs.io/en/latest/transform.html#direct-air-capture)
-  Added new inventories for EPR and SMR nuclear reactors -- see [EPR inventories](https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-nuclear_EPR.xlsx) and [SMR inventories](https://github.com/polca/premise/blob/master/premise/data/additional_inventories/lci-nuclear_SMR.xlsx)
-  Made mapping to new IAM models easier -- see [docs](https://premise.readthedocs.io/en/latest/mapping.html)
-  Better logging of changes made to the ecoinvent database -- see [docs](https://premise.readthedocs.io/en/latest/transform.html#logs)

What's new in 1.3.0?
====================

-   Added support for user-generated scenarios (see [docs](https://premise.readthedocs.io/en/latest/user_scenarios.html) and [notebook](https://github.com/polca/premise/blob/master/examples/examples%20user-defined%20scenarios.ipynb))
-   Updated REMIND scenarios to REMIND v.3.0



Documentation
-------------
[https://premise.readthedocs.io/en/latest/](https://premise.readthedocs.io/en/latest/)

Objective
---------

The objective is to produce life cycle inventories under future energy policies, by modifying the inventory database
ecoinvent 3 to reflect projected energy policy trajectories.

Requirements
------------
* **Python 3.10 or 3.11**
* License for [ecoinvent 3][1]. Please note that the ecoinvent database is not included in this package. Also, read ecoinvent's [GDPR & EULA](https://ecoinvent.org/gdpr-eula/).
* Some IAM output files come with the library and are located by default in the subdirectory "/data/iam_output_files". **If you wish to use
 those files, you need to request (by [email](mailto:romain.sacchi@psi.ch)) an encryption key from the developers**.
 A file path can be specified to fetch IAM output files elsewhere on your computer.
 * [brightway2][2] (optional). If you want to use the results in the Brightway 2 framework 8and Activity Browser), you need `bw2data <4.0.0`. To produce Brightway 2.5-compatible databases, you need `bw2data >=4.0.0`.

How to install this package?
----------------------------

Two options:

From Pypi:

    pip install premise

will install the package and the required dependencies.

A development version with the latest advancements (but with the risks of unseen bugs),
is available from Anaconda Cloud:

    conda install -c conda-forge premise


How to use it?
--------------

The best way is to follow [the examples from the Jupyter Notebook](https://github.com/polca/premise/blob/master/examples/examples.ipynb). 

# Support

Do not hesitate to contact [romain.sacchi@psi.ch](mailto:romain.sacchi@psi.ch).

## Contributors

* [Romain Sacchi](https://github.com/romainsacchi)
* [Alois Dirnaichner](https://github.com/Loisel)
* [Tom Mike Terlouw](https://github.com/tomterlouw)
* [Laurent Vandepaer](https://github.com/lvandepaer)
* [Chris Mutel](https://github.com/cmutel/)


## Maintainers

* [Romain Sacchi](https://github.com/romainsacchi)
* [Chris Mutel](https://github.com/cmutel/)

## Contributing

See [contributing](https://github.com/polca/premise/blob/master/CONTRIBUTING.md).

## References

[1]:https://www.ecoinvent.org/
[2]:https://brightway.dev/

## License

[BSD-3-Clause](https://github.com/polca/premise/blob/master/LICENSE).
Copyright 2020 Potsdam Institute for Climate Impact Research, Paul Scherrer Institut.
