Frequently Asked Questions
""""""""""""""""""""""""""

Here are some frequently asked questions about ``premise``.
If you have a question that is not answered here, please contact us.

IAM models
----------

I use a different IAM than REMIND or IMAGE ... Can I still use ``premise``?
___________________________________________________________________________

There is a MAPPING section in the documentation
that explains how to link to a new IAM. The YAML files under ````premise``/iam_variables_mapping``
are the main body of files that needs to
be changed, to properly establish a correspondence between your IAM variables
and the variables used in ``premise``. It is also necessary to provide ``premise``
with the geographical definitions of the regions used in your IAM. This is done
by providing a .json file with the regions and their corresponding ecoinvent regions.
The rest of the code is generic and should work with any IAM.

What columns are necessary in the IAM files?
____________________________________________

The code has been refactored since.
Any column other than:

* Region
* Variable
* Unit
* and the variable values for each time step

is ignored.

How big an effort would it be to link to a new IAM? As simple as an extension of the mapping files? What difficulties can be anticipated?
_________________________________________________________________________________________________________________________________________

In principle, it is easy. Linking to a new IAM model is a matter of:

* providing the IAM variable for each ``premise`` variable listed in the .yaml mapping files
* and the geographical definitions of the regions used in the IAM.

In practice, it may not always be that simple.
The IAM variables are not always available in the IAM output files (e.g., efficiency or land use-related variables).
In that case, they need to be calculated from other variables or skipped.
Also, some IAM models may represent a technology not yet considered in ``premise`` (e.g., nuclear fusion).
In some cases, ``premise``'s code needs to be extended.

IAM data collection
-------------------

How was the list of variables in the mapping files established?
_______________________________________________________________

The list of IAM variables and mapping with ``premise`` variables has been established
through collaboration with developers of IAM models, to ensure that the meaning between
each IAM variable corresponds with that of ``premise``.

Is it possible to expand this list? (e.g. agriculture crops for energy)
_______________________________________________________________________

It is certainly possible to extend this list. You would however need to extend
``premise``'s code to tell it what to do with these additional variables. For example, if you want to
use the IAM output for integrating projections that relate to agriculture crops for energy,
you would need to write a module in ``premise`` (e.g., energy_crops.py) that would perform a series
of modifications on the LCA datasets, just like other modules do.

Is the unit and the description of these parameters documented? Or are they necessarily the same as the ones of the ecoinvent datasets they refer to?
_____________________________________________________________________________________________________________________________________________________

They are now documented, under the MAPPING section.
There are essentially two types of variables:

* variables that relate to production volumes of technologies, which units must represent a production volume over time (e.g., GWh/year)
* variables that relate to the efficiency of technologies, which is unitless, or represented by an efficiency ratio (e.g., %)

What if a variable in ``premise`` corresponds to several variables in the IAM?
______________________________________________________________________________

We have not really seen that case yet. In any case, mapping one IAM variable
to two ``premise`` variables is possible (whether it is methodologically correct
is a question left to your appreciation).

Regionalization
---------------

Are datasets regionalized on the basis of the IAM scenario only, or does it come from other sources?
____________________________________________________________________________________________________

``premise`` tries to limit the use of external sources of data.
At the moment, the only sources of data, other than those from the IAM scenario, used for projections are:

- efficiency values for different photovoltaic panels (taken from the Fraunhofer ISE database)
- emissions factors for local air pollution (taken from the GAINS-EU and GAINS-IAM databases)

Hence, the regionalization of datasets is based on the IAM scenario only.

Does ``premise`` generate more regionalised datasets than in original EI3.x database?
_________________________________________________________________________________

Yes. ``premise`` generates regionalized datasets for all regions in the IAM model, for
each technology for which a IAM-to-``premise`` correspondence is provided, if not already existing in the Ecoinvent database.
For example, if the IAM model
considers technology A over 10 regions, ``premise`` collects datasets in the ecoinvent database
(or imported inventories) that represent technology A and duplicates it for each region. Sometimes,
only one dataset is available in the ecoinvent database, in which case ``premise`` duplicates it 10 times.
Other times, several datasets are available (ie.g., in FR, CN and RoW), in which case ``premise`` uses the French
dataset for the European region, the Chinese dataset for the Chinese region, and the RoW dataset for the other IAM regions.
Then, ``premise`` proceeds to regionalize these datasets by finding the most
appropriate inputs suppliers for each duplicated dataset.


How does ``premise`` handle the different granularities between the IAM regions and the Ecoinvent regions?
______________________________________________________________________________________________________

``premise`` simply uses the correspondence between IAM regions and Ecoinvent regions (which are, most of the time
defined by ISO alpha-2 country codes), often provided by the IAM developers.

For example, the REMIND ``REF`` region is associated with the following ecoinvent regions:

- AM
- AZ
- BY
- GE
- KZ
- KG
- MD
- RU
- TJ
- TM
- UA
- UZ

If a technology needs to be included within a market for that region (e.g., coal-based electricity),
``premise`` looks for datasets for that technology (e.g., ``electricity production, hard coal``)
in the ecoinvent database that are located in any
of these above-listed locations, and calculates supply shares based on the
production volumes information provided in each of these datasets (i.e., under the ``production volumes`` field).
Hence, coal-based electricity in the ``REF`` electricity market is supplied
by several coal-based electricity datasets, each of which is located in a different country (see list above)
according to their current production volumes. This approach highlights
a limitation, where current production volumes are used to calculate
supply mix for a given technology within a given IAM region.


Consistency with climate targets
--------------------------------

How do we ensure consistency between IAM scenario and pLCA results (in terms of global warming / temperature increase)?
_______________________________________________________________________________________________________________________

In theory, there is consistency between the IAM scenario and pLCA database
when 100% of the IAM variables and related projections are integrated
into the pLCA database.

This is not the case today, as ``premise`` only integrates a subset of IAM variables, notably those that relate to:

- power production
- steel production
- cement production
- fuel production
- transport

Hence, important sectors are still left out, such as:

- agriculture
- heat
- chemicals
- paper

Also, sectors that are considered by ``premise`` are not fully
or perfectly integrated, as:

- some IAM variables are sometimes not available (e.g., efficiency).
- some IAM variables are sometimes not considered by ``premise`` (e.g., fuel mix for cement production)

Hence, ``premise``-generated databases are not fully consistent with the IAM scenario, including
its climate target. If an ambitious climate target is considered, the use of ``premise``-generated
databases probably leads to an overestimate of GHG emissions, since sectors
that are expected to under mitigation measures are left unchanged. It will however
mostly depend on the product system you analyze.


Additional inventories
----------------------

Can additional inventories be modelled with parameters? If so, how are they used?
_________________________________________________________________________________

Additional inventories (imported as such or via data packages) can be modelled with
(brightway2) parameters, but those will not be considered by ``premise``.

Can some parameters of the additional inventories be made scenario- and time-dependant?
_______________________________________________________________________________________

Yes, via the use of data packages. Data packages allow to package additional scenarios
to be considered in addition to the global IAM scenario. With data packages,
it is possible to map the efficiency of processes to a variable. That variable
can vary over time and across scenarios. Besides efficiency, it is also possible
to change a market mix, distribution losses or any other aspects, of a
product's supply chain, via the use of variables in data packages.

Can ``premise`` manage an efficiency evolution for the additional inventories?
______________________________________________________________________________

Yes, via the use of data packages (see User-defined scenarios section). It is possible to map
the efficiency of processes to a variable. That variable can vary over time and across scenarios.

Efficiency adjustments
----------------------

Is the calculated scaling factor (ratio of efficiencies in year 20XX vs 2020) applied to all inputs of the transformed dataset, or only to the energy feedstock input?
______________________________________________________________________________________________________________________________________________________________________

It depends on the nature of the process. For energy conversion processes (e.g., power generation),
all inputs are scaled up or down. For processes that convert energy and material (e.g., cement or steel production),
only the inputs that relate to energy (e.g., fuel, electricity) inputs are scaled up or down, the input of material
remaining unchanged.

What happens if the IAM does not provide efficiencies for certain processes?
____________________________________________________________________________

They will be ignored and the efficiency of said process wil not be adjusted.

Why use external data sources for PV efficiency, rather than the output of IAM?
_______________________________________________________________________________

Efficiency values for photovoltaic panels are not always provided by IAM scenarios.
When they are, they are often constant (i.e., the efficiency does not increase over time).
This can become an issue when they represent a significant share
of the electricity mix. Hence, at the moment, we use external sources
to document the projected efficiency of photovoltaic modules.
A venue of improvement may be to use IAM efficiency variables for
photovoltaic panels when available, and fall back on external sources if not.