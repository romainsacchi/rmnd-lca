Consequential modelling
=======================

The premise module allows users to import and adjust
the consequential system model of the ecoinvent database
v3.8 and 3.9, with a focus on electricity and fuel markets.
This work is based on a publication with a pre-print available
at https://chemrxiv.org/engage/chemrxiv/article-details/63ee10cdfcfb27a31fe227df

Currently, the identification of marginal supplying
technologies is limited
to the electricity and fuel sectors.

Some technologies are excluded from the marginal markets
due to constraints on their feedstock availability.
This typically applies to waste-to-energy (e.g., waste-based CHP)
or waste-to-fuel (e.g., residue-based biofuel) plants.
For steel markets, only the BF-BOF route is considered.


Some imported inventories cannot be
directly linked to the ecoinvent consequential database.
To address this, a mapping file is provided under
https://github.com/polca/premise/blob/master/premise/data/consequential/blacklist.yaml
which proposes alternative candidates to link to the ecoinvent consequential database.


How does it work?
-----------------

From the user viewpoint, the process is as follows:

* prepare a set of parameters that condition the identification of the marginal electricity suppliers
* supply the parameters to `NewDatabase()`
* point to the your local ecoinvent consequential database

The parameters used to identify marginal suppliers that make up
a market are:

* range time
* duration
* foresight
* lead time
* capital replacement rate
* measurement
* weighted slope start
* weighted slope end

Range time
^^^^^^^^^^

Integer. Years. To measure the trend around the point where the additional
capital will be installed, a range of n years before and after the point
is taken as the time interval. Note that if set to a value other than 0,
the duration argument must be set to 0.

Duration
^^^^^^^^

Integer. Years. Duration over which the change in demand should be measured.
Note that if set to a value other than 0, the range time argument must be set to 0.

Foresight
^^^^^^^^^

True or False. In the myopic approach (False), also called a recursive dynamic
approach, the agents have no foresight on relevant parameters (e.g., energy demand,
policy changes and prices) and will only act based on the information they can observe.
In this case, the suppliers can answer to a change in demand only after it has occurred.
In the perfect foresight approach, the future (within the studied time period) is fully
known to all agents. In this case, the decision to invest can be made ahead of the change
in demand. For suppliers with no foresight, capital
will show up a lead time later.

Lead time
^^^^^^^^^

True or False. If False, the market average lead time is taken for all technologies.
If True, technology-specific lead times are used.
If Range and Duration are both set to False, then the lead time is taken as the
time interval (just as with ecoinvent v.3.4).

If you wish to modify the default lead time values used for the different
technologies, you can do so by modifying the file:

https://github.com/polca/premise/blob/master/premise/data/consequential/leadtimes.yaml

Capital replacement rate
^^^^^^^^^^^^^^^^^^^^^^^^

True or False. If False, a horizontal baseline is used.
If True, the capital replacement rate is used as baseline.
The capital replacement rate is equal to -1 divided by
the lifetime (in years) of the technology.

If you wish to modify the default lifetime values used for the different
technologies, you can do so by modifying the file:

https://github.com/polca/premise/blob/master/premise/data/consequential/lifetimes.yaml

Measurement method
^^^^^^^^^^^^^^^^^

0 to 4.

* 0 = slope,
* 1 = linear regression,
* 2 = area under the curve,
* 3 = weighted slope,
* 4 = time interval is split in individual years and measured


Weighted slope start
^^^^^^^^^^^^^^^^^^^^

Weighted slope start is needed for measurement method 3.
The number indicates where the short slope starts
and is given as the fraction of the total time interval.

Weighted slope end
^^^^^^^^^^^^^^^^^^^

Weighted slope end is needed for measurement method 3.
The number indicates where the short slope ends
and is given as the fraction of the total time interval.

Database creation
^^^^^^^^^^^^^^^^^

The user needs to specify the arguments presented above.
If not, the following default arguments value are used:

.. code-block:: python

    args = {
        "range time":0,
        "duration":0,
        "foresight":False,
        "lead time":False,
        "capital replacement rate":False,
        "measurement": 0,
        "weighted slope start": 0.75,
        "weighted slope end": 1.00
    }

.. code-block:: python

    ndb = NewDatabase(
        scenarios = scenarios,
        source_db="ecoinvent 3.8 consequential",
        source_version="3.8",
        key='xxxxxxxxx',
        system_model="consequential",
        system_args=args
    )

    ndb.update_electricity()

    ndb.write_db_to_brightway()