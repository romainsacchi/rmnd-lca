from numbers import Number
import math
import copy


def rescale_exchange(
    exc: dict,
    value: float,
    sector: str,
    remove_uncertainty: bool = True,
    uncertainty_type: int = None,
    loc: float = None,
    scale: float = None,
    minimum: float = None,
    maximum: float = None,
):
    """Function to rescale exchange amount and uncertainty.

    * ``exc`` is an exchange dataset.
    * ``value`` is a number, to be multiplied by the existing amount.
    * ``remove_uncertainty``: Remove (unscaled) uncertainty data, default is ``True``.
    If ``False``, uncertainty data is scaled by the same factor as the amount
    (except for lognormal distributions, where the ``loc`` parameter is scaled by the log of the factor).
    Currently, does not rescale for Bernoulli, Discrete uniform, Weibull, Gamma, Beta, Generalized Extreme value
    and Student T distributions.

    Returns the modified exchange."""
    assert isinstance(exc, dict), "Must pass exchange dictionary"
    assert isinstance(value, Number), "Constant factor ``value`` must be a number"

    # Scale the amount
    exc["amount"] *= value

    flag_change(exchange=exc, factor=float(value), sector=sector)

    # Scale the uncertainty fields if uncertainty is not being removed
    if not remove_uncertainty:
        if not uncertainty_type:
            uncertainty_type = exc.get("uncertainty type", 0)

        # No uncertainty, do nothing
        if uncertainty_type in {0, 6, 7, 8, 9, 10, 11, 12}:
            pass
        elif uncertainty_type in {1, 2, 3, 4, 5}:
            # Scale "loc" by the log of value for lognormal distribution
            if loc:
                exc["loc"] = loc
            else:
                if "loc" in exc and uncertainty_type == 2:
                    exc["loc"] += math.log(value)
                elif "loc" in exc:
                    exc["loc"] *= value

            # "scale" stays the same for lognormal
            # For other distributions, scale "scale" by the absolute value
            if scale:
                exc["scale"] = scale
            else:
                if "scale" in exc and uncertainty_type not in {2}:
                    exc["scale"] *= abs(value)

            # Scale "minimum" and "maximum" by value
            if minimum and maximum:
                exc["minimum"] = minimum
                exc["maximum"] = maximum
            else:
                for bound in ("minimum", "maximum"):
                    if bound in exc:
                        exc[bound] *= value

    # If remove_uncertainty is True, then remove all uncertainty info
    elif remove_uncertainty:
        FIELDS = (
            "scale",
            "minimum",
            "maximum",
        )
        exc["uncertainty type"] = 0
        exc["loc"] = exc["amount"]
        for field in FIELDS:
            if field in exc:
                del exc[field]

    return exc


def flag_change(
    exchange: dict, sector: str, value: float = None, factor: float = None
) -> dict:
    """
    Flag an exchange with the factor applied on the amount.
    """

    if "original value" not in exchange:
        exchange["original value"] = copy.deepcopy(exchange["amount"])

    if value:
        exchange[sector] = value
        return exchange

    if factor:
        exchange[sector] = exchange["original value"] * factor
        return exchange

    raise ValueError("`value` or `factor` must be passed.")
