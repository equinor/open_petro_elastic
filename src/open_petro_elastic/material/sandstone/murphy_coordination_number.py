def murphy_coordination_number(porosity):
    """
    Estimation of coordination number for sandstone from porosity by
    Murphy, W.  F. (1982) "Effects of microstructure and pore fluids on the
    acoustic properties of granular sedimentary materials".
    See also:
    * https://www.sciencedirect.com/topics/earth-and-planetary-sciences/coordination-number
    """
    # sciencedirect uses the following fit: 24.0 * exp(-2.547 * critical_porosity) - 0.3731
    # The polynomial below is fitted from values given in the rock physics handbook.
    return 25.98805 * porosity ** 2 - 43.7622 * porosity + 21.6719
