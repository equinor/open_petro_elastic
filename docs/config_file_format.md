## Config File Format

Input to `open_petro_elastic` consists of a yaml
config file plus an optional csv data file.

The yaml config consists of three sections, minerals, fluids and dryrock. Elastic
values for dry rock, having the mineral as its pore-material and saturated with the fluid
is outputed.

### Minerals

Consists only of a list of constituents which will be mixed using Hashin-Shtrikman,
like so:

```yaml
minerals:
    constituents:
        - material: #sand
            bulk_modulus: 36.0E+9 #Pa
            shear_modulus: 44.0E+9 #Pa
            density: 2650.0 # kg/m3
          fraction: 0.5
        - material: # shale
            bulk_modulus: 11.0E+9 #Pa
            shear_modulus: 4.5E+9 #Pa
            density: 1.0 # kg/m3
          # fraction of shale is 1 - sum(other materials)
```

### Fluids

Same as minerals, but a choice of mixing methods is given, and temperature and
pressure can be given for costituents from theoretical models:

```yaml
fluids:
    mix_method: "wood" # either brie or wood
    temperature: 75.0 # celsius
    pressure: 1000.0 # Pa
    constituents:
        - material:
            type: "oil"
            gas_gravity: 0.5 # ratio to air (air gas gravity is 1)
            reference_density: 800.0 # kg/m3
            gas_oil_ratio: 153.0 # ratio of gas to oil
          #fraction set to 1 - sum(other_fluid_fractions)
        - material:
            type: "brine"
            salinity: 45000 # ppm
          fraction: 0.3 # Ratio of brine to remaining fluid
        - material:
            type: "gas"
            gas_gravity: 0.5 # ratio to air (air gas gravity is 1)
          fraction: 0.3 # Ratio of gas to remaining fluid
```

The list can be any number of arbitrary materials (see example in mineral section above) or
given one of the types `["oil", "brine", "gas", "condensate"]`.


### Dry Rock

This section describes the porous dry rock to the corresponding dense mineral.

```yaml
dry_rock:
    model:
        type: "polyfit"
        coefficients:
           density: [[0.0, 0.0], [1.0, -1.0]]
           bulk_modulus: [[2900.0, -1300.0]]
           shear_modulus: [[1700.0, -800.0]]
    porosity: 0.36
    adjustments:
       - model: "powerfit"
         overburden_pressure: 100.0E6
         pressure: 70.0E+6 Pa
         reference_pressure: 30.0E+6 # Pa
         coefficients:
            bulk_modulus: [20.0E+9, -46.0]
            vp_over_vs: [12.0, -26.0]
            density: [1.0, 0.0]
       - depth: 1000.0
         reference_depth: 0.0
         coefficients:
           bulk_modulus: [[0.0, 1.5], [1.0, 0.0]]
           shear_modulus: [[0.0, 1.2], [1.0, 0.0]]
           density: [[0.0, 1.35], [1.0, 0.0]]
```

In the config above dry rock is created using a polynomial fit to
porosity, however, it is also possible to use the theoretical models
`"friable_sand"` and `"patchy_cement"`, see example configs for useage.

The polynomial fit coefficients are such that the terms of the polynoimal,
for e.g. density, are `dry_rock_density = sum_ij mineral_density^i * porosity^j * coeff[i][j]`.
So the polynomials in the example above are

```
dry_rock_density = mineral_density * (1 - porosity)
dry_rock_bulk_modulus = 2900 - 1300.0 * porosity
dry_rock_shear_modulus = 1700 - 800.0 * porosity

```

Adjustments (pressure dependency and depth trend) are applied in the order they
are given in the adjustments list. In the example above, pressure dependency using
the powerfit model is applied first, then depth trend which uses a polynomial fit.
See `examples/example2.yaml` for a detailed description of adjustments.


### Inserting from datafile
Values from the csv file are inserted into the config file as columns of vectors,
`open_petro_elastic` then maps its computation over the columns row-wise. The insertion
into the yaml file of columns is done by giving the column name to the value you
want to insert it into. For instance

```yaml
minerals:
    constituents:
        - material: # shale
            bulk_modulus:
                column: "shale_k"
            shear_modulus:
                column: "shale_s"
            density:
                column: "shale_r"
```

This specifies that the minerals first constituent should have a material whos
`bulk_modulus`, `shear_modulus`, and `density` can be found in the data files
columns with headers `shale_k`, `shale_s`, and `shale_r` respectively. Lets say
the contents of the data file is as follows:

```csv
shale_k,shale_s,shale_r
10E+9,4E+9,2000
8E+9,4.5E+9,3000
9E+9,5E+9,2666
```

Then the output from `open_petro_elastic` will have three rows corresponding to
the three rows in the input.

If a position is given only one value, as opposed to a vector, that same
value is used for all rows. So in our previous example, if the config
looks like this

```yaml
minerals:
    constituents:
        - material: # shale
            bulk_modulus:
                column: "shale_k"
            shear_modulus:
                column: "shale_s"
            density: 1000
```

A density of 1000 kg/m³ is applied in all three output rows.
