(adjustment_section)=
# Adjusting with Depth and Pressure

The `dry_rock` model can be adjusted according to depth and pressure:

```yaml
dry_rock:
    model:
        # model type is one of:
        #  * polyfit
        #  * friable_sand
        #  * patchy_cement
        type: polyfit
    porosity:
        # porosity of dry rock gotten from column in csv file
        column: porosity
    # Default polyfit coefficients is used which is:
    # density: [[0.0, 0.0], [1.0, -1.0]]
    # bulk_modulus: [[2900.0, -1300.0]]
    # shear_modulus: [[1700.0, -800.0]]
    # In  dry_rock, effective reference pressure is used:
    # pressure = min(max_pressure, overburden_pressure - reference_pressure)
    adjustments:
       - type: pressure_dependency
         # In pressure dependency effective rock pressure is used:
         # pressure = min(max_pressure, overburden_pressure - rock_pressure)
         # With the exception of the powerfit model which uses rock pore pressure:
         # pressure = rock_pressure
         model:
              # model type is one of:
              # * polyfit
              # * expfit
              # * logfit
              # * powerfit
              # * patchy_cement
              # * friable_sand
             type: powerfit
             coefficients:
                # adjusted_bulk_modulus = original_bulk_modulus +
                #         20.0 * delta_pressure^-40.0
                bulk_modulus: [20.0E+9, -46.0]
                # adjusted_shear_modulus = original_shear_modulus +
                #         12.0 * delta_pressure^-20.0
                vp_over_vs: [12.0, -26.0]
                # adjusted_density = original_density
                density: [1.0, 0.0]
       - type: depth_trend
         depth:
           column: depth
         reference_depth: 0.0
         max_depth: 4000.0
         # delta_depth = min(depth - reference_depth, max_depth)
         # coefficients in a 2d polynomial which
         # adjusts dry rock properties to depth
         coefficients:
           # adjusted_bulk_modulus = c[i][j] * delta_depth^i * original_bulk_modulus^j
           bulk_modulus: [[0.0, 1.5], [1.0, 0.0]]
           # adjusted_shear_modulus = c[i][j] * delta_depth^i * original_shear_modulus^j
           shear_modulus: [[0.0, 1.2], [1.0, 0.0]]
           # adjusted_density = c[i][j] * delta_depth^i * original_density^j
           density: [[0.0, 1.35], [1.0, 0.0]]
```

In the above example, the `dry_rock` material changes both with depth and pressure according to
the `depth_trend` and `pressure_dependency` adjustment. Depth trend simply adjusts the material 
according to a two dimensional polynomial (see the above example for details).

## Pressure Dependency

The models for pressure dependency is one of `polyfit`, `expfit`, `logfit`,
`powerfit`, `patchy_cement`, and `friable_sand`. These are used to calculate
three factors: `f_density`, `f_1` and `f_2`.

First the dry rock properties are calculated. Then the pressure depenency at
reference pressure and at effective pressure is calculated and the final
adjust dry rock properties are as follows, except for `powerfit` which we will
discuss last:

```
adjusted_density = density * (fd(effective_pressure) / fd(reference_pressure))
adjusted_bulk_modulus = bulk_modulus * (f1(effective_pressure) / f1(reference_pressure))
adjusted_shear_modulus = shear_modulus * (f2(effective_pressure) / f2(reference_pressure))
```

If velocity coefficients are given to the model then the above calculation is carried
out with `primary_velocity` and `secondary_velocity` instead of `bulk_modulus` and `shear_modulus`. Except
for the `friable_sand` and `patchy_cement` models which always use moduli.

The models `polyfit`, `friable_sand`, and `patchy_cement` are the same models used
in the [`dry_rock` section](dry_rock). The factors calculated from these models
are the `density`, `bulk_modulus` and `shear_modulus` for the corresponding `dry_rock` model, ie.

```
f_d=density
f_1=bulk_modulus
f_2=shear_modulus
```


Given that you have the following coefficient section:

```yaml
 coefficients:
    bulk_modulus: [cb1, cb2, cb3]
    vp_over_vs: [cv1, cv2, cv3]
    density: [cd1, cd2, cd3]
```

The `expfit` model is:

```
f_d = cd1 + cd2*e^(pressure/cd3)
f_1 = cb1 + cb2*e^(pressure/cb3)
f_2 = cv1 + cv2*e^(pressure/cv3)
```

The `logfit` model is:

```
f_d = cd1 + cd2*log_10(pressure)
f_1 = cb1 + cb2*log_10(pressure)
f_2 = cv1 + cv2*log_10(pressure)
```

Finally, `powerfit` calculates adjusted `dry_rock` in a quite different way from the other
models. It uses `vp_over_vs` coefficients:

```
coefficients:
  bulk_modulus: [cb1, cb2]
  vp_over_vs: [cv1, cv2]
  density: [cd1, cd2]
```

Which is used to calculate the following factors:

```
fd = cd1 + pressure ^ cd2
f1 = cb1 + pressure ^ fb2
f2 = cv1 + pressure ^ fv2
```

Then an adjusted `vp_over_vs` is calculated:

```
adjusted_vp_over_vs  = vp_over_vs + f2(effective) - f2(ref)
```

And the final adjusted values are calculated from that:

```
adjusted_density = density + fd
adjusted_bulk_modulus = bulk_modulus + f1(eff) - f(reff)
adjusted_shear_modulus = bulk_modulus / ((vp_over_vs + f2(effective) - f2(ref))^2 - 4/3)
```
