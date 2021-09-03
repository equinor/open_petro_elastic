import pytest
from numpy.testing import assert_allclose, assert_array_equal
import numpy as np
from open_petro_elastic.material.span_wagner.carbon_dioxide import carbon_dioxide_density, carbon_dioxide_pressure,\
    carbon_dioxide, array_carbon_dioxide_density, carbon_dioxide_bulk_modulus


# Excerpt from Table 34 of Span & Wagner [2]
table_34_header = ["temperature", "pressure", "density", "vapor", "speed_of_sound"]
table_34_data = np.array([
    [216.592, 0.51796, 1178.46, False, 975.85],
    [216.592, 0.51796, 13.761, True, 222.78],
    [218, 0.55042, 1173.40, False, 965.66],
    [218, 0.55042, 14.584, True, 222.94],
    [220, 0.59913, 1166.14, False, 951.21],
    [220, 0.59913, 15.817, True, 223.15],
    [222, 0.65102, 1158.81, False, 936.79],
    [222, 0.65102, 17.131, True, 223.31],
    [224, 0.70621, 1151.40, False, 922.37],
    [224, 0.70621, 18.530, True, 223.44],
    [226, 0.76484, 1143.92, False, 907.95],
    [226, 0.76484, 20.016, True, 223.52],
    [228, 0.82703, 1136.34, False, 893.53],
    [228, 0.82703, 21.595, True, 223.57],
    [230, 0.89291, 1128.68, False, 879.09],
    [230, 0.89291, 23.271, True, 223.57],
    [232, 0.96262, 1120.93, False, 864.63],
    [232, 0.96262, 25.050, True, 223.54],
    [234, 1.0363, 1113.08, False, 850.14],
    [234, 1.0363, 26.936, True, 223.46],
    [236, 1.1141, 1105.12, False, 835.61],
    [236, 1.1141, 28.935, True, 223.33],
    [238, 1.1961, 1097.05, False, 821.02],
    [238, 1.1961, 31.052, True, 223.17],
    [240, 1.2825, 1088.87, False, 806.38],
    [240, 1.2825, 33.295, True, 222.96],
    [242, 1.3734, 1080.56, False, 791.67],
    [242, 1.3734, 35.670, True, 222.70],
    [244, 1.4690, 1072.13, False, 776.87],
    [244, 1.4690, 38.184, True, 222.40],
    [246, 1.5693, 1063.56, False, 761.97],
    [246, 1.5693, 40.845, True, 222.06],
    [248, 1.6746, 1054.84, False, 746.95],
    [248, 1.6746, 43.662, True, 221.66],
    [250, 1.7850, 1045.97, False, 731.78],
    [250, 1.7850, 46.644, True, 221.22],
    [252, 1.9007, 1036.93, False, 716.44],
    [252, 1.9007, 49.801, True, 220.72],
    [254, 2.0217, 1027.72, False, 700.88],
    [254, 2.0217, 53.144, True, 220.17],
    [256, 2.1483, 1018.32, False, 685.08],
    [256, 2.1483, 56.685, True, 219.56],
    [258, 2.2806, 1008.71, False, 668.99],
    [258, 2.2806, 60.438, True, 218.90],
    [260, 2.4188, 998.89, False, 652.58],
    [260, 2.4188, 64.417, True, 218.19],
    [262, 2.5630, 988.83, False, 635.84],
    [262, 2.5630, 68.640, True, 217.41],
    [264, 2.7134, 978.51, False, 618.75],
    [264, 2.7134, 73.124, True, 216.59],
    [266, 2.8701, 967.92, False, 601.31],
    [266, 2.8701, 77.891, True, 215.70],
    [268, 3.0334, 957.04, False, 583.54],
    [268, 3.0334, 82.965, True, 214.76],
    [270, 3.2033, 945.83, False, 565.46],
    [270, 3.2033, 88.374, True, 213.75],
    [272, 3.3802, 934.26, False, 547.11],
    [272, 3.3802, 94.140, True, 212.68],
    [274, 3.5642, 922.30, False, 528.51],
    [274, 3.5642, 100.32, True, 211.55],
    [276, 3.7555, 909.90, False, 509.71],
    [276, 3.7555, 106.95, True, 210.35],
    [278, 3.9542, 897.02, False, 490.72],
    [278, 3.9542, 114.07, True, 209.07],
    [280, 4.1607, 883.58, False, 471.54],
    [280, 4.1607, 121.74, True, 207.72],
    [282, 4.3752, 869.52, False, 452.19],
    [282, 4.3752, 130.05, True, 206.28],
    [284, 4.5978, 854.74, False, 432.63],
    [284, 4.5978, 139.09, True, 204.74],
    [286, 4.8289, 839.12, False, 412.81],
    [286, 4.8289, 148.98, True, 203.10],
    [288, 5.0688, 822.50, False, 392.63],
    [288, 5.0688, 159.87, True, 201.34],
    [290, 5.3177, 804.67, False, 371.95],
    [290, 5.3177, 171.96, True, 199.45],
    [292, 5.5761, 785.33, False, 350.49],
    [292, 5.5761, 185.55, True, 197.38],
    [294, 5.8443, 764.09, False, 327.85],
    [294, 5.8443, 201.06, True, 195.09],
    [296, 6.1227, 740.28, False, 303.44],
    [296, 6.1227, 219.14, True, 192.49],
    [298, 6.4121, 712.77, False, 276.42],
    [298, 6.4121, 240.90, True, 189.38],
    [300, 6.7131, 679.24, False, 245.67],
    [300, 6.7131, 268.58, True, 185.33],
    [301, 6.8683, 658.69, False, 228.18],
    [301, 6.8683, 286.15, True, 182.61],
    [302, 7.0268, 633.69, False, 208.08],
    [302, 7.0268, 308.15, True, 178.91],
    [303, 7.1890, 599.86, False, 182.14],
    [303, 7.1890, 339.00, True, 172.71],
    [304, 7.3555, 530.30, False, 134.14],
    [304, 7.3555, 406.42, True, 147.62],
    [304.1282, 7.3773, 467.60, True, np.nan],  # Speed of sound not provided
], dtype='object')


# Values are extracted from Table 35 of Span & Wagner [2]
table_35_header = ["temperature", "pressure", "density"]
table_35_data = np.array([
    # 0.05 MPa
    [190, 0.05, 1.4089],
    [200, 0.05, 1.3359],
    [210, 0.05, 1.2704],
    [220, 0.05, 1.2112],
    [230, 0.05, 1.1575],
    [240, 0.05, 1.1084],
    [250, 0.05, 1.0634],
    [260, 0.05, 1.0219],
    [270, 0.05, 0.98360],
    [280, 0.05, 0.94810],
    [290, 0.05, 0.91510],
    [300, 0.05, 0.88434],
    [325, 0.05, 0.81585],
    [350, 0.05, 0.75726],
    [375, 0.05, 0.70656],
    [400, 0.05, 0.66224],
    [425, 0.05, 0.62317],
    [450, 0.05, 0.58847],
    [475, 0.05, 0.55743],
    [500, 0.05, 0.52951],
    [525, 0.05, 0.50426],
    [550, 0.05, 0.48130],
    [575, 0.05, 0.46035],
    [600, 0.05, 0.44115],
    [700, 0.05, 0.37809],
    [800, 0.05, 0.33081],
    [900, 0.05, 0.29404],
    [1000, 0.05, 0.26463],
    [1100, 0.05, 0.24057],
    # 0.1 MPa
    [200, 0.1, 2.6980],
    [210, 0.1, 2.5617],
    [220, 0.1, 2.4394],
    [230, 0.1, 2.3288],
    [240, 0.1, 2.2282],
    [250, 0.1, 2.1363],
    [260, 0.1, 2.0519],
    [270, 0.1, 1.9741],
    [280, 0.1, 1.9021],
    [290, 0.1, 1.8352],
    [300, 0.1, 1.7730],
    [325, 0.1, 1.6348],
    [350, 0.1, 1.5167],
    [375, 0.1, 1.4147],
    [400, 0.1, 1.3257],
    [425, 0.1, 1.2472],
    [450, 0.1, 1.1776],
    [475, 0.1, 1.1154],
    [500, 0.1, 1.0594],
    [525, 0.1, 1.0088],
    [550, 0.1, 0.96283],
    [575, 0.1, 0.92087],
    [600, 0.1, 0.88242],
    [700, 0.1, 0.75619],
    [800, 0.1, 0.66158],
    [900, 0.1, 0.58803],
    [1000, 0.1, 0.52921],
    [1100, 0.1, 0.48109],
    # 1.00 MPa
    [220, 1, 1167.03],
    [225, 1, 1148.32],
    [230, 1, 1128.97],
    [235, 1, 25.665],
    [240, 1, 24.857],
    [245, 1, 24.117],
    [250, 1, 23.435],
    [255, 1, 22.803],
    [260, 1, 22.215],
    [265, 1, 21.664],
    [270, 1, 21.147],
    [275, 1, 20.660],
    [280, 1, 20.199],
    [285, 1, 19.763],
    [290, 1, 19.349],
    [295, 1, 18.955],
    [300, 1, 18.579],
    [305, 1, 18.221],
    [310, 1, 17.878],
    [315, 1, 17.549],
    [320, 1, 17.234],
    [325, 1, 16.932],
    [330, 1, 16.641],
    [335, 1, 16.361],
    [340, 1, 16.092],
    [345, 1, 15.832],
    [350, 1, 15.581],
    [360, 1, 15.105],
    [370, 1, 14.659],
    [380, 1, 14.241],
    [390, 1, 13.848],
    [400, 1, 13.477],
    [410, 1, 13.127],
    [420, 1, 12.796],
    [430, 1, 12.482],
    [440, 1, 12.183],
    [450, 1, 11.899],
    [460, 1, 11.629],
    [470, 1, 11.371],
    [480, 1, 11.125],
    [490, 1, 10.889],
    [500, 1, 10.664],
    [525, 1, 10.141],
    [550, 1, 9.6675],
    [575, 1, 9.2375],
    [600, 1, 8.8449],
    [625, 1, 8.4849],
    [650, 1, 8.1535],
    [675, 1, 7.8474],
    [700, 1, 7.5638],
    [800, 1, 6.6102],
    [900, 1, 5.8718],
    [1000, 1, 5.2826],
    [1100, 1, 4.8014],
    # 10 MPa
    [220, 10, 1185.63],
    [225, 10, 1168.59],
    [230, 10, 1151.15],
    [235, 10, 1133.28],
    [240, 10, 1114.92],
    [245, 10, 1095.99],
    [250, 10, 1076.42],
    [255, 10, 1056.11],
    [260, 10, 1034.95],
    [265, 10, 1012.80],
    [270, 10, 989.46],
    [280, 10, 938.22],
    [285, 10, 909.56],
    [290, 10, 878.06],
    [295, 10, 842.67],
    [300, 10, 801.62],
    [305, 10, 751.67],
    [310, 10, 685.77],
    [315, 10, 586.02],
    [320, 10, 448.28],
    [325, 10, 358.04],
    [330, 10, 310.25],
    [335, 10, 280.11],
    [340, 10, 258.62],
    [345, 10, 242.11],
    [350, 10, 228.80],
    [360, 10, 208.25],
    [370, 10, 192.74],
    [380, 10, 180.38],
    [390, 10, 170.18],
    [400, 10, 161.53],
    [410, 10, 154.05],
    [420, 10, 147.48],
    [430, 10, 141.63],
    [440, 10, 136.39],
    [450, 10, 131.63],
    [460, 10, 127.30],
    [470, 10, 123.32],
    [480, 10, 119.65],
    [490, 10, 116.24],
    [500, 10, 113.07],
    [525, 10, 106.01],
    [550, 10, 99.930],
    [575, 10, 94.626],
    [600, 10, 89.941],
    [625, 10, 85.761],
    [650, 10, 81.999],
    [675, 10, 78.592],
    [700, 10, 75.486],
    [800, 10, 65.349],
    [900, 10, 57.759],
    [1000, 10, 51.825],
    [1100, 10, 47.040],
    # 100 MPa
    [240, 100, 1257.21],
    [245, 100, 1246.37],
    [250, 100, 1235.57],
    [255, 100, 1224.80],
    [260, 100, 1214.06],
    [265, 100, 1203.36],
    [270, 100, 1192.69],
    [275, 100, 1182.06],
    [280, 100, 1171.46],
    [285, 100, 1160.90],
    [290, 100, 1150.37],
    [295, 100, 1139.89],
    [300, 100, 1129.45],
    [305, 100, 1119.04],
    [310, 100, 1108.69],
    [315, 100, 1098.38],
    [320, 100, 1088.12],
    [325, 100, 1077.91],
    [335, 100, 1057.66],
    [340, 100, 1047.62],
    [345, 100, 1037.64],
    [350, 100, 1027.73],
    [360, 100, 1008.12],
    [370, 100, 988.80],
    [380, 100, 969.80],
    [390, 100, 951.13],
    [400, 100, 932.81],
    [410, 100, 914.87],
    [420, 100, 897.31],
    [430, 100, 880.14],
    [440, 100, 863.37],
    [450, 100, 847.00],
    [460, 100, 831.05],
    [470, 100, 815.51],
    [480, 100, 800.39],
    [490, 100, 785.67],
    [500, 100, 771.37],
    [525, 100, 737.38],
    [550, 100, 705.86],
    [575, 100, 676.71],
    [600, 100, 649.77],
    [625, 100, 624.90],
    # 800 MPa
    [330, 800, 1493.71],
    [335, 800, 1489.46],
    [340, 800, 1485.25],
    [345, 800, 1481.08],
    [350, 800, 1476.96],
    [360, 800, 1468.83],
    [370, 800, 1460.85],
    [380, 800, 1453.03],
    [390, 800, 1445.35],
    [400, 800, 1437.82],
    [410, 800, 1430.42],
    [420, 800, 1423.15],
    [430, 800, 1416.00],
    [440, 800, 1408.98],
    [450, 800, 1402.08],
    [460, 800, 1395.29],
    [470, 800, 1388.60],
    [480, 800, 1382.03],
    [490, 800, 1375.55],
    [500, 800, 1369.17],
    [525, 800, 1353.64],
    [550, 800, 1338.65],
    [575, 800, 1324.16],
    [600, 800, 1310.14],
    [625, 800, 1296.54],
    [650, 800, 1283.34],
    [675, 800, 1270.52],
    [700, 800, 1258.04],
    [800, 800, 1211.18],
    [900, 800, 1168.47],
    [1000, 800, 1129.16],
    [1100, 800, 1092.77],
])


@pytest.mark.parametrize(table_35_header, table_35_data)
def test_gas_density_and_pressure(temperature, pressure, density):
    """
    Tests that values calculated match that of table 35 from Span & Wagner.
    """
    assert carbon_dioxide_pressure(temperature, density) == pytest.approx(
        pressure, rel=0.002
    )
    assert carbon_dioxide_density(temperature, pressure) == pytest.approx(
        density, rel=0.002
    )


def test_vectorized_gas_density_and_pressure():
    """
    Same as above, except for the vectorized version of the function
    """
    data = table_35_data
    calc_pressure = carbon_dioxide_pressure(data[:, 0], data[:, 2])
    assert_allclose(calc_pressure, data[:, 1], rtol=0.002)
    calc_density = carbon_dioxide_density(data[:, 0], data[:, 1])
    assert_allclose(calc_density, data[:, 2], rtol=0.002)


def test_newton_carbon_dioxide_density():
    """
    Tests the alternative vectorized variant for density
    """
    data = table_35_data
    calc_density = array_carbon_dioxide_density(data[:, 0], data[:, 1], 'auto')
    assert_allclose(calc_density, data[:, 2], rtol=0.002)


@pytest.mark.parametrize(
    table_34_header[:4],
    np.vstack((
        table_34_data[:, :4],
        np.array([
            # From Table 35 of [2]
            [186.436, 0.05, 1.4370, True],
            [194.525, 0.1, 2.7796, True],
            [216.695, 1, 1179.10, False],
            [233.028, 1, 1116.90, False],
            [233.028, 1, 26.006, True],
            [218.600, 10, 1190.34, False],
            [236.031, 100, 1265.83, False],
            [327.673, 800, 1495.70, False],
        ], dtype='object')
    ))
)
def test_density_close_to_phase_boundaries(temperature, pressure, density, vapor):
    calc_pressure = carbon_dioxide_pressure(temperature, density)
    assert calc_pressure == pytest.approx(pressure, rel=0.005)
    calc_density = carbon_dioxide_density(temperature, pressure, force_vapor=vapor)
    assert calc_density == pytest.approx(density, rel=0.005)


def test_carbon_dioxide_pressure_derivative():
    t, p, r = table_35_data.T
    eps = 0.01
    cd_estimate = carbon_dioxide_pressure(t, r + eps) - carbon_dioxide_pressure(t, r - eps)
    cd_estimate /= (2 * eps)
    calculated = carbon_dioxide_pressure(t, r, d_density=1)
    assert_allclose(calculated, cd_estimate, rtol=0.01)


def test_carbon_dioxide_primary_velocity():
    t, r, s = table_34_data[:-1, [0, 2, 4]].astype(np.float).T
    co2 = carbon_dioxide(t, None, r)
    calc_s = co2.primary_velocity
    assert_allclose(s, calc_s, rtol=0.02)


def test_bulk_modulus_array_shapes():
    b0 = carbon_dioxide_bulk_modulus(280, 122)
    assert b0.ndim == 0
    b1 = carbon_dioxide_bulk_modulus(np.array([280]), np.array([122]))
    assert b1.ndim == 1
    assert b0 == b1[0]
    b2 = carbon_dioxide_bulk_modulus(np.array([280, 280]), np.array([122, 122]))
    assert b2.ndim == 1
    assert b2.size == 2
    assert b0 == b2[0]
    assert b0 == b2[1]


def test_interpolate_density():
    t, p, d = table_35_data.T
    outside_bounds = (t < 274) | (t > 473) | (p < 0.5) | (p > 99.5)
    di = carbon_dioxide_density(t, p, interpolate=True)
    assert_array_equal(outside_bounds, np.isnan(di))
    assert_allclose(di[~outside_bounds], d[~outside_bounds], rtol=0.0001)


def test_partially_vectorized_bulk_modulus():
    t = np.array([348, 350])
    p = np.array([29.24, 29.27])
    vbm1 = carbon_dioxide(t[0], p, None).bulk_modulus
    nbm1 = [carbon_dioxide(t[0], _p, None).bulk_modulus for _p in p]
    assert_allclose(vbm1, nbm1)
    vbm2 = carbon_dioxide(t, p[0], None).bulk_modulus
    nbm2 = [carbon_dioxide(_t, p[0], None).bulk_modulus for _t in t]
    assert_allclose(vbm2, nbm2)
