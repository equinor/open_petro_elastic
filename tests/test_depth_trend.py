from open_petro_elastic.config import DepthTrend, Pressure


def test_dryrock_depthtrend_vp_vs(
    snapshot,
    max_depth,
    depth,
    reference_depth,
    mineral_mix_material,
    modulifit_dryrock,
    depth_coefficients,
):
    modulifit_dryrock.adjustments = [
        DepthTrend(
            max_depth=max_depth,
            depth=depth,
            reference_depth=reference_depth,
            coefficients=depth_coefficients,
        )
    ]

    snapshot.assert_match(modulifit_dryrock.material(mineral_mix_material, Pressure()))


def test_dryrock_depthtrend_k_my(
    snapshot,
    mineral_mix_material,
    max_depth,
    depth,
    reference_depth,
    modulifit_dryrock,
    depth_coefficients,
):
    modulifit_dryrock.adjustments = [
        DepthTrend(
            depth=depth,
            max_depth=max_depth,
            reference_depth=reference_depth,
            coefficients=depth_coefficients,
        )
    ]

    snapshot.assert_match(modulifit_dryrock.material(mineral_mix_material, Pressure()))
