from open_petro_elastic.material.sandstone import murphy_coordination_number


def test_murphy_coordination_number(snapshot):
    snapshot.assert_match(murphy_coordination_number(0.1))
