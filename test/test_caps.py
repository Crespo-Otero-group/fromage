import pytest

def capital_case(x):
    if not isinstance(x, str):
        raise TypeError('Gimme string pls')
    return x.capitalize()

def test_caps_case():
    assert capital_case('bleh') == 'Bleh'

def test_with_exception():
    with pytest.raises(TypeError):
        capital_case(9)
